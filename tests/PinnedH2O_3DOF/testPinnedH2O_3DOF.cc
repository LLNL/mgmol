// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Control.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "MGmol.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"
#include "mgmol_run.h"

#ifdef MGMOL_HAS_LIBROM
#include "librom.h"

#include <cassert>
#include <iostream>
#include <time.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstring>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

double calculate_bondlength(const double atom1[3], const double atom2[3])
{
    return sqrt(pow(atom1[0] - atom2[0], 2) + pow(atom1[1] - atom2[1], 2) + pow(atom1[2] - atom2[2], 2));
}

double calculate_bondangle(const double atom1[3], const double atom2[3], const double atom3[3])
{
    double vector1[3] = {atom1[0] - atom2[0], atom1[1] - atom2[1], atom1[2] - atom2[2]};
    double vector2[3] = {atom3[0] - atom2[0], atom3[1] - atom2[1], atom3[2] - atom2[2]};

    double dot_product = vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2];
    double magnitude_product = sqrt(pow(vector1[0], 2) + pow(vector1[1], 2) + pow(vector1[2], 2)) *
                               sqrt(pow(vector2[0], 2) + pow(vector2[1], 2) + pow(vector2[2], 2));
    double angle = acos(dot_product / magnitude_product);

    return angle;
}

void rotation_matrix(const double axis[3], double angle, double matrix[3][3])
{
    double cos_theta = cos(angle);
    double sin_theta = sin(angle);
    double ux = axis[0], uy = axis[1], uz = axis[2];

    matrix[0][0] = cos_theta + ux * ux * (1 - cos_theta);
    matrix[0][1] = ux * uy * (1 - cos_theta) - uz * sin_theta;
    matrix[0][2] = ux * uz * (1 - cos_theta) + uy * sin_theta;

    matrix[1][0] = uy * ux * (1 - cos_theta) + uz * sin_theta;
    matrix[1][1] = cos_theta + uy * uy * (1 - cos_theta);
    matrix[1][2] = uy * uz * (1 - cos_theta) - ux * sin_theta;

    matrix[2][0] = uz * ux * (1 - cos_theta) - uy * sin_theta;
    matrix[2][1] = uz * uy * (1 - cos_theta) + ux * sin_theta;
    matrix[2][2] = cos_theta + uz * uz * (1 - cos_theta);
}

void normalize(double vec[3])
{
    double norm = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    vec[0] /= norm;
    vec[1] /= norm;
    vec[2] /= norm;
}

void cross(const double a[3], const double b[3], double result[3])
{
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

void apply_rotation(const double matrix[3][3], const double vec[3], double result[3])
{
    result[0] = matrix[0][0] * vec[0] + matrix[0][1] * vec[1] + matrix[0][2] * vec[2];
    result[1] = matrix[1][0] * vec[0] + matrix[1][1] * vec[1] + matrix[1][2] * vec[2];
    result[2] = matrix[2][0] * vec[0] + matrix[2][1] * vec[1] + matrix[2][2] * vec[2];
}

void apply_transpose_rotation(const double matrix[3][3], const double vec[3], double result[3])
{
    result[0] = matrix[0][0] * vec[0] + matrix[1][0] * vec[1] + matrix[2][0] * vec[2];
    result[1] = matrix[0][1] * vec[0] + matrix[1][1] * vec[1] + matrix[2][1] * vec[2];
    result[2] = matrix[0][2] * vec[0] + matrix[1][2] * vec[1] + matrix[2][2] * vec[2];
}

void rotate_PinnedH2O(std::vector<double>& positions, std::vector<short>& anumbers,
                      double out_of_plane_rotation_matrix[3][3], double& planar_rotation_angle, bool& flipped_bond)
{
    int O1_idx = -1;
    for (int i = 0; i < 3; i++)
    {
        if (positions[3*i] == 0.0 && positions[3*i+1] == 0.0 && positions[3*i+2] == 0.0)
        {
            O1_idx = i;
            break;
        }
    }
    if (O1_idx == -1) return;

    int H1_idx = (O1_idx + 1) % 3;
    int H2_idx = (O1_idx + 2) % 3;

    double O1[3] = {positions[3*O1_idx], positions[3*O1_idx+1], positions[3*O1_idx+2]};
    double H1[3] = {positions[3*H1_idx], positions[3*H1_idx+1], positions[3*H1_idx+2]};
    double H2[3] = {positions[3*H2_idx], positions[3*H2_idx+1], positions[3*H2_idx+2]};

    double bondlength1 = calculate_bondlength(H1, O1);
    double bondlength2 = calculate_bondlength(H2, O1);
    double bondangle = calculate_bondangle(H1, O1, H2);

    double H1_temp[3], H2_temp[3];
    double H1_rotated[3], H2_rotated[3];

    double plane_normal[3];
    cross(H2, H1, plane_normal);
    normalize(plane_normal);
    double target_plane_normal[3] = {0, 0, 1};
    double dot_product = plane_normal[0] * target_plane_normal[0] +
                         plane_normal[1] * target_plane_normal[1] +
                         plane_normal[2] * target_plane_normal[2];
    double angle_to_align = std::acos(std::min(std::max(dot_product, -1.0), 1.0));
    double axis_to_align[3];
    if (abs(dot_product) > 1.0 - 1e-8)
    {
        axis_to_align[0] = 1.0; 
        axis_to_align[1] = 0.0;
        axis_to_align[2] = 0.0;
    }
    else
    {
        cross(plane_normal, target_plane_normal, axis_to_align);
        normalize(axis_to_align);
    }
    rotation_matrix(axis_to_align, angle_to_align, out_of_plane_rotation_matrix);
    apply_rotation(out_of_plane_rotation_matrix, H1, H1_temp);
    apply_rotation(out_of_plane_rotation_matrix, H2, H2_temp);

    double theta1 = std::atan2(H1_temp[1], H1_temp[0]);
    planar_rotation_angle = -theta1 + bondangle / 2.0;
    double planar_rotation_matrix[3][3];
    rotation_matrix(target_plane_normal, planar_rotation_angle, planar_rotation_matrix);
    apply_rotation(planar_rotation_matrix, H1_temp, H1_rotated);
    apply_rotation(planar_rotation_matrix, H2_temp, H2_rotated);

    flipped_bond = (bondlength1 < bondlength2);
    if (flipped_bond)
    {
        H1_rotated[1] *= -1.0;
        H2_rotated[1] *= -1.0;
        std::swap(H1_rotated, H2_rotated);
    }

    positions[0] = H2_rotated[0];
    positions[1] = H2_rotated[1];
    positions[2] = H2_rotated[2];
    positions[3] = 0.0;
    positions[4] = 0.0;
    positions[5] = 0.0;
    positions[6] = H1_rotated[0];
    positions[7] = H1_rotated[1];
    positions[8] = H1_rotated[2];

    anumbers[0] = 1;
    anumbers[1] = 8;
    anumbers[2] = 1;
}

void transpose_rotate_PinnedH2O(std::vector<double>& positions, std::vector<double>& forces,
                                const double out_of_plane_rotation_matrix[3][3], 
                                const double planar_rotation_angle, const bool flipped_bond)
{
    double H2_rotated[3] = {positions[0], positions[1], positions[2]};
    double O1_rotated[3] = {positions[3], positions[4], positions[5]};
    double H1_rotated[3] = {positions[6], positions[7], positions[8]};

    double f_H2_rotated[3] = {forces[0], forces[1], forces[2]};
    double f_O1_rotated[3] = {forces[3], forces[4], forces[5]};
    double f_H1_rotated[3] = {forces[6], forces[7], forces[8]};

    if (flipped_bond)
    {
        H1_rotated[1] *= -1.0;
        H2_rotated[1] *= -1.0;
        f_O1_rotated[1] *= -1.0;
        f_H1_rotated[1] *= -1.0;
        f_H2_rotated[1] *= -1.0;
    }

    double planar_rotation_matrix[3][3];
    double target_plane_normal[3] = {0, 0, 1};
    rotation_matrix(target_plane_normal, planar_rotation_angle, planar_rotation_matrix);

    double H1_temp[3], H2_temp[3];
    apply_transpose_rotation(planar_rotation_matrix, H1_rotated, H1_temp);
    apply_transpose_rotation(planar_rotation_matrix, H2_rotated, H2_temp);

    double f_O1_temp[3], f_H1_temp[3], f_H2_temp[3];
    apply_transpose_rotation(planar_rotation_matrix, f_O1_rotated, f_O1_temp);
    apply_transpose_rotation(planar_rotation_matrix, f_H1_rotated, f_H1_temp);
    apply_transpose_rotation(planar_rotation_matrix, f_H2_rotated, f_H2_temp);

    double H1_restored[3], H2_restored[3];
    apply_transpose_rotation(out_of_plane_rotation_matrix, H1_temp, H1_restored);
    apply_transpose_rotation(out_of_plane_rotation_matrix, H2_temp, H2_restored);

    double f_O1_restored[3], f_H1_restored[3], f_H2_restored[3];
    apply_transpose_rotation(out_of_plane_rotation_matrix, f_O1_temp, f_O1_restored);
    apply_transpose_rotation(out_of_plane_rotation_matrix, f_H1_temp, f_H1_restored);
    apply_transpose_rotation(out_of_plane_rotation_matrix, f_H2_temp, f_H2_restored);

    positions[0] = H2_restored[0];
    positions[1] = H2_restored[1];
    positions[2] = H2_restored[2];
    positions[6] = H1_restored[0];
    positions[7] = H1_restored[1];
    positions[8] = H1_restored[2];

    forces[0] = f_H2_restored[0];
    forces[1] = f_H2_restored[1];
    forces[2] = f_H2_restored[2];
    forces[3] = f_O1_restored[0];
    forces[4] = f_O1_restored[1];
    forces[5] = f_O1_restored[2];
    forces[6] = f_H1_restored[0];
    forces[7] = f_H1_restored[1];
    forces[8] = f_H1_restored[2];
}

int main(int argc, char** argv)
{
    int mpirc = MPI_Init(&argc, &argv);
    if (mpirc != MPI_SUCCESS)
    {
        std::cerr << "MPI Initialization failed!!!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    MPI_Comm comm = MPI_COMM_WORLD;

    /*
     * Initialize general things, like magma, openmp, IO, ...
     */
    mgmol_init(comm);

    /*
     * read runtime parameters
     */
    std::string input_filename("");
    std::string lrs_filename;
    std::string constraints_filename("");

    float total_spin = 0.;
    bool with_spin   = false;

    po::variables_map vm;

    // read from PE0 only
    if (MPIdata::onpe0)
    {
        read_config(argc, argv, vm, input_filename, lrs_filename,
            constraints_filename, total_spin, with_spin);
    }

    MGmol_MPI::setup(comm, std::cout, with_spin);
    MGmol_MPI& mmpi      = *(MGmol_MPI::instance());
    MPI_Comm global_comm = mmpi.commGlobal();

    /*
     * Setup control struct with run time parameters
     */
    Control::setup(global_comm, with_spin, total_spin);
    Control& ct = *(Control::instance());

    ct.setOptions(vm);

    int ret = ct.checkOptions();
    if (ret < 0) return ret;

    mmpi.bcastGlobal(input_filename);
    mmpi.bcastGlobal(lrs_filename);

    // Enter main scope
    {
        if (MPIdata::onpe0)
        {
            std::cout << "-------------------------" << std::endl;
            std::cout << "Construct MGmol object..." << std::endl;
            std::cout << "-------------------------" << std::endl;
        }

        MGmolInterface* mgmol = new MGmol<ExtendedGridOrbitals>(global_comm,
            *MPIdata::sout, input_filename, lrs_filename, constraints_filename);

        if (MPIdata::onpe0)
        {
            std::cout << "-------------------------" << std::endl;
            std::cout << "MGmol setup..." << std::endl;
            std::cout << "-------------------------" << std::endl;
        }
        mgmol->setup();

        if (MPIdata::onpe0)
        {
            std::cout << "-------------------------" << std::endl;
            std::cout << "Setup done..." << std::endl;
            std::cout << "-------------------------" << std::endl;
        }

        // here we just use the atomic positions read in and used
        // to initialize MGmol
        std::vector<double> positions;
        mgmol->getAtomicPositions(positions);
        std::vector<short> anumbers;
        mgmol->getAtomicNumbers(anumbers);
        if (MPIdata::onpe0)
        {
            std::cout << "Positions:" << std::endl;
            std::vector<short>::iterator ita = anumbers.begin();
            for (std::vector<double>::iterator it = positions.begin();
                 it != positions.end(); it += 3)
            {
                std::cout << *ita;
                for (int i = 0; i < 3; i++)
                    std::cout << "    " << *(it + i);
                std::cout << std::endl;
                ita++;
            }
        }

        // rotate the molecule to the reference coordinate system
        int atom_order[3];
        double out_of_plane_rotation_matrix[3][3];
        double planar_rotation_angle;
        bool flipped_bond;
        rotate_PinnedH2O(positions, anumbers, out_of_plane_rotation_matrix, planar_rotation_angle, flipped_bond);

        Mesh* mymesh             = Mesh::instance();
        const pb::Grid& mygrid   = mymesh->grid();
        const pb::PEenv& myPEenv = mymesh->peenv();

        // compute energy and forces again with projected problem onto ROM subspace
        const int rdim = ct.getROMOptions().num_orbbasis;
        if (rdim != ct.numst)
        {
            std::cerr << "The number of functions in the ROM basis file, "
                      << rdim << " is not equal to ct.numst, " << ct.numst
                      << std::endl;
            MPI_Abort(mmpi.commSameSpin(), 0);
        }

        std::shared_ptr<ProjectedMatricesInterface> projmatrices
            = mgmol->getProjectedMatrices();

        ExtendedGridOrbitals orbitals("new_orbitals", mygrid, mymesh->subdivx(),
            ct.numst, ct.bcWF, projmatrices.get(), nullptr, nullptr, nullptr,
            nullptr);

        orbitals.set(ct.getROMOptions().basis_file, ct.numst); 
        orbitals.orthonormalizeLoewdin();
        orbitals.setDataWithGhosts(true);

        // set the iterative index to 1 to differentiate it from first instance
        // in MGmol initial() function. This is not very clean and could be
        // better designed, but works for now
        orbitals.setIterativeIndex(10);

        // set initial DM with uniform occupations
        projmatrices->setDMuniform(ct.getNelSpin(), 0);
        projmatrices->printDM(std::cout);

        //
        // evaluate energy and forces with ROM bases just read
        //
        std::vector<double> forces;
        double eks = mgmol->evaluateDMandEnergyAndForces(
            &orbitals, positions, anumbers, forces);

        // print out results
        if (MPIdata::onpe0)
        {
            std::cout << "Eks: " << eks << std::endl;
            std::cout << "Positions in the reference domain:" << std::endl;
            std::vector<short>::iterator ita = anumbers.begin();
            for (std::vector<double>::iterator it = positions.begin();
                 it != positions.end(); it += 3)
            {
                std::cout << *ita;
                for (int i = 0; i < 3; i++)
                    std::cout << "    " << *(it + i);
                std::cout << std::endl;
                ita++;
            }

            std::cout << "Forces in the reference domain:" << std::endl;
            ita = anumbers.begin();
            for (std::vector<double>::iterator it = forces.begin();
                 it != forces.end(); it += 3)
            {
                std::cout << *ita;
                for (int i = 0; i < 3; i++)
                    std::cout << "    " << *(it + i);
                std::cout << std::endl;
                ita++;
            }
        }

        // rotate the forces to the original coordinate system
        transpose_rotate_PinnedH2O(positions, forces, out_of_plane_rotation_matrix, planar_rotation_angle, flipped_bond);

        // print out results
        if (MPIdata::onpe0)
        {
            std::cout << "Positions:" << std::endl;
            std::vector<short>::iterator ita = anumbers.begin();
            for (std::vector<double>::iterator it = positions.begin();
                 it != positions.end(); it += 3)
            {
                std::cout << *ita;
                for (int i = 0; i < 3; i++)
                    std::cout << "    " << *(it + i);
                std::cout << std::endl;
                ita++;
            }
            std::cout << "Forces:" << std::endl;
            ita = anumbers.begin();
            for (std::vector<double>::iterator it = forces.begin();
                 it != forces.end(); it += 3)
            {
                std::cout << *ita;
                for (int i = 0; i < 3; i++)
                    std::cout << "    " << *(it + i);
                std::cout << std::endl;
                ita++;
            }
        }

        delete mgmol;

    } // close main scope

    mgmol_finalize();

    mpirc = MPI_Finalize();
    if (mpirc != MPI_SUCCESS)
    {
        std::cerr << "MPI Finalize failed!!!" << std::endl;
    }

    time_t tt;
    time(&tt);
    if (onpe0) std::cout << " Run ended at " << ctime(&tt) << std::endl;

    return 0;
}
#endif  // MGMOL_HAS_LIBROM
