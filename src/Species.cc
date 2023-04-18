// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Species.h"
#include "Control.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"
#include "read.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#define Ha2Ry 2.

void Species::read_1species(const std::string& filename)
{
    Control& ct     = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    const int size_buf = 64;
    char buf_name[size_buf];
    char buf_ignored[size_buf];
    memset(buf_name, 0, size_buf * sizeof(char));

    std::ifstream* tfile = nullptr;
    if (mmpi.instancePE0())
    {
        tfile = new std::ifstream(filename.data(), std::ios::in);
        if (tfile->fail())
        {
            (*MPIdata::serr)
                << " Cannot open file " << filename.data() << std::endl;
            mmpi.abort();
        }
        else
        {
            if (ct.verbose > 0)
                (*MPIdata::sout)
                    << " read_1species: Open " << filename << std::endl;
        }

        // Description
        read_comments(*tfile);
        tfile->get(buf_name, size_buf);
        read_comments(*tfile);

#ifdef DEBUG
        (*MPIdata::sout) << " read_1species: Species = " << buf_name
                         << std::endl;
#endif
        tfile->get(buf_ignored, size_buf); // unused
        read_comments((*tfile));

        double ignored1, ignored2;
        (*tfile) >> ignored1 >> ignored2;
        read_comments((*tfile));

        // Read in the nlcc flag
        short nlccflag;
        read_data(nlccflag, (*tfile));
        if (nlccflag != 0)
        {
            (*MPIdata::serr) << " NLCC Option not implemented!!!" << std::endl;
            mmpi.abort();
        }

        // Read in the atomic number
        read_data(atomic_number_, (*tfile));

        // Read in the atomic mass
        read_data(mass_, (*tfile));

        // Read in the number of valence electrons
        read_data(zion_, (*tfile));
        assert(zion_ < 50);

        /* Gaussian charge parameter */
        read_data(rc_, (*tfile));
        assert(rc_ > 0.);

        /* Number of potentials */
        read_data(num_potentials_, (*tfile));

        if (num_potentials_ <= 0)
        {
            (*MPIdata::serr) << " Species: num_potentials_=" << num_potentials_
                             << " Need potential functions" << std::endl;
            mmpi.abort();
        }
        if (num_potentials_ > 5)
        {
            (*MPIdata::serr) << " Species: num_potentials_=" << num_potentials_
                             << " Too many potential functions" << std::endl;
            mmpi.abort();
        }

        // L-value for the local potential
        (*tfile) >> llocal_ >> type_flag_;
        read_comments((*tfile));
#ifdef DEBUG
        (*MPIdata::sout) << " Potential " << llocal_ << " is local"
                         << std::endl;
        if (type_flag_ == 1)
            (*MPIdata::sout) << " Divide ref. state by radius" << std::endl;
        if (type_flag_ == 2)
            (*MPIdata::sout) << " Divide projector by radius" << std::endl;
#endif
        assert(llocal_ < num_potentials_);
        assert(llocal_ >= 0);

        /* Local potential radius */
        read_data(lradius_, (*tfile));

        /* Non-Local potential radius */
        read_data(nlradius_, (*tfile));

        /* Number of points in the radial grid */
        read_data(n_rad_points_, (*tfile));
    }

    // Broadcast information on pseudopotentials
    MPI_Bcast(&atomic_number_, 1, MPI_SHORT, 0, comm_);
    MPI_Bcast(&zion_, 1, MPI_SHORT, 0, comm_);
    MPI_Bcast(&n_rad_points_, 1, MPI_INT, 0, comm_);
    MPI_Bcast(&rc_, 1, MPI_DOUBLE, 0, comm_);
    MPI_Bcast(&lradius_, 1, MPI_DOUBLE, 0, comm_);
    MPI_Bcast(&nlradius_, 1, MPI_DOUBLE, 0, comm_);
    MPI_Bcast(&mass_, 1, MPI_DOUBLE, 0, comm_);
    MPI_Bcast(&num_potentials_, 1, MPI_SHORT, 0, comm_);
    MPI_Bcast(&llocal_, 1, MPI_SHORT, 0, comm_);
    MPI_Bcast(&type_flag_, 1, MPI_SHORT, 0, comm_);

    setRcDependentData();

    input_kbp_.resize(num_potentials_);
    multiplicity_.resize(num_potentials_);
    kb_coeff_.resize(num_potentials_);
    ekb_.resize(num_potentials_);
    kb_sign_.resize(num_potentials_);

    if (mmpi.instancePE0())
    {
        if (type_flag_ == 2 || type_flag_ == 3) // multiple projectors
        {
            for (short l = 0; l < num_potentials_; l++)
            {
                assert(l < static_cast<int>(multiplicity_.size()));
                if (l != llocal_)
                {
                    std::string query;
                    getline(*tfile, query);

                    std::stringstream ss(query);
                    short ll, nproj;
                    ss >> ll >> nproj;
                    assert(nproj > 0);
                    multiplicity_[ll] = nproj;

                    assert(ll < static_cast<int>(ekb_.size()));

                    ekb_[ll].resize(multiplicity_[ll]);
                    for (short p = 0; p < nproj; ++p)
                    {
                        assert(p < static_cast<int>(ekb_[ll].size()));
                        ss >> ekb_[ll][p];
                        std::cout << "Read ekb_[" << ll << "][" << p
                                  << "]=" << ekb_[ll][p] << std::endl;
                    }
                }
                else
                {
                    multiplicity_[l] = 0;
                }
            }

            h1s_ = 0.; // not GTH potential
        }
        else // single KB projectors
        {
            for (short l = 0; l < num_potentials_; l++)
            {
                assert(l < static_cast<int>(multiplicity_.size()));
                if (l != llocal_)
                    multiplicity_[l] = 1;
                else
                    multiplicity_[l] = 0;
            }
            // read Goedecker's potentials info
            (*tfile) >> h1s_;
            if (h1s_ > 1.e-12) (*tfile) >> h2s_ >> h1p_;
        }

        read_comments((*tfile));

#ifdef DEBUG
        (*MPIdata::sout) << " h1s=" << h1s_ << std::endl;
        (*MPIdata::sout) << " atomic number " << atomic_number_ << std::endl;
        (*MPIdata::sout) << " atomic mass " << mass_ << std::endl;
        (*MPIdata::sout) << " rc=" << rc_ << std::endl;
        (*MPIdata::sout) << num_potentials_ << " potentials\n";
        (*MPIdata::sout) << " type_flag=" << type_flag_ << std::endl;
        (*MPIdata::sout) << " Local potential radius=" << lradius_ << std::endl;
        (*MPIdata::sout) << " Non-local potential radius=" << nlradius_
                         << std::endl;
#endif
    }

    int mpirc = MPI_Bcast(buf_name, size_buf, MPI_CHAR, 0, comm_);
    if (mpirc != MPI_SUCCESS)
    {
        (*MPIdata::serr) << "Species, MPI Bcast of buf_name failed!!!"
                         << std::endl;
        mmpi.abort();
    }

    MPI_Bcast(&h1s_, 1, MPI_DOUBLE, 0, comm_);
    if (h1s_ > 1.e-12)
    {
        MPI_Bcast(&h2s_, 1, MPI_DOUBLE, 0, comm_);
        MPI_Bcast(&h1p_, 1, MPI_DOUBLE, 0, comm_);
    }
    assert(n_rad_points_ > 0);
    assert(num_potentials_ > 0);
    assert(llocal_ < num_potentials_);

#ifdef DEBUG
    if (printFlag)
        (*MPIdata::sout) << " n_rad_points_=" << n_rad_points_ << std::endl;
#endif
    name_ = buf_name;
    MPI_Bcast(&multiplicity_[0], num_potentials_, MPI_SHORT, 0, comm_);
    for (short l = 0; l < num_potentials_; l++)
    {
        kb_coeff_[l].resize(multiplicity_[l]);
        ekb_[l].resize(multiplicity_[l]);
        kb_sign_[l].resize(multiplicity_[l]);
    }

    if (type_flag_ == 2 || type_flag_ == 3)
    {
        for (short l = 0; l < num_potentials_; l++)
        {
            if (l != llocal_)
            {
                assert(multiplicity_[l] > 0);
                MPI_Bcast(&ekb_[l][0], multiplicity_[l], MPI_DOUBLE, 0, comm_);
                assert(fabs(ekb_[l][0]) > 0.);
            }
        }
    }

    max_l_ = num_potentials_ - 1;

    if (h1s_ == 0.)
        readRadialKBpotentials(tfile);
    else
        readRadialGTHpotentials(tfile);

    if (mmpi.instancePE0())
    {
        if (tfile != nullptr)
        {
            tfile->close();
            delete tfile;
        }
    }
#ifdef DEBUG
    (*MPIdata::sout) << "Potential read" << std::endl;
#endif

    lradius_  = std::min(lradius_, input_kbp_[0].x(n_rad_points_ - 1));
    nlradius_ = std::min(nlradius_, input_kbp_[0].x(n_rad_points_ - 1));

    kbp_.resize(num_potentials_);
    for (short l = 0; l < num_potentials_; l++)
    {
        if (l != llocal_) kbp_[l].resize(multiplicity_[l]);
    }

    if (mmpi.instancePE0()) checkLRadius();

    assert(zion_ < 50);
    assert(llocal_ < num_potentials_);
    assert(num_potentials_ > 0);
    assert(rc_ > 0.);
}

void Species::readRadialKBpotentials(std::ifstream* tfile)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    if (mmpi.instancePE0())
    {
        // Next read in the radial grid, potential, and reference state
        for (short l = 0; l < num_potentials_; l++)
        {
            if (l != llocal_) assert(multiplicity_[l] > 0);
#ifdef DEBUG
            (*MPIdata::sout) << "Read potential " << l << std::endl;
#endif
            short npot = multiplicity_[l] > 1 ? multiplicity_[l] : 2;
            short ncol = (l != llocal_) ? 1 + npot : 2;

            input_kbp_[l].read(n_rad_points_, ncol, tfile);

            if (l < (num_potentials_ - 1)) read_comments((*tfile));

            if (type_flag_ == 1 && (num_potentials_ > 1))
            {
                (*MPIdata::sout) << "Divide ref. states by radius" << std::endl;
                const std::vector<double>& rps(input_kbp_[l].x());
                std::vector<double>& rpsi(ref_psi(l));
                int kmin = 0;
                if (rps[0] < 1.e-10)
                {
                    rpsi[0] = 2. * rpsi[1] - rpsi[2];
                    kmin++;
                }
                for (int k = kmin; k < n_rad_points_; k++)
                {
                    assert(rps[k] >= 1.e-10);
                    rpsi[k] /= rps[k];
                }
            }
        }
    }

    for (short l = 0; l < num_potentials_; l++)
        input_kbp_[l].bcast(comm_, 0);
}

void Species::readRadialGTHpotentials(std::ifstream* tfile)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    if (mmpi.instancePE0())
    {
        // Next read in the radial grid, potential, and reference state
        for (short j = 0; j < num_potentials_; j++)
        {
#ifdef DEBUG
            (*MPIdata::sout) << "Read potential " << j << std::endl;
#endif
            short ncol  = 2;
            short index = j;

            if (j == 0) index = num_potentials_ - 1;
            if (j == 1) index = 0;

            input_kbp_[index].read(n_rad_points_, ncol, tfile);

            if (j < (num_potentials_ - 1)) read_comments((*tfile));

            if (type_flag_ == 1 && (num_potentials_ > 1))
            {
                (*MPIdata::sout) << "Divide ref. states by radius" << std::endl;
                const std::vector<double>& rps = input_kbp_[j].x();
                std::vector<double>& rpsi      = ref_psi(j);
                int kmin                       = 0;
                if (rps[0] < 1.e-10)
                {
                    rpsi[0] = 2. * rpsi[1] - rpsi[2];
                    kmin++;
                }
                for (int k = kmin; k < n_rad_points_; k++)
                {
                    assert(rps[k] >= 1.e-10);
                    rpsi[k] /= rps[k];
                }
            }
        }
    }

    for (short j = 0; j < num_potentials_; j++)
        input_kbp_[j].bcast(comm_, 0);
}

void Species::print(std::ostream& os) const
{
    os << " Atom: " << name_ << std::endl;
    os << " Atomic mass     = " << mass_ << std::endl;
    os << " Zion            = " << zion_ << std::endl;
    os << " Core rc         = " << rc_ << std::endl;
    os << " # of potentials = " << num_potentials_ << std::endl;
    os << " Local pot.      = " << llocal_ << std::endl;

    for (short ll = 0; ll < num_potentials_; ll++)
        for (short pp = 0; pp < multiplicity_[ll]; pp++)
        {
            os << "\t EKB[l=" << ll << "][" << pp << "]= ";

            if (ll == llocal_)
            {
                os << "N/A" << std::endl;
            }
            else
            {
                os << ekb_[ll][pp] << std::endl;
            }
        }
}

void Species::setKBcoeff(const short l, const short p, const double norm)
{
    assert(l < static_cast<int>(ekb_.size()));

    ekb_[l][p] = norm;
    if (norm < 0.)
        kb_sign_[l][p] = -1;
    else
        kb_sign_[l][p] = 1;

    if (l != llocal_)
    {
        //(*MPIdata::sout)<<" kbnorm["<<l<<"]["<<p<<"]= "<<ekb_[l][p]<<endl;
        assert(fabs(ekb_[l][p]) > 0.00001);
        kb_coeff_[l][p] = sqrt(fabs(ekb_[l][p]) * Ha2Ry);
    }
}

void Species::gauss_filter_local_pot(const double rcut, const bool printFlag)
{
    if (printFlag) (*MPIdata::sout) << "Filter local potential" << std::endl;
    local_pot_.gauss_filter(rcut);
}

void Species::gauss_filter_kbp(
    const short l, const short p, const double rcut, const bool printFlag)
{
    if (printFlag)
        (*MPIdata::sout) << "Filter non-local potential" << std::endl;
    assert(l < static_cast<int>(kbp_.size()));
    assert(p < static_cast<int>(kbp_[l].size()));

    kbp_[l][p].gauss_filter(rcut);
}

void Species::initLocalPotential(const char flag_filter, const double hmax,
    std::ofstream* tfile, const bool printFlag)
{
    Control& ct = *(Control::instance());

    const double lrcut    = lradius_ * 0.75;
    const bool print_flag = (printFlag && !ct.restart_run);

    const std::vector<double>& potl(input_kbp_[llocal_].y(0));
    const std::vector<double>& rps(input_kbp_[llocal_].x());

    // output local pseudopotential
    if (printFlag && (tfile != nullptr))
    {
        input_kbp_[llocal_].printLocalPot(name_, llocal_, tfile);
    }

    // Generate the local neutralized pseudopotential
    RadialInter rwork(rps);
    std::vector<double>& work(rwork.y(0));
    assert(work.size() > 0);

    for (int idx = 0; idx < n_rad_points_; idx++)
    {
        double r  = rps[idx];
        work[idx] = potl[idx] + getVcomp(r);
    }

    // output local projector
    if (printFlag && (tfile != nullptr))
    {
        (*tfile) << "\"" << name_ << ", neutralized potential, l=" << llocal_
                 << "\"" << std::endl;

        for (int idx = 0; idx < n_rad_points_; idx += 3)
        {
            (*tfile) << rps[idx] << "\t" << work[idx] << std::endl;
        }
        (*tfile) << std::endl;
        (*tfile) << std::endl;
    }

    if (flag_filter == 'f')
    {
        // Filter it
        RadialInter rfunc(rps);
        rwork.rft(rfunc, 0, hmax, print_flag, 0, printFlag);
        assignLocalPot(rfunc.x(), rfunc.y());
        // Now cut it off smoothly in real space
        gauss_filter_local_pot(lrcut, (printFlag && ct.verbose > 2));
        if (printFlag && (tfile != nullptr))
        {
            if (ct.verbose > 0)
                std::cout << "Write local potential in file..." << std::endl;
            const RadialInter& plocal_pot = local_pot();
            (*tfile) << "\"" << name_ << " filtered: " << llocal_ << "\""
                     << std::endl;
            plocal_pot.print(*tfile);
            (*tfile) << std::endl;
            (*tfile) << std::endl;
        }
    }
    else
    {
        assignLocalPot(rps, work);
    }
}

void Species::initPotentials(
    const char flag_filter, const double hmax, const bool printFlag)
{
    assert(flag_filter == 'n' || flag_filter == 's' || flag_filter == 'f');
    assert(n_rad_points_ < 10000);
    assert(n_rad_points_ > 0);

    Control& ct          = *(Control::instance());
    std::ofstream* tfile = nullptr;
    if (printFlag)
    {
        int mpi_rank;
        MPI_Comm_rank(comm_, &mpi_rank);
        if (ct.verbose > 0)
            (*MPIdata::sout)
                << "Initialize radial potentials for species " << name_
                << " on MPI rank " << mpi_rank << std::endl;
        //(*MPIdata::sout)<<" initPotentials: potential local:"<<llocal_<<endl;
        //(*MPIdata::sout)<<"Max. l="<<max_l_<<endl;
    }
    if (printFlag && !ct.restart_run)
    {
        std::string name     = "projectors" + name_ + ".dat";
        std::string filename = ct.getFullFilename(name);
        tfile = new std::ofstream(filename.data(), std::ios::out);
    }

    initLocalPotential(flag_filter, hmax, tfile, printFlag);

    if (type_flag_ == 2)
        initNonlocalMultiProjectorsPotentials(
            flag_filter, hmax, tfile, true, printFlag);
    else if (type_flag_ == 3)
        initNonlocalMultiProjectorsPotentials(
            flag_filter, hmax, tfile, false, printFlag);
    else if (h1s_ == 0.)
        initNonlocalKBPotentials(flag_filter, hmax, tfile, printFlag);
    else
        initNonlocalGTHPotentials(flag_filter, hmax, tfile, printFlag);

    if (printFlag && !ct.restart_run)
    {
        assert(tfile != nullptr);
        tfile->close();
        delete tfile;
    }
}

void Species::initNonlocalKBPotentials(const char flag_filter,
    const double hmax, std::ofstream* tfile, const bool printFlag)
{
    assert(max_l_ >= 0);

    Control& ct = *(Control::instance());

    const bool print_flag = (printFlag && !ct.restart_run);
    const double nlrcut   = nlradius_ * 0.75;

    const std::vector<double>& potloc(input_kbp_[llocal_].y(0));

    // Loop over over angular momentum states
    for (short ll = 0; ll <= max_l_; ll++)
    {
        if (ll != llocal_)
        {
            const std::vector<double>& potl(input_kbp_[ll].y(0));
            const std::vector<double>& rps(input_kbp_[ll].x());
            RadialInter rwork(rps);
            std::vector<double>& work(rwork.y(0));
            assert(work.size() > 0);

            const std::vector<double>& refl(ref_psi(ll));
            assert(n_rad_points_ == (int)refl.size());

            if (printFlag && ct.verbose > 2)
                (*MPIdata::sout) << "Generate radial KB projector" << std::endl;
            // Loop over radial grid points to build radial KB projector
            for (int idx = 0; idx < n_rad_points_; idx++)
            {
                work[idx] = refl[idx] * (potl[idx] - potloc[idx]);
            }
            if (printFlag && (tfile != nullptr))
            {
                (*tfile) << "\"" << name_ << ", projector l=" << ll << "\""
                         << std::endl;
                // output non-local projector
                for (int idx = 0; idx < n_rad_points_; idx += 3)
                {
                    (*tfile) << rps[idx] << "\t" << work[idx] << std::endl;
                }
                (*tfile) << std::endl;
                (*tfile) << std::endl;
            }

            if (flag_filter == 'f')
            {
                // Filter it
                RadialInter rfunc(rps);
                rwork.rft(rfunc, ll, hmax, print_flag, 0, printFlag);
                assign_kbp(ll, 0, rfunc.x(), rfunc.y());
                // Now cut it off smoothly in real space
                gauss_filter_kbp(ll, 0, nlrcut, (printFlag && ct.verbose > 2));
            }
            else
            {
                assign_kbp(ll, 0, rps, work);
            }
            // output non-local projector
            if (printFlag && (tfile != nullptr))
            {
                (*tfile) << "\"" << name_ << " filtered projector, l=" << ll
                         << "\"" << std::endl;
                print_kbp(*tfile, ll);
            }

            // Evaluate the normalization coefficient for the projector
            // <phi, delta V*phi>
            for (int idx = 0; idx < n_rad_points_; idx++)
            {
                work[idx] *= refl[idx];
            }
            double factor = 1. / (rwork.radint());
            setKBcoeff(ll, 0, factor);
        }
    }
}

void Species::initNonlocalMultiProjectorsPotentials(const char flag_filter,
    const double hmax, std::ofstream* tfile, const bool divide_by_r,
    const bool printFlag)
{
    assert(max_l_ >= 0);

    Control& ct = *(Control::instance());

    if (printFlag && ct.verbose > 2)
    {
        (*MPIdata::sout) << "Species::initNonlocalMultiProjectorsPotentials()"
                         << std::endl;
    }

    const bool print_flag = (printFlag && !ct.restart_run);
    const double nlrcut   = nlradius_ * 0.75;

    // Loop over over angular momentum states
    for (short ll = 0; ll <= max_l_; ll++)
    {
        if (ll != llocal_)
        {
            assert(multiplicity_[ll] > 0);
            for (short p = 0; p < multiplicity_[ll]; ++p)
            {
                assert(ll < static_cast<int>(input_kbp_.size()));

                const std::vector<double>& pot(input_kbp_[ll].y(p));
                const std::vector<double>& rps(input_kbp_[ll].x());
                RadialInter rwork(rps);
                std::vector<double>& work(rwork.y(0));
                assert(work.size() > 0);

                // if( printFlag ) (*MPIdata::sout)<<"Copy radial
                // projector"<<endl;
                for (int idx = 0; idx < n_rad_points_; idx++)
                {
                    work[idx] = pot[idx];
                }

                // divide radial function by r
                if (divide_by_r)
                {
                    int kmin = 0;
                    if (rps[0] < 1.e-10)
                    {
                        work[0] = 2. * work[1] - work[2];
                        kmin++;
                    }
                    for (int k = kmin; k < n_rad_points_; k++)
                    {
                        assert(rps[k] >= 1.e-10);
                        work[k] /= rps[k];
                    }
                }

                if (printFlag && (tfile != nullptr))
                {
                    (*tfile) << "\"" << name_ << ", projector, l=" << ll
                             << ", p=" << p << "\"" << std::endl;
                    // output non-local projector
                    for (int idx = 0; idx < n_rad_points_; idx += 3)
                    {
                        (*tfile) << rps[idx] << "\t" << work[idx] << std::endl;
                    }
                    (*tfile) << std::endl;
                    (*tfile) << std::endl;
                }

                if (flag_filter == 'f')
                {
                    // Filter it
                    RadialInter rfunc(rps);
                    rwork.rft(rfunc, ll, hmax, print_flag, 0, printFlag);
                    assign_kbp(ll, p, rfunc.x(), rfunc.y());
                    // Now cut it off smoothly in real space
                    gauss_filter_kbp(
                        ll, p, nlrcut, (printFlag && ct.verbose > 2));
                }
                else
                {
                    assign_kbp(ll, p, rps, work);
                }

                // output non-local projector
                if (printFlag && (tfile != nullptr) && flag_filter == 'f')
                {
                    (*tfile) << "\"" << name_ << " filtered projector, l=" << ll
                             << ", p=" << p << "\"" << std::endl;
                    print_kbp(*tfile, ll, p);
                }

                // Evaluate the normalization coefficient for the projector
                // double  norm2 = rwork.radintf2();
                // if( printFlag )cout<<"Normalization factor for projector
                // l="<<ll<<", p="<<p<<": "<<norm2<<endl;
                // kbp_[ll][p].scale(sqrt(1./norm2));
                if (printFlag)
                    std::cout
                        << "Norm2 radial KB projector=" << kbp_[ll][p].norm2()
                        << std::endl;
                setKBcoeff(ll, p, ekb_[ll][p]);
            }
        }
    }
}

void Species::initNonlocalGTHPotentials(const char flag_filter,
    const double hmax, std::ofstream* tfile, const bool printFlag)
{
    assert(max_l_ >= 0);

    Control& ct = *(Control::instance());

    const bool print_flag = (printFlag && !ct.restart_run);
    const double nlrcut   = nlradius_ * 0.75;

    // Loop over over angular momentum states
    for (short ll = 0; ll <= max_l_; ll++)
    {
        if (ll != llocal_)
        {
            const std::vector<double>& potl(input_kbp_[ll].y(0));
            const std::vector<double>& rps(input_kbp_[ll].x());
            RadialInter rwork(rps);
            std::vector<double>& work(rwork.y(0));
            assert(work.size() > 0);

            if (printFlag && ct.verbose > 2)
                (*MPIdata::sout) << "Copy radial KB projector" << std::endl;
            for (int idx = 0; idx < n_rad_points_; idx++)
            {
                work[idx] = potl[idx];
            }

            if (flag_filter == 'f')
            {
                // Filter it
                RadialInter rfunc(rps);
                rwork.rft(rfunc, ll, hmax, print_flag, 0, printFlag);
                assign_kbp(ll, 0, rfunc.x(), rfunc.y());
                // Now cut it off smoothly in real space
                gauss_filter_kbp(ll, 0, nlrcut, (printFlag && ct.verbose > 2));
            }
            else
            {
                assign_kbp(ll, 0, rps, work);
            }

            // output non-local projector
            if (printFlag && (tfile != nullptr))
            {
                (*tfile) << "\"" << name_ << " filtered, l=" << ll << "\""
                         << std::endl;
                print_kbp(*tfile, ll);
            }

            double kbnorm = h1s_;
            setKBcoeff(ll, 0, kbnorm);
        }
    }
}

// set KB projectors on all tasks based on values in root
void Species::syncKBP(const int root)
{
    for (short l = 0; l <= max_l_; l++)
        if (l != llocal_)
        {
            MPI_Bcast(&ekb_[l][0], multiplicity_[l], MPI_DOUBLE, root, comm_);
            MPI_Bcast(
                &kb_coeff_[l][0], multiplicity_[l], MPI_DOUBLE, root, comm_);
            MPI_Bcast(
                &kb_sign_[l][0], multiplicity_[l], MPI_SHORT, root, comm_);
            for (short p = 0; p < multiplicity_[l]; ++p)
                kbp_[l][p].bcast(comm_, root);
        }
    local_pot_.bcast(comm_, root);
}

void Species::setRcDependentData()
{
    assert(rc_ > 0.);

    invrc_ = 1. / rc_;

    const double pi3half = M_PI * sqrt(M_PI);
    const double rcnorm  = rc_ * rc_ * rc_ * pi3half;
    assert(rcnorm > 1.e-8);

    comp_charge_factor_ = (double)zion_ / rcnorm;
}

void Species::checkLRadius() const
{
    const double rtol = 1.e-4;

    double rhoc0 = getRhoComp(0.);
    double rhocr = getRhoComp(lradius_);
    if (fabs(rhocr / rhoc0) > rtol)
    {
        std::cout << "WARNING: radius for species " << name_
                  << " is too small\n";
        std::cout << "WARNING: rhoc(" << lradius_ << ") = " << rhocr
                  << std::endl;
    }

    // check neutralizing potential close to Z/r at lradius_
    double vcr = getVcomp(lradius_);
    double zor = (double)zion_ / lradius_;
    if (fabs((vcr - zor) / zor) > rtol)
    {
        std::cout << "WARNING: radius for species " << name_
                  << " is too small\n";
        std::cout << "WARNING: vc(" << lradius_ << ") = " << vcr
                  << " z/r = " << zor << std::endl;
    }
}
