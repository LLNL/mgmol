// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#include "IonicStepper.h"
#include "Control.h"
#include "MGmol_MPI.h"
#include "MGmol_blas1.h"
#include <stdlib.h>

using namespace std;

IonicStepper::IonicStepper(const double dt, const std::vector<short>& atmove,
    std::vector<double>& tau0, std::vector<double>& taup)
    : atmove_(atmove), tau0_(tau0), taup_(taup)
{
    assert(3 * atmove.size() == tau0.size());
    assert(taup.size() == tau0.size());

    Control& ct(*(Control::instance()));

    dt_ = dt;

    ndofs_                                = 0;
    int num_movable                       = 0;
    int na                                = 0;
    std::vector<short>::const_iterator ia = atmove_.begin();
    while (ia != atmove_.end())
    {
        if (*ia) // atmove_[ia] )
        {
            ndofs_ += 3;
            num_movable++;
        }
        na++;
        ia++;
    }
    // get total dofs_
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(&ndofs_, 1, MPI_SUM);
    mmpi.allreduce(&na, 1, MPI_SUM);
    mmpi.allreduce(&num_movable, 1, MPI_SUM);

    // assumes periodic BC
    if (num_movable == na && ct.enforceVmass0) ndofs_ -= 3;
}

int IonicStepper::writeAtomicFields(HDFrestart& h5f_file,
    const vector<double>& data, const string& name, const bool create) const
{
    hid_t file_id = h5f_file.file_id();
    if (file_id < 0) return 0;

    hid_t dataspace_id = -1;
    hid_t dataset_id   = -1;

    if (create)
    {
        // Create the data space for new datasets
        hsize_t dims[2] = { (hsize_t)data.size() / 3, 3 };
        dataspace_id    = H5Screate_simple(2, dims, nullptr);
        if (dataspace_id < 0)
        {
            (*MPIdata::serr)
                << "IonicStepper::writeAtomicFields, H5Screate_simple failed!!!"
                << endl;
            return -1;
        }

        // Create the dataset
        dataset_id = H5Dcreate2(file_id, name.c_str(), H5T_NATIVE_DOUBLE,
            dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (dataset_id < 0)
        {
            (*MPIdata::serr) << "IonicStepper::writeAtomicFields, H5Dcreate2 "
                                "failed for dataset "
                             << name << "!!!" << endl;
            return -1;
        }
    }
    else
    {
        // Open the dataset
        dataset_id = H5Dopen2(file_id, name.c_str(), H5P_DEFAULT);
        if (dataset_id < 0)
        {
            (*MPIdata::serr) << "IonicStepper::writeAtomicFields, H5Dopen2 "
                                "failed for dataset "
                             << name << "!!!" << endl;
            return -1;
        }
    }

    herr_t status = H5Dwrite(
        dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);
    if (status < 0)
    {
        (*MPIdata::serr) << "IonicStepper::writeAtomicFields: H5Dwrite " << name
                         << " failed!!!" << endl;
        return -1;
    }
    else
    {
        if (onpe0)
            (*MPIdata::sout)
                << "IonicStepper::writeAtomicFields, Data written into file "
                << h5f_file.filename() << endl;
    }

    status = H5Dclose(dataset_id);
    if (status < 0)
    {
        (*MPIdata::serr)
            << "IonicStepper::writeAtomicFields, H5Dclose failed!!!" << endl;
        return -1;
    }

    if (create)
    {
        status = H5Sclose(dataspace_id);
        if (status < 0)
        {
            (*MPIdata::serr)
                << "IonicStepper::writeAtomicFields, H5Sclose failed!!!"
                << endl;
            return -1;
        }
    }

    return 0;
}

int IonicStepper::writePositions(HDFrestart& h5f_file, const string& name) const
{
    return writeAtomicFields(h5f_file, tau0_, name, false);
}

// Write velocities defined as: (taup-tau0)/dt
int IonicStepper::writeVelocities(HDFrestart& h5f_file) const
{
    hid_t file_id = h5f_file.file_id();
    if (file_id < 0) return 0;

    // Create the data space for new datasets
    hsize_t dims[2] = { (hsize_t)tau0_.size() / 3, 3 };

    hid_t dataspace_id = H5Screate_simple(2, dims, nullptr);
    if (dataspace_id < 0)
    {
        (*MPIdata::serr) << "IonicStepper: H5Screate_simple failed!!!" << endl;
        return -1;
    }

    // Create dataset
    hid_t dataset_id = H5Dcreate2(file_id, "/Ionic_velocities",
        H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0)
    {
        (*MPIdata::serr)
            << "IonicStepper:: H5Dcreate2 /Ionic_velocities failed!!!" << endl;
        return -1;
    }

    vector<double> data(taup_);
    double minus = -1.;
    int n = (int)tau0_.size(), ione = 1;
    DAXPY(&n, &minus, &tau0_[0], &ione, &data[0], &ione);
    if (dt_ > 0.)
    {
        double invdt = 1. / dt_;
        DSCAL(&n, &invdt, &data[0], &ione);
    }

    herr_t status = H5Dwrite(
        dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);
    if (status < 0)
    {
        (*MPIdata::serr) << "IonicStepper::H5Dwrite velocities failed!!!"
                         << endl;
        return -1;
    }
    else
    {
        if (onpe0)
            (*MPIdata::sout) << "Ionic velocities written into "
                             << h5f_file.filename() << endl;
    }

    status = H5Dclose(dataset_id);
    if (status < 0)
    {
        (*MPIdata::serr) << "H5Dclose failed!!!" << endl;
        return -1;
    }
    H5Sclose(dataspace_id);

    return 0;
}

int IonicStepper::readAtomicFields(
    HDFrestart& h5f_file, vector<double>& data, const string& name)
{
    hid_t file_id = h5f_file.file_id();

    // Open the dataset
    if (file_id >= 0)
    {
        hid_t dataset_id = H5Dopen2(file_id, name.c_str(), H5P_DEFAULT);
        if (dataset_id < 0)
        {
            (*MPIdata::serr) << "IonicStepper, H5Dopen2 failed for " << name
                             << " !!!" << endl;
            return -1;
        }

        assert(H5Dget_storage_size(dataset_id) == data.size() * sizeof(double));
        herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, &data[0]);
        if (status < 0)
        {
            (*MPIdata::serr) << "IonicStepper, H5Dread failed!!!" << endl;
            return -1;
        }
        // close dataset
        status = H5Dclose(dataset_id);
        if (status < 0)
        {
            (*MPIdata::serr) << "H5Dclose failed!!!" << endl;
            return -1;
        }
    }

    return 0;
}

int IonicStepper::readPositions_hdf5(HDFrestart& h5f_file, const string& name)
{
    return readAtomicFields(h5f_file, tau0_, name);
}

int IonicStepper::writeRandomStates(HDFrestart& h5f_file,
    vector<unsigned short>& data, const string& name) const
{
    hid_t file_id = h5f_file.file_id();
    bool create   = false;
    hid_t dataspace_id;

    // try to open the dataset
    hid_t dataset_id = H5Dopen2(file_id, name.c_str(), H5P_DEFAULT);
    if (dataset_id < 0)
    {
        (*MPIdata::serr) << "IonicStepper:: H5Dopen2 " << name
                         << " failed!!! Creating new data file " << endl;

        // Create the data space for the dataset
        hsize_t dims[2] = { (hsize_t)data.size() / 3, 3 };
        dataspace_id    = H5Screate_simple(2, dims, nullptr);
        if (dataspace_id < 0)
        {
            (*MPIdata::serr)
                << "Ions::writeRandomStates: H5Screate_simple failed!!!"
                << endl;
            return -1;
        }
        // Create the dataset
        hid_t dataset_id = H5Dcreate2(file_id, name.c_str(), H5T_NATIVE_USHORT,
            dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (dataset_id < 0)
        {
            (*MPIdata::serr)
                << "Ions::writeRandomStates: H5Dcreate2 failed!!!" << endl;
            return -1;
        }
        create = true;
    }

    herr_t status = H5Dwrite(
        dataset_id, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

    if (status < 0)
    {
        (*MPIdata::serr) << "IonicStepper::writeRandomStates: H5Dwrite " << name
                         << " failed!!!" << endl;
        return -1;
    }
    else
    {
        if (onpe0)
            (*MPIdata::sout)
                << "IonicStepper::writeRandomStates, Data written into file "
                << h5f_file.filename() << endl;
    }

    status = H5Dclose(dataset_id);
    if (status < 0)
    {
        (*MPIdata::serr)
            << "IonicStepper::writeRandomStates, H5Dclose failed!!!" << endl;
        return -1;
    }

    if (create)
    {
        status = H5Sclose(dataspace_id);
        if (status < 0)
        {
            (*MPIdata::serr)
                << "IonicStepper::writeAtomicFields, H5Sclose failed!!!"
                << endl;
            return -1;
        }
    }

    return 0;
}

int IonicStepper::readRandomStates(
    HDFrestart& h5f_file, vector<unsigned short>& data, const string& name)
{
    hid_t file_id = h5f_file.file_id();

    // open dataset
    if (file_id >= 0)
    {
        hid_t dataset_id = H5Dopen2(file_id, name.c_str(), H5P_DEFAULT);
        if (dataset_id < 0)
        {
            if (onpe0)
            {
                (*MPIdata::sout)
                    << "H5Dopen failed for /Ionic_RandomStates" << endl;
                (*MPIdata::sout)
                    << "Set random states to default computed in Ion.cc"
                    << endl;
            }
        }
        else
        {
            assert(H5Dget_storage_size(dataset_id)
                   == data.size() * sizeof(unsigned short));

            if (onpe0)
                (*MPIdata::sout) << "Read Ionic random states from "
                                 << h5f_file.filename() << endl;
            herr_t status = H5Dread(dataset_id, H5T_NATIVE_USHORT, H5S_ALL,
                H5S_ALL, H5P_DEFAULT, &data[0]);
            if (status < 0)
            {
                (*MPIdata::serr) << "IonicStepper: H5Dread failed!!!" << endl;
                return -1;
            }
            // close dataset
            status = H5Dclose(dataset_id);
            if (status < 0)
            {
                (*MPIdata::serr) << "H5Dclose failed!!!" << endl;
                return -1;
            }
        }
    }

    return 0;
}
