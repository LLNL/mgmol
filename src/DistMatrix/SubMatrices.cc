// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "SubMatrices.h"
#include "SubMatricesIndexing.h"

#include <limits.h>

//#define DEBUG
using namespace std;

namespace dist_matrix
{

// constructor based on vectors of indexes corresponding
// to "local" matrix elements
template <class T>
SubMatrices<T>::SubMatrices(const string& name,
    const vector<vector<int>>& indexes, MPI_Comm comm, const DistMatrix<T>& mat,
    const SubMatricesIndexing<T>& submat_indexing)
    : object_name_(name), comm_(comm), submat_indexing_(submat_indexing)
{
#if USE_MPI
    MPI_Comm_size(comm_, &npes_);
    MPI_Comm_rank(comm_, &mype_);
#else
    npes_ = 1;
    mype_ = 0;
#endif
    npes_distmat_ = mat.nprow() * mat.npcol();

#ifdef DEBUG
    if (mype_ == 0)
        cout << "SubMatrices with " << indexes[0].size() << " indexes" << endl;
#endif

    const int nprow  = mat.nprow();
    const int mypr   = mat.myrow();
    const int mypc   = mat.mycol();
    const int my_ind = mypr + mypc * nprow;
    send_            = (my_ind > -1);

    nb_local_matrices_ = (int)indexes.size();
    n_                 = (int)indexes[0].size();
    assert(n_ < SHRT_MAX);
    assert(n_ <= mat.m());
    assert(n_ <= mat.n());
    assert(n_ >= 0);
    assert(nb_local_matrices_ > 0);

    // allocate containers for data (matrices)
    int n2 = n_ * n_;
    local_data_.resize(nb_local_matrices_);
    if (n_ > 0)
    {
        storage_ = new T[nb_local_matrices_ * n2];
        memset(storage_, 0., nb_local_matrices_ * n2 * sizeof(T));
        for (int im = 0; im < nb_local_matrices_; im++)
            local_data_[im] = &storage_[im * n2];
    }
    else
    {
        storage_ = 0;
        for (int im = 0; im < nb_local_matrices_; im++)
            local_data_[im] = 0;
    }
}

template <class T>
void SubMatrices<T>::scal(const double alpha)
{
    int n2 = nb_local_matrices_ * n_ * n_;
    //   int incx=1;
    //   dscal(&n2, &alpha, storage_, &incx);
    MPscal(n2, alpha, storage_);
}

template <>
void SubMatrices<double>::sendrecv(double* sendbuf, double* recvbuf,
    const vector<type_displ>& my_displ, const vector<type_displ>& remote_displ)
{
    MPI_Request req;
    bool irecv = false;
    int src = mype_;

    for (int p = 0; p < npes_ - 1; p++)
    {
        // next destination/source
        int dst = (mype_ + p + 1) % npes_;
        src--;
        if (src < 0) src += npes_;

        assert(dst < npes_);
        assert(dst >= 0);

        // post receive data
        int mysize = submat_indexing_.getMySize(src);
        if (mysize)
        {
            MPI_Irecv(recvbuf + my_displ[src], mysize, MPI_DOUBLE, src, p,
                comm_, &req);
            irecv = true;
        }
        else
        {
            irecv = false;
        }

        // send data
        int rsize = submat_indexing_.getRemoteSize(dst);
        if (rsize)
        {
            MPI_Send(
                sendbuf + remote_displ[dst], rsize, MPI_DOUBLE, dst, p, comm_);
        }
        if (irecv) MPI_Wait(&req, MPI_STATUS_IGNORE);
    }
}

template <>
void SubMatrices<float>::sendrecv(float* sendbuf, float* recvbuf,
    const vector<type_displ>& my_displ, const vector<type_displ>& remote_displ)
{
    MPI_Request req;
    bool irecv = false;
    int src = mype_;

    for (int p = 0; p < npes_ - 1; p++)
    {
        // next destination/source
        int dst = (mype_ + p + 1) % npes_;
        src--;
        if (src < 0) src += npes_;

        assert(dst < npes_);
        assert(dst >= 0);

        // post receive data
        int mysize = submat_indexing_.getMySize(src);
        if (mysize)
        {
            MPI_Irecv(recvbuf + my_displ[src], mysize, MPI_FLOAT, src, p,
                comm_, &req);
            irecv = true;
        }
        else
        {
            irecv = false;
        }

        // send data
        int rsize = submat_indexing_.getRemoteSize(dst);
        if (rsize)
        {
            MPI_Send(
                sendbuf + remote_displ[dst], rsize, MPI_FLOAT, dst, p, comm_);
        }
        if (irecv) MPI_Wait(&req, MPI_STATUS_IGNORE);
    }
}

// from mat to local_data_
template <class T>
void SubMatrices<T>::gather(const DistMatrix<T>& mat)
{
#ifdef DEBUG
    if (mype_ == 0) cout << "SubMatrices<T>::gather()" << endl;
#endif
    gather_tm_.start();

    // sizes
    int tot_remote_size = 0;
    for (int p = 0; p < npes_; p++)
        tot_remote_size += submat_indexing_.getRemoteSize(p);
    int my_size = 0;
    for (int p = 0; p < npes_distmat_; p++)
        my_size += submat_indexing_.getMySize(p);

    // compute displacements
    vector<type_displ> remote_displ(npes_, 0);
    for (int pe = 1; pe < npes_; pe++)
        remote_displ[pe] = (type_displ)(
            submat_indexing_.getRemoteSize(pe - 1) + remote_displ[pe - 1]);
    vector<type_displ> my_displ(npes_distmat_, 0);
    for (int pe = 1; pe < npes_distmat_; pe++)
        my_displ[pe] = (type_displ)(
            submat_indexing_.getMySize(pe - 1) + my_displ[pe - 1]);

    // build buffer for data array to send
    vector<T> buf_remote_val(tot_remote_size, 0.);
    if (send_)
        for (int pe = 0; pe < npes_; pe++)
        {
            const int displ = remote_displ[pe];
            assert(displ >= 0);
            const int size = submat_indexing_.getRemoteSize(pe);
            assert(size >= 0);
            // initialize buffer to send to PE pe
            for (int k = 0; k < size; k++)
            {
                const int ii = 2 * (displ + k);
                const int i  = submat_indexing_.getRemoteDoubleIndex(ii);
                const int j  = submat_indexing_.getRemoteDoubleIndex(ii + 1);
#ifdef DEBUG
                cout << "PE: " << mype_ << ", (i,j)=(" << i << "," << j << ")"
                     << endl;
#endif
                assert(displ + k < (int)buf_remote_val.size());
                assert(j < mat.n());
                assert(i < mat.m());
                assert(i >= 0);
                assert(j >= 0);
                buf_remote_val[displ + k] = mat.getGlobalVal(i, j, true);
            }
        }

    gather_comm_tm_.start();

    vector<T> buf_my_val(my_size, 0.);

#ifdef USE_MPI

#if 0
  // send/recv data 
  MPI_Alltoallv(&buf_remote_val[0],&remote_sizes_[0],&remote_displ[0],MPI_DOUBLE,
                &buf_my_val[0],    &my_sizes_[0],    &my_displ[0],    MPI_DOUBLE,
                comm_);
#else
    // The main interest of the MPI_Isend operation is the possible overlap
    // of computation with communications,
    // or the possible overlap between multiple communications.
    // For large messages, in order to keep the memory usage on the receiver
    // at a reasonable level, a rendezvous protocol is used. The sender
    //[after sending a small packet] wait until the receiver confirm the
    // message exchange (i.e. the corresponding receive operation has been
    // posted) to send the large data. Using MPI_Isend can lead to longer
    // execution times, as the real transfer will be delayed until the
    // program enter in the next MPI call.

    T* recvbuf = &buf_my_val[0];
    T* sendbuf = &buf_remote_val[0];

    int nelements = submat_indexing_.getMySize(mype_);
    if (nelements > 0)
        memcpy(recvbuf + my_displ[mype_], sendbuf + remote_displ[mype_],
            nelements * sizeof(T));

    int src = mype_;

    sendrecv(sendbuf,recvbuf,my_displ,remote_displ);

#endif

#else
    assert((int)buf_remote_val.size() == nb_local_matrices_ * n_ * n_);
    assert((int)buf_my_val.size() == nb_local_matrices_ * n_ * n_);
    buf_my_val = buf_remote_val;
#endif // USE_MPI

    gather_comm_tm_.stop();
    gather_comp_tm_.start();

    // loop over matrices to assign received data
    if (buf_my_val.size() > 0)
        for (int im = 0; im < nb_local_matrices_; im++)
        {
            assert(im < (int)local_data_.size());
            T* plocal_data = local_data_[im];
            // loop over indexes (2 int)
            // index ii and ii+1
            for (int i = 0; i < n_; i++)
            {
                for (int j = 0; j < n_; j++)
                {
                    const int il = i + n_ * j;
                    const int ii = submat_indexing_.getVectorIndex(im, il);
                    // if( !(ii<(int)buf_my_val.size()) ){
                    //  cout<<"ii="<<ii<<",
                    //  buf_my_val.size()="<<buf_my_val.size()<<endl;
                    //}
                    assert(ii < (int)buf_my_val.size());
                    if (ii > -1) plocal_data[il] = buf_my_val[ii];
                }
            }
        }

    gather_comp_tm_.stop();

    gather_tm_.stop();
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
void SubMatrices<T>::print(ostream& os)
{
    for (int im = 0; im < nb_local_matrices_; im++)
    {
        os << "SubMatrices " << n_ << " x " << n_ << ", name=" << object_name_
           << ", matrix " << im << endl;
        os << "Indexes: ";
        for (int i = 0; i < n_; i++)
            os << submat_indexing_.getLocalIndex(im, i) << "\t";
        os << endl;

        for (int i = 0; i < n_; i++)
        {
            for (int j = 0; j < n_; j++)
                os << local_data_[im][i + n_ * j] << "\t";
            os << endl;
        }
        os << endl;
    }
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
void SubMatrices<T>::get_array(const int im, vector<T>& val)
{
    assert(im < nb_local_matrices_);
    const map<int, short>& mapind = submat_indexing_.mapIndexes(im);
    val.resize(n_ * n_);

    const T* plocal_data = local_data_[im];

    for (int i = 0; i < n_; i++)
    {
        map<int, short>::const_iterator p0
            = mapind.find(submat_indexing_.getLocalIndex(im, i));
        assert(p0 != mapind.end());
        short il = p0->second;
        assert(il < n_);
        for (int j = 0; j < n_; j++)
        {
            map<int, short>::const_iterator p1
                = mapind.find(submat_indexing_.getLocalIndex(im, j));
            assert(p1 != mapind.end());
            short jl = p1->second;
            assert(jl < n_);

            val[il + n_ * jl] = plocal_data[i + n_ * j];
        }
    }
}

template class SubMatrices<double>;
template class SubMatrices<float>;

} // namespace
