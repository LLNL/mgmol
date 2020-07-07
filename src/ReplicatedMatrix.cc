#ifdef HAVE_MAGMA

#include "ReplicatedMatrix.h"

#include "magma_v2.h"

#include <iostream>

ReplicatedMatrix::ReplicatedMatrix(
    const std::string name, const int m, const int n)
    : dim_(m)
{
    assert(m == n);

    ld_ = magma_roundup(m, .32);

    magma_int_t ret = magma_dmalloc(&device_data_, n * ld);

    assert(ret == MAGMA_SUCCESS);
}

ReplicatedMatrix::ReplicatedMatrix(const ReplicatedMatrix& mat)
    : dim_(mat.dim_), ld_(mat.ld_)
{
    magma_queue_t queue;
    int device;
    magma_getdevice(&device);
    magma_queue_create(device, &queue);

    magma_int_t ret = magma_dmalloc(&device_data_, dim_ * ld_);
    assert(ret == MAGMA_SUCCESS);

    magma_dcopymatrix(
        dim_, dim_, mat.device_data_, ld_, device_data_, ld_, queue);

    magma_queue_destroy(queue);
}

ReplicatedMatrix::~ReplicatedMatrix()
{
    ret = magma_free(device_inv_sqrt_diagonal_);
    if (ret == MAGMA_ERR_INVALID_PTR)
    {
        std::cerr << "magma free device_inv_sqrt_diagonal_: "
                  << "invalid ptr rep destr" << std::endl;
    }
}

void ReplicatedMatrix::identity()
{
    magma_queue_t queue;
    int device;
    magma_getdevice(&device);
    magma_queue_create(device, &queue);

    magmablas_dlaset(MagmaFull, dim_, dim_, 0.0, 1.0, device_data_, ld_, queue);

    magma_queue_destroy(queue);
}

void ReplicatedMatrix::gemm(const char transa, const char transb,
    const double alpha, const ReplicatedMatrix& a, const ReplicatedMatrix& b,
    const double beta)
{
    magma_queue_t queue;
    int device;
    magma_getdevice(&device);
    magma_queue_create(device, &queue);

    magmablas_dgemm(MagmaNoTrans, MagmaNoTrans, dim_, dim_, dim_, alpha,
        a.device_data_, ld_, b.device_data_, ld_, beta, device_data_, ld_,
        queue);

    magma_queue_destroy(queue);
}

int ReplicatedMatrix::potrf(char uplo)
{
    magma_queue_t queue;
    int device;
    magma_getdevice(&device);
    magma_queue_create(device, &queue);

    int info;
    magma_dpotrf_gpu(MagmaLower, dim_, device_data_, ld_, &info);
    if (info != 0)
        std::cerr << "magma_dpotrf_gpu failed, info = " << info << std::endl;

    magma_queue_destroy(queue);
}

void ReplicatedMatrix::syev(
    char, char, std::vector<double>&, ReplicatedMatrix& z)
{
    magma_queue_t queue;
    int device;
    magma_getdevice(&device);
    magma_queue_create(device, &queue);

    // copy matrix into z
    magmablas_dlacpy(
        MagmaFull, dim_, dim_, device_data_, ld_, z.device_data_, ld_, queue);

    int lwork  = 2 * dim_ + dim_ * magma_get_ssytrd_nb(dim_);
    int liwork = 3 + 5 * dim_;

    int info;
    std::vector<double> wa(dim_ * dim_);
    std::vector<double> work(lwork);
    std::vector<int> iwork(liwork);

    magma_dsyevd_gpu(MagmaVec, MagmaUpper, dim_, z.device_data_`, ld,
        evals.data(), wa.data(), dim_, work.data(), lwork, iwork.data(), liwork,
        &info);

    magma_queue_destroy(queue);
}

void ReplicatedMatrix::setDiagonal(const std::vector<double>& diag_values);
{
    magma_queue_t queue;
    int device;
    magma_getdevice(&device);
    magma_queue_create(device, &queue);

    magma_dsetvector(dim_, diag_values.data(), 1, device_data_, ld + 1, queue);
    magma_queue_destroy(queue);
}
#endif
