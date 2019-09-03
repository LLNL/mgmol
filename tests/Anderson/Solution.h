#ifndef MGMOL_SOLUTIONH
#define MGMOL_SOLUTIONH

#include "MGmol_blas1.h"

#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "boost/shared_ptr.hpp"
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

class Solution
{
    std::string name_;
    std::vector<double> u_;
    int n_;

    static double invs_;

public:
    Solution() { n_ = 0; }

    ~Solution() { n_ = 0; }

    Solution(const int n)
    {
        n_ = n;
        u_.resize(n);
        invs_ = 1.;
    }

    Solution(std::string name, const Solution& y)
    {
        name_ = std::move(name);
        n_    = y.n_;
        u_    = y.u_;
    }

    Solution& operator=(const Solution& y)
    {
        n_ = y.n_;
        u_ = y.u_;
        return *this;
    }

    void assign(const Solution& y) { u_ = y.u_; }

    void setInvS(const double invs) { invs_ = invs; }

    double norm()
    {
        int ione     = 1;
        double norm2 = DDOT(&n_, &u_[0], &ione, &u_[0], &ione);
        return sqrt(norm2);
    }

    void normalize()
    {
        const double my_norm = norm();
        scal(1. / my_norm);
    }

    double dotProduct(const Solution& v)
    {
        int ione = 1;
        return invs_ * DDOT(&n_, &u_[0], &ione, &v.u_[0], &ione);
    }

    Solution& operator-=(const Solution& y)
    {
        for (int i = 0; i < n_; i++)
            u_[i] -= y.u_[i];
        return *this;
    }

    void axpy(const double alpha, const Solution& y)
    {
        int ione = 1;
        DAXPY(&n_, &alpha, &y.u_[0], &ione, &u_[0], &ione);
    }

    void scal(const double alpha)
    {
        int ione = 1;
        DSCAL(&n_, &alpha, &u_[0], &ione);
    }

    double operator[](const int i) { return u_[i]; }

    void setVal(const int i, const double val) { u_[i] = val; }
    void initRand(const int imax)
    {
        typedef boost::minstd_rand rng_type;
        typedef boost::uniform_real<> distribution_type;
        typedef boost::variate_generator<rng_type, distribution_type> gen_type;

        rng_type rng(0);
        distribution_type nd(-0.5, 0.5);

        std::unique_ptr<boost::variate_generator<rng_type, distribution_type>>
            gen(new gen_type(rng, nd));

        for (int i = 0; i < imax; i++)
            u_[i] = (*gen)();
        for (int i = imax; i < n_; i++)
            u_[i] = 0.;
    }
};

#endif
