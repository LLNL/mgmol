#include "ChebyshevApproximationInterface.h"

Timer ChebyshevApproximationInterface::compute_coeffs_tm_(
    "ChebyshevApproximationInterface::compChebCoeffs");
Timer ChebyshevApproximationInterface::double_loop_tm_(
    "ChebyshevApproximationInterface::double_loop");

void ChebyshevApproximationInterface::evaluateFunction(
    const std::vector<double>& fnodes, std::vector<double>& fvals)
{
    fvals.clear();
    fvals = chebfunc_->eval(fnodes);
}

void ChebyshevApproximationInterface::computeChebyshevCoeffs()
{
    compute_coeffs_tm_.start();
    order_ = max_order_;
    int n  = order_;
    // define scaling variable here
    const double fac = 2.0 / (double)n;
    // evaluate function at interpolation points
    std::vector<double> fvals;
    evaluateFunction(interp_points_, fvals);

    // build coefficients
    coeffs_.clear();
    coeffs_.reserve(order_);

    // initialize first entry of coeffs_
    double fsum = 0.;
    for (std::vector<double>::iterator it = fvals.begin(); it != fvals.end();
         ++it)
        fsum += *it;
    // scale and store first coeff
    fsum /= (double)n;
    coeffs_.push_back(fsum);

    // compute remaining coeffs_ values
    double_loop_tm_.start();

    LocalVector<double> fvec(fvals);
    coeffs_.resize(n, 0.);
    LocalVector<double> cvec(coeffs_);

    cmat_->matvec(fvec, cvec);
    (cvec.data())[0] = fsum;
    cvec.swap(coeffs_);

    double_loop_tm_.stop();

    // scale only last n-1 coefficients
    const int n1 = n - 1;
    MPscal(n1, fac, &coeffs_[1]);

    //    for (int i = 1; i < n; i++)assert(coeffs_[i] == coeffs_[i]);

    compute_coeffs_tm_.stop();
}

// Compute Chebyshev approximation given function points.
// Assume that points are sorted from left to right on the real axis.
// Example: (sorted) eigenvalues of some matrix.
void ChebyshevApproximationInterface::computeChebyshevApproximation(
    std::vector<double>& points, std::vector<double>& vals)
{
    assert((int)points.size() == order_);
    vals.clear();
    // get extents of points
    const double a = points[0];
    const double b = points[order_ - 1];
    // scale points to be in the range [-1, 1]
    std::vector<double> scaled_points;
    scalePointsToChebyshevInterval(points, scaled_points);

    for (int j = 0; j < order_; j++)
        assert(fabs(scaled_points[j]) <= 1.);

    // compute Chebyshev coefficients
    computeChebyshevCoeffs();

    // initialize vals
    for (int j = 0; j < order_; j++)
        vals.push_back(coeffs_[0]);

    // loop to compute approximation values
    for (int k = 1; k < order_; k++)
    {
        for (int j = 0; j < order_; j++)
        {
            double prod = k * acos(scaled_points[j]);
            double t_kj = cos(prod);
            vals[j] += coeffs_[k] * t_kj;
        }
    }

    for (int j = 0; j < order_; j++)
        assert(vals[j] == vals[j]);
}

void ChebyshevApproximationInterface::computeInterpolationPoints()
{
    assert(order_ > 0);

    interp_points_.clear();
    angles_.clear();
    const double a  = extents_[0];
    const double b  = extents_[1];
    const double a1 = -1.;
    const double b1 = 1.;

    const double alpha = M_PI / (2. * order_);
    const double beta  = (b - a) / (b1 - a1);
    for (int i = 0; i < order_; i++)
    {
        double ang = (2. * (i + 1.) - 1.) * alpha;
        angles_.push_back(ang);
        double x_i = cos(ang);
        // map node to interval [a, b]
        const double scaled_node = a + beta * (x_i - a1);

        assert(scaled_node == scaled_node);
        interp_points_.push_back(scaled_node);
    }

    double_loop_tm_.start();
    if (cmat_) delete cmat_;
    cmat_ = new LocalMatrices<double>(1, order_, order_);
    for (int i = 1; i < order_; i++)
    {
        for (int k = 0; k < order_; k++)
        {
            double iang = i * angles_[k];
            double val  = cos(iang);
            cmat_->setVal(i, k, val);
        }
    }
    double_loop_tm_.stop();
}

void ChebyshevApproximationInterface::scalePointsToChebyshevInterval(
    const std::vector<double>& points, std::vector<double>& scaled_points)
{
    const double a1 = -1;
    const double b1 = 1;
    // get extents of points
    const double a = points[0];
    const double b = points[order_ - 1];

    const double alpha = (b1 - a1) / (b - a);

    scaled_points.clear();
    scaled_points.reserve(points.size());

    for (std::vector<double>::const_iterator it = points.begin();
         it != points.end(); ++it)
    {
        double sp = a1 + alpha * (*it - a);
        assert(sp == sp);

        scaled_points.push_back(sp);
    }
}
