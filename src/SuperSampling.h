#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

const double pi = 3.14159265358979323846;

double sincFunction(double point);

double setSProjector(const double currentValue, const double radius);

std::array<double, 3> setPProjector(const double currentValue,
    const double radius, const double constant,
    const std::array<double, 3> xyz);

std::array<double, 5> setDProjector(const double currentValue,
    const double radius, const std::array<double, 3> constants,
    const std::array<double, 3> xyz, const std::array<double, 3> xyz2);

std::array<double, 7> setFProjector(const double currentValue,
    const double radius, const std::array<double, 5> constants,
    const std::array<double, 3> xyz, const std::array<double, 3> xyz2);

std::vector<double> computeWeights(const int numExtraPts,
    const std::array<double, 3> fineGridSpace,
    const std::array<double, 3> coarGridSpace, double filter(double));

template <int lMax>
class SuperSampling
{
    const bool harmonics_;
    const std::array<double, 3> atomicCenter_;
    const std::array<double, 3> botMeshCorner_;
    const std::array<double, 3> topMeshCorner_;
    int convNumPts_;
    std::array<int, 3> coarNumPts_;
    std::array<int, 3> fineNumPts_;
    std::array<int, 3> funcNumPts_;

    static int sampleRate_;
    static int numExtraPts_;
    static std::array<double, 3> coarGridSpace_;
    static std::array<double, 3> fineGridSpace_;
    static std::vector<double> weights_;

    // Constants

    const double constLMAX0_                = 1.0 / (2.0 * std::sqrt(pi));
    const double constLMAX1_                = std::sqrt(3.0 / (4.0 * pi));
    const std::array<double, 3> constLMAX2_ = { std::sqrt(15.0 / pi) / 2.0,
        std::sqrt(5.0 / pi) / 2.0, std::sqrt(15.0 / pi) / 4.0 };
    const std::array<double, 5> constLMAX3_
        = { std::sqrt(35.0 / (2.0 * pi)) / 4.0, std::sqrt(105.0 / pi) / 2.0,
              std::sqrt(21.0 / (2.0 * pi)) / 4.0, std::sqrt(7.0 / (pi)) / 4.0,
              std::sqrt(105.0 / pi) / 4.0 };

public:
    std::array<std::vector<double>, 2 * lMax + 1> values_;

    // setup must always have been called before SuperSampling, but only needs
    // to be called once. Note that using different lMax values means that setup
    // needs to be called for each of them.
    static void setup(int sampleRate, int numExtraPts,
        std::array<double, 3> coarGridSpace,
        double filter(double) = sincFunction);

    // Constructor does the supersampling. top and bot MeshCorner allow us to
    // construct the subdomain mesh to supersample, atomicCenter is obvious,
    // harmonics needs to be true if spherical harmonics is needed. If no
    // harmonics, then should always call SuperSampling<0>. If harmonics, then
    // lmax value is called in constructor SuperSampling<lMax>
    SuperSampling(std::array<double, 3> atomicCenter,
        std::array<double, 3> botMeshCorner,
        std::array<double, 3> topMeshCorner, bool harmonics,
        std::function<double(double)> const& Func)
        : harmonics_(harmonics),
          atomicCenter_(atomicCenter),
          botMeshCorner_(botMeshCorner),
          topMeshCorner_(topMeshCorner)
    {
        convNumPts_    = 2 * numExtraPts_ + 1;
        coarNumPts_[0] = std::round((topMeshCorner_[0] - botMeshCorner_[0])
                                    / coarGridSpace_[0])
                         + 1;
        fineNumPts_[0] = std::round((topMeshCorner_[0] - botMeshCorner_[0])
                                    / fineGridSpace_[0])
                         + 1;
        funcNumPts_[0] = convNumPts_ + fineNumPts_[0] - 1;
        coarNumPts_[1] = std::round((topMeshCorner_[1] - botMeshCorner_[1])
                                    / coarGridSpace_[1])
                         + 1;
        fineNumPts_[1] = std::round((topMeshCorner_[1] - botMeshCorner_[1])
                                    / fineGridSpace_[1])
                         + 1;
        funcNumPts_[1] = convNumPts_ + fineNumPts_[1] - 1;
        coarNumPts_[2] = std::round((topMeshCorner_[2] - botMeshCorner_[2])
                                    / coarGridSpace_[2])
                         + 1;
        fineNumPts_[2] = std::round((topMeshCorner_[2] - botMeshCorner_[2])
                                    / fineGridSpace_[2])
                         + 1;
        funcNumPts_[2] = convNumPts_ + fineNumPts_[2] - 1;
        std::array<std::vector<double>, 2 * lMax + 1> fineMeshFuncValues
            = getFuncValues(Func);
        for (int ll = 0; ll < 2 * lMax + 1; ++ll)
        {
            values_[ll] = computeSuperSampling(fineMeshFuncValues[ll]);
        }
    }

private:
    // Loop through coarse subdomain mesh and compute convolution at each of
    // those points
    std::vector<double> computeSuperSampling(
        const std::vector<double>& fineMeshFuncValues);

    // Store all the function values to be used in the algorithm in a single
    // vector Accessing values in vector is quicker than evaluating a function
    // every time
    std::array<std::vector<double>, 2 * lMax + 1> getFuncValues(
        const std::function<double(const double)>& Func);

    // At a given point in the coarse subdomain mesh, loop through the nearby
    // points to compute the convolution
    double computeFiltering(const std::array<int, 3> sampleIndex,
        const std::vector<double>& fineMeshFuncValues);
};

template <int lMax>
int SuperSampling<lMax>::sampleRate_ = 0;
template <int lMax>
int SuperSampling<lMax>::numExtraPts_ = 0;
template <int lMax>
std::array<double, 3> SuperSampling<lMax>::coarGridSpace_ = { 0, 0, 0 };
template <int lMax>
std::array<double, 3> SuperSampling<lMax>::fineGridSpace_ = { 0, 0, 0 };
template <int lMax>
std::vector<double> SuperSampling<lMax>::weights_;
