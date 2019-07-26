#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <functional>

#include "RadialInter.h"

const double pi=3.14159265358979323846;

double sincFunction(double point);

double setSProjector(
        const double currentValue,
        const double radius);

std::array<double,3> setPProjector(
        const double currentValue,
        const double radius,
        const double constant,
        const std::array<double,3> xyz);

std::array<double,5> setDProjector(
        const double currentValue,
        const double radius,
        const std::array<double,3> constants,
        const std::array<double,3> xyz,
        const std::array<double,3> xyz2);

std::array<double,7> setFProjector(
        const double currentValue,
        const double radius,
        const std::array<double,5> constants,
        const std::array<double,3> xyz,
        const std::array<double,3> xyz2);

std::vector<double> computeWeights(
        const int numExtraPts,
        const std::array<double,3> fineGridSpace,
        const std::array<double,3> coarGridSpace,
        double filter(double));

template <int lMax>
class superSampling {
    const bool harmonics_;
    const std::array<double,3> atomicCenter_;
    const std::array<double,3> botMeshCorner_;
    const std::array<double,3> topMeshCorner_;
    int convNumPts_;
    std::array<int,3> coarNumPts_;
    std::array<int,3> fineNumPts_;
    std::array<int,3> funcNumPts_;
    
    static int sampleRate_;
    static int numExtraPts_;
    static std::array<double,3> coarGridSpace_;
    static std::array<double,3> fineGridSpace_;
    static std::vector<double> weights_;
    
    // Constants
    
    const double constLMAX0_ = 1.0 / (2.0*std::sqrt(pi));
    const double constLMAX1_ = std::sqrt(3.0 / (4.0*pi));
    const std::array<double,3> constLMAX2_ = {
            std::sqrt(15.0/pi) / 2.0,
            std::sqrt(5.0 /pi) / 2.0,
            std::sqrt(15.0/pi) / 4.0};
    const std::array<double,5> constLMAX3_ = {
            std::sqrt(35.0 /(2.0*pi)) / 4.0,
            std::sqrt(105.0/     pi ) / 2.0,
            std::sqrt(21.0 /(2.0*pi)) / 4.0,
            std::sqrt(7.0  /(    pi)) / 4.0,
            std::sqrt(105.0/     pi ) / 4.0};
public :
    std::array<std::vector<double>, 2*lMax + 1> values_;

    static void setup(
             int sampleRate, 
             int numExtraPts,
             std::array<double,3> coarGridSpace, 
             double filter(double)=sincFunction
         );
    
    superSampling(
                     std::array<double,3> atomicCenter, 
                     std::array<double,3> botMeshCorner,
                     std::array<double,3> topMeshCorner,
                     bool harmonics,
                     std::function<double(double)> const& Func
                     //const RadialInter& objectFunc 
                 ) : 
                     atomicCenter_(atomicCenter), botMeshCorner_(botMeshCorner),
                     topMeshCorner_(topMeshCorner), harmonics_(harmonics) {
        convNumPts_ = 2 * numExtraPts_ + 1;
        coarNumPts_[0] = std::round((topMeshCorner_[0] - botMeshCorner_[0])/coarGridSpace_[0]) + 1;
        fineNumPts_[0] = std::round((topMeshCorner_[0] - botMeshCorner_[0])/fineGridSpace_[0]) + 1;
        funcNumPts_[0] = convNumPts_ + fineNumPts_[0] - 1;
        coarNumPts_[1] = std::round((topMeshCorner_[1] - botMeshCorner_[1])/coarGridSpace_[1]) + 1;
        fineNumPts_[1] = std::round((topMeshCorner_[1] - botMeshCorner_[1])/fineGridSpace_[1]) + 1;
        funcNumPts_[1] = convNumPts_ + fineNumPts_[1] - 1;
        coarNumPts_[2] = std::round((topMeshCorner_[2] - botMeshCorner_[2])/coarGridSpace_[2]) + 1;
        fineNumPts_[2] = std::round((topMeshCorner_[2] - botMeshCorner_[2])/fineGridSpace_[2]) + 1;
        funcNumPts_[2] = convNumPts_ + fineNumPts_[2] - 1;
        //std::array<std::vector<double>,2*lMax + 1> fineMeshFuncValues = getFuncValues(objectFunc);
        std::array<std::vector<double>,2*lMax + 1> fineMeshFuncValues = getFuncValues(Func);
        for(int ll=0; ll<2*lMax + 1; ++ll) {
            values_[ll] = computeSuperSampling(fineMeshFuncValues[ll]);
        }
    }
    
private :
    
    std::vector<double> computeSuperSampling(
                    const std::vector<double>& fineMeshFuncValues
                    );

    std::array<std::vector<double>,2*lMax + 1> getFuncValues(std::function<double(const double)> Func);
    //std::array<std::vector<double>,2*lMax + 1> getFuncValues(const RadialInter& objectFunc);
    
    double computeFiltering(
            const std::array<int,3> sampleIndex,
            const std::vector<double>& fineMeshFuncValues
            );
};

template class superSampling<0>;
template class superSampling<1>;
template class superSampling<2>;
template class superSampling<3>;


template <int lMax>
int superSampling<lMax>::sampleRate_=0;
template <int lMax>
int superSampling<lMax>::numExtraPts_=0;
template <int lMax>
std::array<double,3> superSampling<lMax>::coarGridSpace_={0,0,0};
template <int lMax>
std::array<double,3> superSampling<lMax>::fineGridSpace_={0,0,0};
template <int lMax>
std::vector<double> superSampling<lMax>::weights_;
