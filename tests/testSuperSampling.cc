#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <array>
#include "SuperSampling.h"
using namespace std;


int unitTestSuperSamplingGaussian(void);
int unitTestSuperSamplingConstant(void);
int unitTestSuperSamplingGausHighFreq(void);
int unitTestSuperSamplingScale(void);
int unitTestSuperSamplingDelta(void);
int unitTestSuperSamplingHarmonics(void);
int unitTestSuperSamplingFiltering(void);

double gaussianFunc(double radius) { return exp(-radius*radius/2); }
double constantFunc(double radius) { return 1; }
double gaussianFuncScaled(double radius) { return exp(-radius*radius/200); }
double deltaFilter(double point) { return (abs(point)<1e-10 ? 3 : 0); }
double highFreq(double radius) { return exp(-radius*radius/2) + sin(100*radius); }

int main()
{
    int ret = 0;

    int deltTest = unitTestSuperSamplingDelta();
    if(deltTest > 0)
    {
        cout << "SuperSampling Delta Test failed!!!" << endl;
        ret = 1;
    }

    int freqTest = unitTestSuperSamplingGausHighFreq();
    if(freqTest > 0)
    {
        cout << "SuperSampling Frequency Test failed!!!" << endl;
        ret = 1;
    }

    int gausTest = unitTestSuperSamplingGaussian();
    if(gausTest > 0)
    {
        cout << "SuperSampling Gaussian Test failed!!!" << endl;
        ret = 1;
    }
    //int consTest = unitTestSuperSamplingConstant();
    //cout << "Constant Test: " << consTest << '\n';
    int scalTest = unitTestSuperSamplingScale();
    if(scalTest > 0)
    {
        cout << "SuperSampling Scale Test failed!!!" << endl;
        ret = 1;
    }

    int otherTest = unitTestSuperSamplingHarmonics();
    if(otherTest > 0)
    {
        cout << "SuperSampling Harmonics Test failed!!!" << endl;
        ret = 1;
    }
    int filterTest = unitTestSuperSamplingFiltering();
    if(filterTest > 0)
    {
        cout << "SuperSampling Filtering Test failed!!!" << endl;
        ret = 1;
    }

    return ret;
}


int unitTestSuperSamplingFiltering()
{
    std::cout << "unitTestSuperSamplingFiltering...\n";

    int result=0;
    
    //array<double,3> atomicCenter={.06,-.06,.1};
    array<double,3> atomicCenter={0,0,0};
    //double cutOffRadius=.5;
    int numExtraPts=10;
    int sampleRate=3;
    const int lMax=0;
    //array<double,3> coarGridSpace={.099,.099,.105};
    array<double,3> coarGridSpace={.1,.1,.1};
    array<double,3> botMeshCorner={-.1,-.1,-.1};
    array<double,3> topMeshCorner={.1,.1,.1};
    bool harmonics=false;
    
    SuperSampling<lMax>::setup(sampleRate, numExtraPts, coarGridSpace);
    SuperSampling<lMax> uniTest(atomicCenter, botMeshCorner, topMeshCorner, harmonics, gaussianFunc);
    
    int numPts=round((topMeshCorner[0]-botMeshCorner[0])/coarGridSpace[0]) + 1;
    int offset=0;
    double idx=botMeshCorner[0];
    for(int i=0; i<numPts; ++i) {
        double idy=botMeshCorner[1];
        for(int j=0; j<numPts; ++j) {
            double idz=botMeshCorner[2];
            for(int k=0; k<numPts; ++k) {
                double radius = sqrt((idx-atomicCenter[0])*(idx-atomicCenter[0]) + 
                                     (idy-atomicCenter[1])*(idy-atomicCenter[1]) + 
                                     (idz-atomicCenter[2])*(idz-atomicCenter[2]) );
                //cout << i << '\t' << j << '\t' << k << '\n';
                //cout << uniTest.values_[0][offset] - gaussianFunc(radius) << "\n\n";
                if(abs( uniTest.values_[0][offset] - gaussianFunc(radius)) > .1) {
                        result=1;
                        cout << i << '\t' << j << '\t' << k << '\n';
                        cout << uniTest.values_[0][offset] - gaussianFunc(radius) << "\n\n";
                        cout << uniTest.values_[0][offset] << '\n';
                        cout << gaussianFunc(radius) << '\n';
                        cout << idx << ',' << idy << ',' << idz << '\n';
                }
                offset += 1;
                idz += coarGridSpace[2];
            }
            idy += coarGridSpace[1];
        }
        idx += coarGridSpace[0];
    }
    return result;
}


int unitTestSuperSamplingHarmonics()
{
    std::cout << "unitTestSuperSamplingHarmonics...\n";
    int result=0;
    
    array<double,3> atomicCenter={0,0,0};
    int numExtraPts=10;
    int sampleRate=3;
    const int lMax=2;
    array<double,3> coarGridSpace={1,1,1};
    array<double,3> botMeshCorner={-10,-10,-10};
    array<double,3> topMeshCorner={10,10,10};
    bool harmonics=true;
    
    SuperSampling<lMax>::setup(sampleRate, numExtraPts, coarGridSpace);
    SuperSampling<lMax> uniTestScaled(atomicCenter, botMeshCorner, topMeshCorner, harmonics, gaussianFuncScaled);

    return result;
}

int unitTestSuperSamplingScale()
{
    std::cout << "unitTestSuperSamplingScale...\n";

    int result=0;
    
    array<double,3> atomicCenter={0,0,0};
    //double cutOffRadius=5;
    int numExtraPts=10;
    int sampleRate=3;
    const int lMax=0;
    array<double,3> coarGridSpace={1,1,1};
    array<double,3> botMeshCorner={-10,-10,-10};
    array<double,3> topMeshCorner={10,10,10};
    bool harmonics=false;
    
    SuperSampling<lMax>::setup(sampleRate, numExtraPts, coarGridSpace);
    SuperSampling<lMax> uniTestScaled(atomicCenter, botMeshCorner, topMeshCorner, harmonics, gaussianFuncScaled);
    
    coarGridSpace={.1,.1,.1};
    botMeshCorner={-1,-1,-1};
    topMeshCorner={1,1,1};
    
    SuperSampling<lMax>::setup(sampleRate, numExtraPts, coarGridSpace);
    SuperSampling<lMax> uniTestCompare(atomicCenter, botMeshCorner, topMeshCorner, harmonics, gaussianFunc);
    
    int numPts=round((topMeshCorner[0]-botMeshCorner[0])/coarGridSpace[0]) + 1;
    double idx=botMeshCorner[0];
    for(int i=0; i<numPts; ++i) {
       double idy=botMeshCorner[1];
       for(int j=0; j<numPts; ++j) {
          double idz=botMeshCorner[2];
          for(int k=0; k<numPts; ++k) {
              //cout << idx << ",\t" << idy << ",\t" << idz << '\n';
              //cout << 
              //    uniTestScaled.values_ [0][numPts*numPts*i + numPts*j + k] -
              //    uniTestCompare.values_[0][numPts*numPts*i + numPts*j + k] << "\n\n";
              if(abs( 
                  uniTestScaled.values_ [0][numPts*numPts*i + numPts*j + k] -
                  uniTestCompare.values_[0][numPts*numPts*i + numPts*j + k]
                  ) > 1e-10) {
                  result=1;
                  //cout << 
                  //    uniTestScaled.values_ [0][numPts*numPts*i + numPts*j + k] -
                  //    uniTestCompare.values_[0][numPts*numPts*i + numPts*j + k] << "\n\n";
              }
              idz += coarGridSpace[2];
          }
          idy += coarGridSpace[1];
       }
       idx += coarGridSpace[0];
    }
    return result;
}

int unitTestSuperSamplingGausHighFreq(void)
{
    std::cout << "unitTestSuperSamplingGausHighFreq...\n";
    int result=0;

    //array<double,3> atomicCenter={.601,-.65,.1};
    //array<double,3> atomicCenter={0,0,0};
    array<double,3> atomicCenter={0.001,0.0,-.001};
    //double cutOffRadius=.5;
    int numExtraPts=10;
    int sampleRate=3;
    const int lMax=0;
    //array<double,3> coarGridSpace={.099,.099,.105};
    array<double,3> coarGridSpace={.1,.1,.1};
    array<double,3> botMeshCorner={-1,-1,-1};
    array<double,3> topMeshCorner={1,1,1};
    bool harmonics=false;
    
    SuperSampling<lMax>::setup(sampleRate, numExtraPts, coarGridSpace);
    SuperSampling<lMax> uniTest(atomicCenter, botMeshCorner, topMeshCorner, harmonics, highFreq);
    
    int numPts=round((topMeshCorner[0]-botMeshCorner[0])/coarGridSpace[0]) + 1;
    int offset=0;
    double idx=botMeshCorner[0];
    ofstream myfile, myfile2;
    myfile.open("testSuperSamplingHighFreq.txt");
    myfile2.open("testSuperSamplingHighFreqAcutalFunc.txt");
    for(int i=0; i<numPts; ++i) {
        double idy=botMeshCorner[1];
        for(int j=0; j<numPts; ++j) {
            double idz=botMeshCorner[2];
            for(int k=0; k<numPts; ++k) {
                double radius = sqrt((idx-atomicCenter[0])*(idx-atomicCenter[0]) + 
                                     (idy-atomicCenter[1])*(idy-atomicCenter[1]) + 
                                     (idz-atomicCenter[2])*(idz-atomicCenter[2]) );
                //cout << i << '\t' << j << '\t' << k << '\n';
                //cout << uniTest.values_[0][offset] - gaussianFunc(radius) << "\n\n";
		        myfile << radius << ';' << uniTest.values_[0][offset] << ',';
		        myfile2<< radius << ';' << highFreq(radius) << ',';
                if(abs( uniTest.values_[0][offset] - gaussianFunc(radius)) > .25) {
                        result=1;
                        cout << i << '\t' << j << '\t' << k << '\n';
                        cout << uniTest.values_[0][offset] - gaussianFunc(radius) << "\n\n";
                }
                offset += 1;
                idz += coarGridSpace[2];
            }
            idy += coarGridSpace[1];
        }
        idx += coarGridSpace[0];
    }
    myfile.close();
    return result;
}

int unitTestSuperSamplingConstant()
{
    std::cout << "unitTestSuperSamplingConstant...\n";

    int result=0;
    
    array<double,3> atomicCenter={.6,-.6,.1};
    //array<double,3> atomicCenter={0,0,0};
    //double cutOffRadius=.5;
    int numExtraPts=40;
    int sampleRate=3;
    const int lMax=0;
    //array<double,3> coarGridSpace={.099,.099,.105};
    array<double,3> coarGridSpace={.1,.1,.1};
    array<double,3> botMeshCorner={-1,-1,-1};
    array<double,3> topMeshCorner={0,1,1};
    bool harmonics=false;
    
    SuperSampling<lMax>::setup(sampleRate, numExtraPts, coarGridSpace);
    SuperSampling<lMax> uniTest(atomicCenter, botMeshCorner, topMeshCorner, harmonics, constantFunc);
    
    int numPts=round((topMeshCorner[0]-botMeshCorner[0])/coarGridSpace[0]) + 1;
    int offset=0;
    double idx=botMeshCorner[0];
    for(int i=0; i<numPts; ++i) {
        double idy=botMeshCorner[1];
        for(int j=0; j<numPts; ++j) {
            double idz=botMeshCorner[2];
            for(int k=0; k<numPts; ++k) {
                double radius = sqrt((idx-atomicCenter[0])*(idx-atomicCenter[0]) + 
                                     (idy-atomicCenter[1])*(idy-atomicCenter[1]) + 
                                     (idz-atomicCenter[2])*(idz-atomicCenter[2]) );
                //cout << i << '\t' << j << '\t' << k << '\n';
                //cout << uniTest.values_[0][offset] - gaussianFunc(radius) << "\n\n";
                if(abs( uniTest.values_[0][offset] - constantFunc(radius)) > 1e-10) {
                        result=1;
                        cout << i << '\t' << j << '\t' << k << '\n';
                        cout << uniTest.values_[0][offset] - constantFunc(radius) << "\n\n";
                }
                offset += 1;
                idz += coarGridSpace[2];
            }
            idy += coarGridSpace[1];
        }
        idx += coarGridSpace[0];
    }
    return result;
}

int unitTestSuperSamplingGaussian()
{
    std::cout << "unitTestSuperSamplingGaussian...\n";

    int result=0;
    
    array<double,3> atomicCenter={.6,-.6,.1};
    //array<double,3> atomicCenter={0,0,0};
    //double cutOffRadius=.5;
    int numExtraPts=10;
    int sampleRate=3;
    const int lMax=0;
    //array<double,3> coarGridSpace={.099,.099,.105};
    array<double,3> coarGridSpace={.1,.1,.1};
    array<double,3> botMeshCorner={-1,-1,-1};
    array<double,3> topMeshCorner={0,1,1};
    bool harmonics=false;
    
    SuperSampling<lMax>::setup(sampleRate, numExtraPts, coarGridSpace);
    SuperSampling<lMax> uniTest(atomicCenter, botMeshCorner, topMeshCorner, harmonics, gaussianFunc);
    
    int numPtsX=round((topMeshCorner[0]-botMeshCorner[0])/coarGridSpace[0]) + 1;
    int numPtsY=round((topMeshCorner[1]-botMeshCorner[1])/coarGridSpace[1]) + 1;
    int numPtsZ=round((topMeshCorner[2]-botMeshCorner[2])/coarGridSpace[2]) + 1;
    int offset=0;
    double idx=botMeshCorner[0];
    for(int i=0; i<numPtsX; ++i) {
        double idy=botMeshCorner[1];
        for(int j=0; j<numPtsY; ++j) {
            double idz=botMeshCorner[2];
            for(int k=0; k<numPtsZ; ++k) {
                double radius = sqrt((idx-atomicCenter[0])*(idx-atomicCenter[0]) + 
                                     (idy-atomicCenter[1])*(idy-atomicCenter[1]) + 
                                     (idz-atomicCenter[2])*(idz-atomicCenter[2]) );
                //cout << i << '\t' << j << '\t' << k << '\n';
                //cout << uniTest.values_[0][offset] - gaussianFunc(radius) << "\n\n";
                if(abs( uniTest.values_[0][offset] - gaussianFunc(radius)) > 0.1) {
                        result=1;
                        cout << i << '\t' << j << '\t' << k << '\n';
                        cout << uniTest.values_[0][offset] - gaussianFunc(radius) << "\n\n";
                }
                offset += 1;
                idz += coarGridSpace[2];
            }
            idy += coarGridSpace[1];
        }
        idx += coarGridSpace[0];
    }
    return result;
}

int unitTestSuperSamplingDelta()
{
    std::cout << "unitTestSuperSamplingDelta...\n";

    int result=0;
    
    array<double,3> atomicCenter={0,0,0};
    //array<double,3> atomicCenter={.5,.45,-.43};
    //double cutOffRadius=.5;
    int numExtraPts=10;
    int sampleRate=3;
    const int lMax=0;
    array<double,3> coarGridSpace={.1,.1,.1};
    array<double,3> botMeshCorner={-1,-1,-1};
    array<double,3> topMeshCorner={1,1,1};
    bool harmonics=false;
    
    SuperSampling<lMax>::setup(sampleRate, numExtraPts, coarGridSpace, deltaFilter);
    SuperSampling<lMax> uniTest(atomicCenter, botMeshCorner, topMeshCorner, harmonics, gaussianFunc);
    
    int numPts=round((topMeshCorner[0]-botMeshCorner[0])/coarGridSpace[0]) + 1;
    double idx=botMeshCorner[0];
    for(int i=0; i<numPts; ++i) {
            double idy=botMeshCorner[1];
            for(int j=0; j<numPts; ++j) {
                    double idz=botMeshCorner[2];
                    for(int k=0; k<numPts; ++k) {
                            double radius = sqrt((idx-atomicCenter[0])*(idx-atomicCenter[0]) + 
                                                 (idy-atomicCenter[1])*(idy-atomicCenter[1]) + 
                                                 (idz-atomicCenter[2])*(idz-atomicCenter[2]) );
                            if(std::abs(
                         uniTest.values_[0][numPts*numPts*i + numPts*j + k] - gaussianFunc(radius)) > 1e-10) {
                                result=1;
                                //cout << idx << ",\t" << idy << ",\t" << idz << '\n';
                                //cout << uniTest.values_[0][numPts*numPts*i + numPts*j + k];
                                //cout << gaussianFunc(radius)/(2*sqrt(pi)) << "\n\n";
                            }
                            idz += coarGridSpace[2];
                    }
                    idy += coarGridSpace[1];
            }
            idx += coarGridSpace[0];
    }
    return result;
}
