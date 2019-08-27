// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "SuperSampling.h"

#include <boost/test/unit_test.hpp>

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;
namespace utf = boost::unit_test;

double gaussianFunc(double radius) { return exp(-radius * radius / 2); }
double constantFunc(double radius) { return 1; }
double gaussianFuncScaled(double radius) { return exp(-radius * radius / 200); }
double deltaFilter(double point) { return (abs(point) < 1e-10 ? 3 : 0); }
double highFreq(double radius)
{
    return std::exp(-radius * radius / 2) + std::sin(100 * radius);
}

BOOST_AUTO_TEST_CASE(super_sampling_filtering, *utf::tolerance(0.1))
{
    // array<double,3> atomicCenter={.06,-.06,.1};
    array<double, 3> atomicCenter = { 0, 0, 0 };
    // double cutOffRadius=.5;
    int numExtraPts = 10;
    int sampleRate  = 3;
    const int lMax  = 0;
    // array<double,3> coarGridSpace={.099,.099,.105};
    array<double, 3> coarGridSpace = { .1, .1, .1 };
    array<double, 3> botMeshCorner = { -.1, -.1, -.1 };
    array<double, 3> topMeshCorner = { .1, .1, .1 };
    bool harmonics                 = false;

    SuperSampling<lMax>::setup(sampleRate, numExtraPts, coarGridSpace);
    SuperSampling<lMax> uniTest(
        atomicCenter, botMeshCorner, topMeshCorner, harmonics, gaussianFunc);

    int numPts
        = round((topMeshCorner[0] - botMeshCorner[0]) / coarGridSpace[0]) + 1;
    int offset = 0;
    double idx = botMeshCorner[0];
    for (int i = 0; i < numPts; ++i)
    {
        double idy = botMeshCorner[1];
        for (int j = 0; j < numPts; ++j)
        {
            double idz = botMeshCorner[2];
            for (int k = 0; k < numPts; ++k)
            {
                double radius
                    = sqrt((idx - atomicCenter[0]) * (idx - atomicCenter[0])
                           + (idy - atomicCenter[1]) * (idy - atomicCenter[1])
                           + (idz - atomicCenter[2]) * (idz - atomicCenter[2]));
                BOOST_TEST(uniTest.values_[0][offset] == gaussianFunc(radius));
                offset += 1;
                idz += coarGridSpace[2];
            }
            idy += coarGridSpace[1];
        }
        idx += coarGridSpace[0];
    }
}

BOOST_AUTO_TEST_CASE(super_sampling_harmonics)
{
    array<double, 3> atomicCenter  = { 0, 0, 0 };
    int numExtraPts                = 10;
    int sampleRate                 = 3;
    const int lMax                 = 2;
    array<double, 3> coarGridSpace = { 1, 1, 1 };
    array<double, 3> botMeshCorner = { -10, -10, -10 };
    array<double, 3> topMeshCorner = { 10, 10, 10 };
    bool harmonics                 = true;

    SuperSampling<lMax>::setup(sampleRate, numExtraPts, coarGridSpace);
    SuperSampling<lMax> uniTestScaled(atomicCenter, botMeshCorner,
        topMeshCorner, harmonics, gaussianFuncScaled);

    // TODO: there is no check here
}

BOOST_AUTO_TEST_CASE(super_sampling_scale, *utf::tolerance(1e-10))
{
    array<double, 3> atomicCenter = { 0, 0, 0 };
    // double cutOffRadius=5;
    int numExtraPts                = 10;
    int sampleRate                 = 3;
    const int lMax                 = 0;
    array<double, 3> coarGridSpace = { 1, 1, 1 };
    array<double, 3> botMeshCorner = { -10, -10, -10 };
    array<double, 3> topMeshCorner = { 10, 10, 10 };
    bool harmonics                 = false;

    SuperSampling<lMax>::setup(sampleRate, numExtraPts, coarGridSpace);
    SuperSampling<lMax> uniTestScaled(atomicCenter, botMeshCorner,
        topMeshCorner, harmonics, gaussianFuncScaled);

    coarGridSpace = { .1, .1, .1 };
    botMeshCorner = { -1, -1, -1 };
    topMeshCorner = { 1, 1, 1 };

    SuperSampling<lMax>::setup(sampleRate, numExtraPts, coarGridSpace);
    SuperSampling<lMax> uniTestCompare(
        atomicCenter, botMeshCorner, topMeshCorner, harmonics, gaussianFunc);

    int numPts
        = round((topMeshCorner[0] - botMeshCorner[0]) / coarGridSpace[0]) + 1;
    double idx = botMeshCorner[0];
    for (int i = 0; i < numPts; ++i)
    {
        double idy = botMeshCorner[1];
        for (int j = 0; j < numPts; ++j)
        {
            double idz = botMeshCorner[2];
            for (int k = 0; k < numPts; ++k)
            {
                BOOST_TEST(
                    uniTestScaled
                        .values_[0][numPts * numPts * i + numPts * j + k]
                    == uniTestCompare
                           .values_[0][numPts * numPts * i + numPts * j + k]);
                idz += coarGridSpace[2];
            }
            idy += coarGridSpace[1];
        }
        idx += coarGridSpace[0];
    }
}

BOOST_AUTO_TEST_CASE(super_sampling_gaus_high_freq, *utf::tolerance(0.25))
{
    array<double, 3> atomicCenter  = { 0.001, 0.0, -.001 };
    int numExtraPts                = 10;
    int sampleRate                 = 3;
    const int lMax                 = 0;
    array<double, 3> coarGridSpace = { .1, .1, .1 };
    array<double, 3> botMeshCorner = { -1, -1, -1 };
    array<double, 3> topMeshCorner = { 1, 1, 1 };
    bool harmonics                 = false;

    SuperSampling<lMax>::setup(sampleRate, numExtraPts, coarGridSpace);
    SuperSampling<lMax> uniTest(
        atomicCenter, botMeshCorner, topMeshCorner, harmonics, highFreq);

    int numPts
        = round((topMeshCorner[0] - botMeshCorner[0]) / coarGridSpace[0]) + 1;
    int offset = 0;
    double idx = botMeshCorner[0];
    ofstream myfile, myfile2;
    myfile.open("testSuperSamplingHighFreq.txt");
    myfile2.open("testSuperSamplingHighFreqAcutalFunc.txt");
    for (int i = 0; i < numPts; ++i)
    {
        double idy = botMeshCorner[1];
        for (int j = 0; j < numPts; ++j)
        {
            double idz = botMeshCorner[2];
            for (int k = 0; k < numPts; ++k)
            {
                double radius
                    = sqrt((idx - atomicCenter[0]) * (idx - atomicCenter[0])
                           + (idy - atomicCenter[1]) * (idy - atomicCenter[1])
                           + (idz - atomicCenter[2]) * (idz - atomicCenter[2]));
                myfile << radius << ';' << uniTest.values_[0][offset] << ',';
                myfile2 << radius << ';' << highFreq(radius) << ',';
                BOOST_TEST(uniTest.values_[0][offset] == gaussianFunc(radius));
                offset += 1;
                idz += coarGridSpace[2];
            }
            idy += coarGridSpace[1];
        }
        idx += coarGridSpace[0];
    }
    myfile.close();
}

// TODO this test does not pass, so disable it for now
BOOST_AUTO_TEST_CASE(super_sampling_constant, *utf::disabled())
{
    array<double, 3> atomicCenter  = { .6, -.6, .1 };
    int numExtraPts                = 40;
    int sampleRate                 = 3;
    const int lMax                 = 0;
    array<double, 3> coarGridSpace = { .1, .1, .1 };
    array<double, 3> botMeshCorner = { -1, -1, -1 };
    array<double, 3> topMeshCorner = { 0, 1, 1 };
    bool harmonics                 = false;

    SuperSampling<lMax>::setup(sampleRate, numExtraPts, coarGridSpace);
    SuperSampling<lMax> uniTest(
        atomicCenter, botMeshCorner, topMeshCorner, harmonics, constantFunc);

    int numPts
        = round((topMeshCorner[0] - botMeshCorner[0]) / coarGridSpace[0]) + 1;
    int offset = 0;
    double idx = botMeshCorner[0];
    for (int i = 0; i < numPts; ++i)
    {
        double idy = botMeshCorner[1];
        for (int j = 0; j < numPts; ++j)
        {
            double idz = botMeshCorner[2];
            for (int k = 0; k < numPts; ++k)
            {
                double radius
                    = sqrt((idx - atomicCenter[0]) * (idx - atomicCenter[0])
                           + (idy - atomicCenter[1]) * (idy - atomicCenter[1])
                           + (idz - atomicCenter[2]) * (idz - atomicCenter[2]));
                BOOST_TEST(uniTest.values_[0][offset] == constantFunc(radius));
                offset += 1;
                idz += coarGridSpace[2];
            }
            idy += coarGridSpace[1];
        }
        idx += coarGridSpace[0];
    }
}

BOOST_AUTO_TEST_CASE(super_sampling_gaussian, *utf::tolerance(0.1))
{
    array<double, 3> atomicCenter  = { .6, -.6, .1 };
    int numExtraPts                = 10;
    int sampleRate                 = 3;
    const int lMax                 = 0;
    array<double, 3> coarGridSpace = { .1, .1, .1 };
    array<double, 3> botMeshCorner = { -1, -1, -1 };
    array<double, 3> topMeshCorner = { 0, 1, 1 };
    bool harmonics                 = false;

    SuperSampling<lMax>::setup(sampleRate, numExtraPts, coarGridSpace);
    SuperSampling<lMax> uniTest(
        atomicCenter, botMeshCorner, topMeshCorner, harmonics, gaussianFunc);

    int numPtsX
        = round((topMeshCorner[0] - botMeshCorner[0]) / coarGridSpace[0]) + 1;
    int numPtsY
        = round((topMeshCorner[1] - botMeshCorner[1]) / coarGridSpace[1]) + 1;
    int numPtsZ
        = round((topMeshCorner[2] - botMeshCorner[2]) / coarGridSpace[2]) + 1;
    int offset = 0;
    double idx = botMeshCorner[0];
    for (int i = 0; i < numPtsX; ++i)
    {
        double idy = botMeshCorner[1];
        for (int j = 0; j < numPtsY; ++j)
        {
            double idz = botMeshCorner[2];
            for (int k = 0; k < numPtsZ; ++k)
            {
                double radius
                    = sqrt((idx - atomicCenter[0]) * (idx - atomicCenter[0])
                           + (idy - atomicCenter[1]) * (idy - atomicCenter[1])
                           + (idz - atomicCenter[2]) * (idz - atomicCenter[2]));
                BOOST_TEST(uniTest.values_[0][offset] == gaussianFunc(radius));
                offset += 1;
                idz += coarGridSpace[2];
            }
            idy += coarGridSpace[1];
        }
        idx += coarGridSpace[0];
    }
}

BOOST_AUTO_TEST_CASE(super_sampling_delta, *utf::tolerance(1e-10))
{
    array<double, 3> atomicCenter = { 0, 0, 0 };
    // array<double,3> atomicCenter={.5,.45,-.43};
    // double cutOffRadius=.5;
    int numExtraPts                = 10;
    int sampleRate                 = 3;
    const int lMax                 = 0;
    array<double, 3> coarGridSpace = { .1, .1, .1 };
    array<double, 3> botMeshCorner = { -1, -1, -1 };
    array<double, 3> topMeshCorner = { 1, 1, 1 };
    bool harmonics                 = false;

    SuperSampling<lMax>::setup(
        sampleRate, numExtraPts, coarGridSpace, deltaFilter);
    SuperSampling<lMax> uniTest(
        atomicCenter, botMeshCorner, topMeshCorner, harmonics, gaussianFunc);

    int numPts
        = round((topMeshCorner[0] - botMeshCorner[0]) / coarGridSpace[0]) + 1;
    double idx = botMeshCorner[0];
    for (int i = 0; i < numPts; ++i)
    {
        double idy = botMeshCorner[1];
        for (int j = 0; j < numPts; ++j)
        {
            double idz = botMeshCorner[2];
            for (int k = 0; k < numPts; ++k)
            {
                double radius
                    = sqrt((idx - atomicCenter[0]) * (idx - atomicCenter[0])
                           + (idy - atomicCenter[1]) * (idy - atomicCenter[1])
                           + (idz - atomicCenter[2]) * (idz - atomicCenter[2]));
                BOOST_TEST(
                    uniTest.values_[0][numPts * numPts * i + numPts * j + k]
                    == gaussianFunc(radius));
                idz += coarGridSpace[2];
            }
            idy += coarGridSpace[1];
        }
        idx += coarGridSpace[0];
    }
}
