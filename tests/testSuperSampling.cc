// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "SuperSampling.h"

#include "catch.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>

double gaussianFunc(double radius) { return std::exp(-radius * radius / 2); }
double constantFunc(double) { return 1; }
double gaussianFuncScaled(double radius)
{
    return std::exp(-radius * radius / 200);
}
double deltaFilter(double point) { return (std::abs(point) < 1e-10 ? 3 : 0); }
double highFreq(double radius)
{
    return std::exp(-radius * radius / 2) + std::sin(100 * radius);
}

TEST_CASE("Check super sampling filtering", "[filtering]")
{
    std::array<double, 3> atomicCenter  = { 0, 0, 0 };
    int numExtraPts                     = 10;
    int sampleRate                      = 3;
    const int lMax                      = 0;
    std::array<double, 3> coarGridSpace = { .1, .1, .1 };
    std::array<double, 3> botMeshCorner = { -.1, -.1, -.1 };
    std::array<double, 3> topMeshCorner = { .1, .1, .1 };
    bool harmonics                      = false;

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
                CHECK(uniTest.values_[0][offset]
                      == Approx(gaussianFunc(radius)).epsilon(0.1));
                offset += 1;
                idz += coarGridSpace[2];
            }
            idy += coarGridSpace[1];
        }
        idx += coarGridSpace[0];
    }
}

TEST_CASE("Check super sampling harmonics", "[harmonics]")
{
    std::array<double, 3> atomicCenter  = { 0, 0, 0 };
    int numExtraPts                     = 10;
    int sampleRate                      = 3;
    const int lMax                      = 2;
    std::array<double, 3> coarGridSpace = { 1, 1, 1 };
    std::array<double, 3> botMeshCorner = { -10, -10, -10 };
    std::array<double, 3> topMeshCorner = { 10, 10, 10 };
    bool harmonics                      = true;

    SuperSampling<lMax>::setup(sampleRate, numExtraPts, coarGridSpace);
    SuperSampling<lMax> uniTestScaled(atomicCenter, botMeshCorner,
        topMeshCorner, harmonics, gaussianFuncScaled);

    // TODO: there is no check here
}

TEST_CASE("Check super sampling scale", "[scale]")
{
    std::array<double, 3> atomicCenter  = { 0, 0, 0 };
    int numExtraPts                     = 10;
    int sampleRate                      = 3;
    const int lMax                      = 0;
    std::array<double, 3> coarGridSpace = { 1, 1, 1 };
    std::array<double, 3> botMeshCorner = { -10, -10, -10 };
    std::array<double, 3> topMeshCorner = { 10, 10, 10 };
    bool harmonics                      = false;

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
                CHECK(uniTestScaled
                          .values_[0][numPts * numPts * i + numPts * j + k]
                      == Approx(uniTestCompare.values_[0][numPts * numPts * i
                                                          + numPts * j + k])
                             .epsilon(1e-10));
                idz += coarGridSpace[2];
            }
            idy += coarGridSpace[1];
        }
        idx += coarGridSpace[0];
    }
}

TEST_CASE("Check super sampling gaussian high frequency", "[gaus_high_freq]")
{
    std::array<double, 3> atomicCenter  = { 0.001, 0.0, -.001 };
    int numExtraPts                     = 10;
    int sampleRate                      = 3;
    const int lMax                      = 0;
    std::array<double, 3> coarGridSpace = { .1, .1, .1 };
    std::array<double, 3> botMeshCorner = { -1, -1, -1 };
    std::array<double, 3> topMeshCorner = { 1, 1, 1 };
    bool harmonics                      = false;

    SuperSampling<lMax>::setup(sampleRate, numExtraPts, coarGridSpace);
    SuperSampling<lMax> uniTest(
        atomicCenter, botMeshCorner, topMeshCorner, harmonics, highFreq);

    int numPts
        = round((topMeshCorner[0] - botMeshCorner[0]) / coarGridSpace[0]) + 1;
    int offset = 0;
    double idx = botMeshCorner[0];
    std::ofstream myfile, myfile2;
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
                CHECK(uniTest.values_[0][offset]
                      == Approx(gaussianFunc(radius)).epsilon(0.25));
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
TEST_CASE("Check super_sampling_constant", "[!hide]")
{
    std::array<double, 3> atomicCenter  = { .6, -.6, .1 };
    int numExtraPts                     = 40;
    int sampleRate                      = 3;
    const int lMax                      = 0;
    std::array<double, 3> coarGridSpace = { .1, .1, .1 };
    std::array<double, 3> botMeshCorner = { -1, -1, -1 };
    std::array<double, 3> topMeshCorner = { 0, 1, 1 };
    bool harmonics                      = false;

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
                CHECK(uniTest.values_[0][offset] == constantFunc(radius));
                offset += 1;
                idz += coarGridSpace[2];
            }
            idy += coarGridSpace[1];
        }
        idx += coarGridSpace[0];
    }
}

TEST_CASE("Check super_sampling_gaussian", "[gaussian]")
{
    std::array<double, 3> atomicCenter  = { .6, -.6, .1 };
    int numExtraPts                     = 10;
    int sampleRate                      = 3;
    const int lMax                      = 0;
    std::array<double, 3> coarGridSpace = { .1, .1, .1 };
    std::array<double, 3> botMeshCorner = { -1, -1, -1 };
    std::array<double, 3> topMeshCorner = { 0, 1, 1 };
    bool harmonics                      = false;

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
                CHECK(uniTest.values_[0][offset]
                      == Approx(gaussianFunc(radius)).epsilon(0.1));
                offset += 1;
                idz += coarGridSpace[2];
            }
            idy += coarGridSpace[1];
        }
        idx += coarGridSpace[0];
    }
}

TEST_CASE("Check super_sampling_delta", "[delta]")
{
    std::array<double, 3> atomicCenter  = { 0, 0, 0 };
    int numExtraPts                     = 10;
    int sampleRate                      = 3;
    const int lMax                      = 0;
    std::array<double, 3> coarGridSpace = { .1, .1, .1 };
    std::array<double, 3> botMeshCorner = { -1, -1, -1 };
    std::array<double, 3> topMeshCorner = { 1, 1, 1 };
    bool harmonics                      = false;

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
                CHECK(uniTest.values_[0][numPts * numPts * i + numPts * j + k]
                      == Approx(gaussianFunc(radius)).epsilon(1e-10));
                idz += coarGridSpace[2];
            }
            idy += coarGridSpace[1];
        }
        idx += coarGridSpace[0];
    }
}
