// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef SAVEDATA_H_
#define SAVEDATA_H_

#ifdef SCALAPACK
#include "BlacsContext.h"
#include "DistMatrix.h"
#include "MatricesBlacsContext.h"
#endif

#include <csignal>
#include <fstream>
#include <iostream>
#include <iterator>
#include <mpi.h>
#include <string>

template <class T>
void saveData(std::vector<std::vector<T>> data, const char* filename)
{

    if (onpe0)
    {

        std::ofstream outputFile(filename);

        outputFile.precision(16);

        if (outputFile.is_open())
        {

            for (int i = 0; i < data.size(); ++i)
            {

                copy(data[i].begin(), data[i].end(),
                    std::ostream_iterator<T>(outputFile, " "));

                outputFile << std::endl;
            }
        }
    }
}

template <class T>
void saveData(std::vector<T> data, const char* filename)
{

    if (onpe0)
    {

        std::ofstream outputFile(filename);

        outputFile.precision(16);

        if (outputFile.is_open())
        {

            copy(data.begin(), data.end(),
                std::ostream_iterator<T>(outputFile, " "));

            outputFile << std::endl;
        }
    }
}

#ifdef SCALAPACK
template <class T>
void saveData(dist_matrix::DistMatrix<T> data, const char* filename)
{
    dist_matrix::BlacsContext bc(MPI_COMM_WORLD, 1, 1);

    std::ofstream outputFile(filename);

    outputFile.precision(16);

    if (outputFile.is_open())
    {
        data.print(outputFile);
    }
}
#endif

#endif /* SAVEDATA_H_ */
