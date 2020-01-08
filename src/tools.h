// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef TOOLS_H
#define TOOLS_H

#include "DistMatrix.h"

#include <fstream>
#include <string>

class Vector3D;

void noMoreMemory();
void stripLeadingAndTrailingBlanks(std::string& stringToModify);
void read_comments(std::ifstream& tfile);
void rotateSym(dist_matrix::DistMatrix<DISTMATDTYPE>& mat,
    const dist_matrix::DistMatrix<DISTMATDTYPE>& rotation_matrix,
    dist_matrix::DistMatrix<DISTMATDTYPE>& work);
void setSparseDistMatriConsolidationNumber(const int npes);
void reduceBytes(std::vector<char>& val, const MPI_Comm comm);
void arrayops(const double* const a, const double* const b, const double s,
    const double e, const int dim, double* result);
void printWithTimeStamp(const std::string& string2print, std::ostream& os);
bool isOverlaping(const Vector3D& center, const float radius);
void exitWithErrorMessage(const std::string& name);
bool fileExists(const char* file);
void stripName(std::string& stringToModify);
void appendNumber(std::string& name, const int number);
double ran0();
bool checkValidName(const std::string& name);
double minQuadPolynomial(const double e0, const double e1, const double de0,
    const bool print_flag, std::ostream& os);
double minQuadPolynomialFrom3values(const double e0, const double e1,
    const double e12, const bool print_flag, std::ostream& os);
void getkvector(const int index, const int kmax, int kvector[3]);

#endif
