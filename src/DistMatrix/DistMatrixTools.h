// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "DistMatrix.h"

void getProcrustesTransform(dist_matrix::DistMatrix<DISTMATDTYPE>& q, 
                           dist_matrix::DistMatrix<DISTMATDTYPE>& yyt);
void getAlignFrobeniusTransform(dist_matrix::DistMatrix<DISTMATDTYPE>& omatrix, 
                                dist_matrix::DistMatrix<DISTMATDTYPE>& rotation);
