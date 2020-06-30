// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "DistMatrix.h"

void getProcrustesTransform(dist_matrix::DistMatrix<DISTMATDTYPE>& q,
    dist_matrix::DistMatrix<DISTMATDTYPE>& yyt);
void getAlignFrobeniusTransform(dist_matrix::DistMatrix<DISTMATDTYPE>& omatrix,
    dist_matrix::DistMatrix<DISTMATDTYPE>& rotation);
void rotateSym(dist_matrix::DistMatrix<DISTMATDTYPE>& mat,
    const dist_matrix::DistMatrix<DISTMATDTYPE>& rotation_matrix,
    dist_matrix::DistMatrix<DISTMATDTYPE>& work);
void sqrtDistMatrix(dist_matrix::DistMatrix<DISTMATDTYPE>& u);
