// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

////////////////////////////////////////////////////////////////////////////////
//
// mcstep.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id$

int mcstep(double& stx, double& fx, double& dx, double& sty, double& fy,
    double& dy, double& stp, const double fp, const double dp, bool& brackt,
    const double stpmin, const double stpmax);
