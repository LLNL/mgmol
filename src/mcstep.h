// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
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
