// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_DIELECTRICCONTROL
#define MGMOL_DIELECTRICCONTROL

#include "Electrostatic.h"
#include "Potentials.h"

class DielectricControl
{
private:
    Potentials& pot_;
    Electrostatic* electrostat_;

    std::ostream& os_;
    const bool onpe0_;

    bool pbset_;

public:
    DielectricControl(Potentials& pot, Electrostatic* electrostat,
        std::ostream& os, const bool onpe0);

    void activate(const int it_scf, const double deig2);
};

#endif
