// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_AOMMPROJECTOR_H
#define MGMOL_AOMMPROJECTOR_H

#include "LocGridOrbitals.h"
#include "SquareLocalMatrices.h"
#include "SubspaceProjector.h"
#include "global.h"

class MasksSet;
class ProjectedMatricesInterface;
class LocalizationRegions;

// implements AOMM algorithm ("kernel" function projectors)
// by E. Tsuchida, J. Phys. Soc. Japan 76 (3), 2007, p. 034708
class AOMMprojector
{
private:
    LocGridOrbitals* kernel_phi_;

    SubspaceProjector<LocGridOrbitals>* kernelprojector_;

    MasksSet* kernelMasks_;

    ProjectedMatricesInterface* kernel_proj_matrices_;

    SquareLocalMatrices<MATDTYPE>* matrix_mask_;

    short counter_;

public:
    AOMMprojector(LocGridOrbitals& phi, LocalizationRegions& lrs);
    ~AOMMprojector();

    void projectOut(LocGridOrbitals& phi);

    void resetProjectors(LocGridOrbitals& phi);
};

#endif
