// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
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

// implements AOMM algorithm ("kernel" function projectors)
// by E. Tsuchida, J. Phys. Soc. Japan 76 (3), 2007, p. 034708
class AOMMprojector
{
private:
    LocGridOrbitals<ORBDTYPE>* kernel_phi_;

    SubspaceProjector<LocGridOrbitals<ORBDTYPE>>* kernelprojector_;

    MasksSet* kernelMasks_;

    ProjectedMatricesInterface* kernel_proj_matrices_;

    SquareLocalMatrices<MATDTYPE, MemorySpace::Host>* matrix_mask_;

    short counter_;

public:
    AOMMprojector(LocGridOrbitals<ORBDTYPE>& phi,
        const std::shared_ptr<LocalizationRegions>& lrs);
    ~AOMMprojector();

    void projectOut(LocGridOrbitals<ORBDTYPE>& phi);

    void resetProjectors(LocGridOrbitals<ORBDTYPE>& phi);
};

#endif
