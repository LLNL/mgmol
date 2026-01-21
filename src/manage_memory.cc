// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "BlockVector.h"
#include "Control.h"
#include "global.h"

// Increase memory slots in BlockVector as needed based on runtime
// options
template <typename DataType, typename MemorySpaceType>
void increaseMemorySlotsForOrbitals()
{
    Control& ct = *(Control::instance());

    switch (ct.OuterSolver())
    {
        case OuterSolverType::ABPG:
        {
            // r_k-1, phi_k-1
            BlockVector<DataType, MemorySpaceType>::incMaxAllocInstances(2);
            break;
        }
        case OuterSolverType::PolakRibiere:
        {
            // r_k-1, z_k, z_k-1, p_k
            BlockVector<DataType, MemorySpaceType>::incMaxAllocInstances(4);
            break;
        }
        case OuterSolverType::Davidson:
        {
            BlockVector<DataType, MemorySpaceType>::incMaxAllocInstances(2);
            break;
        }
        default:
            break;
    }

    switch (ct.WFExtrapolation())
    {
        case WFExtrapolationType::Reversible:
        {
            BlockVector<DataType, MemorySpaceType>::incMaxAllocInstances(2);
            break;
        }
        case WFExtrapolationType::Order2:
        {
            BlockVector<DataType, MemorySpaceType>::incMaxAllocInstances(1);
            break;
        }
        case WFExtrapolationType::Order3:
        {
            BlockVector<DataType, MemorySpaceType>::incMaxAllocInstances(2);
            break;
        }
        default:
            break;
    }

    for (short i = 1; i < ct.wf_m; i++)
        BlockVector<DataType, MemorySpaceType>::incMaxAllocInstances(2);
    if (ct.use_kernel_functions)
        BlockVector<DataType, MemorySpaceType>::incMaxAllocInstances(1);

    switch (ct.AtomsDynamic())
    {
        case AtomsDynamicType::LBFGS:
            BlockVector<DataType, MemorySpaceType>::incMaxAllocInstances(1);
            break;
        default:
            break;
    }
}

template void increaseMemorySlotsForOrbitals<ORBDTYPE, MemorySpace::Host>();
#ifdef HAVE_MAGMA
template void increaseMemorySlotsForOrbitals<ORBDTYPE, MemorySpace::Device>();
#endif
