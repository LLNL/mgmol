#include "Control.h"
#include "BlockVector.h"
#include "global.h"

// Increase memory slots in BlockVector as needed based on runtime
// options
void increaseMemorySlotsForOrbitals()
{
    Control& ct     = *(Control::instance());

    if (ct.WFExtrapolation() == WFExtrapolationType::Reversible)
    {
        BlockVector<ORBDTYPE>::incMaxAllocInstances(2);
    }
    if (ct.WFExtrapolation() == WFExtrapolationType::Order3)
        BlockVector<ORBDTYPE>::incMaxAllocInstances(1);
    for (short i = 1; i < ct.wf_m; i++)
        BlockVector<ORBDTYPE>::incMaxAllocInstances(2);
    if (ct.use_kernel_functions) BlockVector<ORBDTYPE>::incMaxAllocInstances(1);
}

