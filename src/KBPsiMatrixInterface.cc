// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "KBPsiMatrixInterface.h"
#include "Ion.h"
#include "KBprojectorSparse.h"
#include "Species.h"

using namespace std;

Timer KBPsiMatrixInterface::computeLocalElement_tm_(
    "KBPsiMatrixInterface::computeLocalElement");

void KBPsiMatrixInterface::computeLocalElement(Ion& ion, const int istate,
    const int iloc, const ORBDTYPE* const psi, const bool flag)
{
    assert(istate != -1);
    assert(iloc >= 0);

    vector<int> gids;
    ion.getGidsNLprojs(gids);
    if (gids.empty()) return;

    std::shared_ptr<KBprojector> ion_kbproj(ion.kbproj());

    if (!ion_kbproj->overlaps(iloc)) return;

    if (omp_get_thread_num() == 0) computeLocalElement_tm_.start();

    ion_kbproj->registerPsi(iloc, psi);

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    const double vel       = mygrid.vel();

    // Loop over projectors for this ion and evaluate
    // <KB|psi>
    const short nprojs = (short)gids.size();
    for (short i = 0; i < nprojs; i++)
    {
        const int gid    = gids[i];
        const double val = vel * ion_kbproj->dotPsi(iloc, i);
        if (flag)
        {
            addKBBPsi(gid, istate, val);
        }
        else
        {
            addKBPsi(gid, istate, val);
        }
    }

    if (omp_get_thread_num() == 0) computeLocalElement_tm_.stop();
}

void KBPsiMatrixInterface::printTimers(ostream& os)
{
    computeLocalElement_tm_.print(os);
}
