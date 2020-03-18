// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Ions.h"
#include "KBPsiMatrixSparse.h"
#include "MGmol.h"
#include "Mesh.h"
#include "Species.h"

#include <cassert>
#include <iostream>
#include <list>
using namespace std;

Timer get_kbpsi_tm("get_kbpsi");
Timer vnlpsi_tm("vnlpsi");

void get_vnlpsi(const Ions& ions, const vector<vector<int>>& subdomain_gids,
    const int color, const KBPsiMatrixSparse* const kbpsi, ORBDTYPE* const vpsi)
{
    vnlpsi_tm.start();

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    const int subdivx      = mymesh->subdivx();

    memset(vpsi, 0, mygrid.size() * sizeof(ORBDTYPE));

    const vector<Ion*>& local_ions(ions.overlappingNL_ions());

    for (int iloc = 0; iloc < subdivx; iloc++)
        if (subdomain_gids[iloc][color] != -1)
        {
            const int st = subdomain_gids[iloc][color];

            // Loop over all the ions if nl proj. overlaps with sub-domain
            vector<Ion*>::const_iterator ion = local_ions.begin();
            while (ion != local_ions.end())
            {
                vector<int> ion_gids;
                (*ion)->getGidsNLprojs(ion_gids);

                vector<short> signs;
                (*ion)->getKBsigns(signs);

                vector<double> kbcoeffs;
                (*ion)->getKBcoeffs(kbcoeffs);

                const std::shared_ptr<KBprojector> ion_kbproj((*ion)->kbproj());

                if (ion_kbproj->onlyOneProjector())
                {
                    assert(ion_gids.size() == 1);
                    const int gid = ion_gids[0];

                    const double coeff = kbpsi->getValIonState(gid, st)
                                         * kbcoeffs[0] * signs[0];
                    ion_kbproj->axpySKet(iloc, coeff, vpsi);
                }
                else
                {
                    vector<double> coeff;
                    const short nprojs = (short)ion_gids.size();
                    for (short i = 0; i < nprojs; i++)
                    {
                        const int gid = ion_gids[i];
                        coeff.push_back(kbpsi->getValIonState(gid, st)
                                        * kbcoeffs[i] * signs[i]);
                    }
                    ion_kbproj->axpyKet(iloc, coeff, vpsi);
                }

                ion++;

            } // end loop over ions

        } // end loop over iloc

    vnlpsi_tm.stop();
}
