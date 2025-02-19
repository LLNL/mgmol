// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "IonData.h"
#include "tools.h"
#include <string.h>

void IonData::unpack(char*& cptr, int*& iptr, double*& dptr)
{
    // get name
    std::string name(cptr, IonData_MaxStrLength);
    stripLeadingAndTrailingBlanks(name);
    ion_name = name;
    cptr += IonData_MaxStrLength;
    // get atomic_num
    atomic_num = *(iptr++);
    // get index
    index = *(iptr++);
    // get nlproj_id
    nlproj_id = *(iptr++);
    // get atmove
    atmove = *(iptr++);
    // get rand_states
    for (int st = 0; st < 3; st++)
        rand_state[st] = *(iptr++);

    // get pmass
    pmass = *(dptr++);
    // get initial_position
    for (short st = 0; st < 3; st++)
        initial_position[st] = *(dptr++);
    // get old_position
    for (short st = 0; st < 3; st++)
        old_position[st] = *(dptr++);
    // get current_position
    for (short st = 0; st < 3; st++)
        current_position[st] = *(dptr++);
    // get velocity
    for (short st = 0; st < 3; st++)
        velocity[st] = *(dptr++);
    // get force
    for (short st = 0; st < 3; st++)
        force[st] = *(dptr++);
}

// pack Ions data for communication
void IonData::packIonData(
    char* cbuff, int* ibuff, double* dbuff, std::vector<IonData>& data)
{
    // pack ion_names buffer
    int idx = 0;
    for (std::vector<IonData>::iterator idata = data.begin();
         idata != data.end(); ++idata)
    {
        std::string s = (*idata).ion_name;
        FixedLengthString t;
        strncpy(t.mystring, s.c_str(), IonData_MaxStrLength - 1);
        memcpy(&cbuff[idx], t.mystring, IonData_MaxStrLength);
        idx += IonData_MaxStrLength;
    }

    // pack integer datatypes
    int* iptr = &ibuff[0];
    // first pack local data size
    *(iptr++)                            = data.size();
    std::vector<IonData>::iterator idata = data.begin();
    while (idata != data.end())
    {
        // pack atomic_num
        *(iptr++) = (*idata).atomic_num;
        // pack ion_index
        *(iptr++) = (*idata).index;
        // pack ion nlproj_id
        *(iptr++) = (*idata).nlproj_id;
        // pack atmove
        *(iptr++) = (*idata).atmove;
        // pack rand_states
        for (short i = 0; i < 3; i++)
            *(iptr++) = (*idata).rand_state[i];

        idata++;
    }

    // pack double datatypes
    double* dptr = &dbuff[0];
    idata        = data.begin();
    while (idata != data.end())
    {
        // pack pmass
        *(dptr++) = (*idata).pmass;
        // pack initial_position
        for (short i = 0; i < 3; i++)
            *(dptr++) = (*idata).initial_position[i];
        // pack old_position
        for (short i = 0; i < 3; i++)
            *(dptr++) = (*idata).old_position[i];
        // pack current_position
        for (short i = 0; i < 3; i++)
            *(dptr++) = (*idata).current_position[i];
        // pack velocities
        for (short i = 0; i < 3; i++)
            *(dptr++) = (*idata).velocity[i];
        // pack forces
        for (short i = 0; i < 3; i++)
            *(dptr++) = (*idata).force[i];

        idata++;
    }
}
