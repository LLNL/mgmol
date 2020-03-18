// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef IonData_H
#define IonData_H

#include <string>
#include <vector>

const unsigned short IonData_MaxStrLength = 10;
struct FixedLengthString
{
    char mystring[IonData_MaxStrLength];
};

class IonData
{
    friend class Ions;
    friend class Ion;

private:
    std::string ion_name;
    int atomic_num;
    unsigned int index;
    unsigned int nlproj_id;
    unsigned int rand_state[3];
    int atmove;
    double initial_position[3];
    double old_position[3];
    double current_position[3];
    double force[3];
    double velocity[3];
    double pmass;

public:
    void unpack(char*& cptr, int*& iptr, double*& dptr);

    static void packIonData(
        char* cbuff, int* ibuff, double* dbuff, std::vector<IonData>& data);
};

#endif
