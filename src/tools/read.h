// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <fstream>
#include <iostream>

template <class T>
void read_data(T& data, std::ifstream& tfile)
{
    tfile >> data;
    while (tfile.get() != '\n')
        ; // read end of line
    char cc = (char)tfile.peek();
    while ((cc) == ('#') || (cc == '\n') || (cc == ' '))
    {
        while (tfile.get() != '\n')
            ;
        cc = (char)tfile.peek();
    }
}

static void read_comments(std::ifstream& tfile)
{
    while (tfile.get() != '\n')
        ;
    char cc = (char)tfile.peek();
    while ((cc) == ('#') || (cc == '\n'))
    {
        while (tfile.get() != '\n')
            ;
        cc = (char)tfile.peek();
    }
}
