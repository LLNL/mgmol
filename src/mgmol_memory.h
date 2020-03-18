// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef INCLUDE_MGMOL_MEMORY_H_
#define INCLUDE_MGMOL_MEMORY_H_

#include "MPIdata.h"

#include <string>
#include <vector>

namespace MGmolMem
{

class Stamp
{
public:
    static std::string filename_;
    static int line_;

    Stamp(const std::string& filename, const int line)
    {
        line_     = line;
        filename_ = filename;
    }
};

template <class T>
inline T* operator*(const Stamp& ppair, T* p)
{
    return p;
}

} // namespace MGmolMem

#define MEMTRACK_NEW MGmolMem::Stamp(__FILE__, __LINE__) * new
#define new MEMTRACK_NEW

#endif
