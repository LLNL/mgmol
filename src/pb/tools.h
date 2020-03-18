// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef PB_TOOLS_H
#define PB_TOOLS_H

#include <stdio.h>
#include <string>

namespace pb
{

FILE* open_file(const std::string&, const std::string&);
double timer(void);
}

#endif
