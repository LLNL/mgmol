// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef OPTION_DESCRIPTION_H
#define OPTION_DESCRIPTION_H

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;

void setupGenericOption(po::options_description &generic, string &config_file, string &lrs_filename);

void setupConfigOption(po::options_description &config, string &constraints_filename);

void setupHiddenOption(po::options_description &hidden);

void setupROMConfigOption(po::options_description &rom_cfg);

#endif
