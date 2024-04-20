// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "rom_workflows.h"
#include <memory>
#include <string>
#include <stdexcept>

template<typename ... Args>
std::string string_format( const std::string& format, Args ... args )
{
    int size_s = std::snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
    if( size_s <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
    auto size = static_cast<size_t>( size_s );
    std::unique_ptr<char[]> buf( new char[ size ] );
    std::snprintf( buf.get(), size, format.c_str(), args ... );
    return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}

template <class OrbitalsType>
void readRestartFiles(MGmolInterface *mgmol_)
{
    Control& ct              = *(Control::instance());
    ROMPrivateOptions rom_options = ct.getROMOptions();
    assert(rom_options.restart_file_minidx >= 0);
    assert(rom_options.restart_file_maxidx >= 0);
    const int minidx = rom_options.restart_file_minidx;
    const int maxidx = rom_options.restart_file_maxidx;
    const int num_restart = maxidx - minidx + 1;

    MGmol<OrbitalsType> *mgmol = static_cast<MGmol<OrbitalsType> *>(mgmol_);
    OrbitalsType *orbitals = nullptr;
    std::string filename;

    /* Read the first snapshot to determin dimension and number of snapshots */
    filename = string_format(rom_options.restart_file_fmt, minidx);
    orbitals = mgmol->loadOrbitalFromRestartFile(filename);
    const int dim = orbitals->getLocNumpt();
    const int chrom_num = orbitals->chromatic_number();
    const int totalSamples = orbitals->chromatic_number() * num_restart;
    delete orbitals;

    /* Initialize libROM classes */
    CAROM::Options svd_options(dim, totalSamples, 1);
    CAROM::BasisGenerator basis_generator(svd_options, false, rom_options.basis_file);
    
    /* Collect the restart files */
    for (int k = minidx; k <= maxidx; k++)
    {
        filename = string_format(rom_options.restart_file_fmt, k);
        orbitals = mgmol->loadOrbitalFromRestartFile(filename);
        assert(dim == orbitals->getLocNumpt());
        assert(chrom_num == orbitals->chromatic_number());

        for (int i = 0; i < chrom_num; ++i)
            basis_generator.takeSample(orbitals->getPsi(i));

        delete orbitals;
    }
    basis_generator.writeSnapshot();
    basis_generator.endSamples();
}

template void readRestartFiles<LocGridOrbitals>(MGmolInterface *mgmol_);
template void readRestartFiles<ExtendedGridOrbitals>(MGmolInterface *mgmol_);