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
    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();

    ROMPrivateOptions rom_options = ct.getROMOptions();
    /* type of variable we intend to run POD */
    ROMVariable rom_var = rom_options.variable;

    /* number of restart files, start/end indices */
    assert(rom_options.restart_file_minidx >= 0);
    assert(rom_options.restart_file_maxidx >= 0);
    const int minidx = rom_options.restart_file_minidx;
    const int maxidx = rom_options.restart_file_maxidx;
    const int num_restart = maxidx - minidx + 1;

    MGmol<OrbitalsType> *mgmol = static_cast<MGmol<OrbitalsType> *>(mgmol_);
    OrbitalsType *orbitals = mgmol->getOrbitals();
    Potentials& pot = mgmol->getHamiltonian()->potential();
    std::string filename;

    /* Determine basis prefix, dimension, and sample size */
    std::string basis_prefix = rom_options.basis_file;
    int dim;
    int totalSamples = num_restart;
    const int chrom_num = orbitals->chromatic_number();
    switch (rom_var)
    {
    case ROMVariable::ORBITALS:
        basis_prefix += "_orbitals";
        dim = orbitals->getLocNumpt();
        /* if orbitals, each sample have chromatic number of wave functions */
        totalSamples *= orbitals->chromatic_number();
        break;

    case ROMVariable::POTENTIAL:
        basis_prefix += "_potential";
        dim = pot.size();
        break;
    
    default:
        ct.global_exit(-1);
        break;
    }

    /* Initialize libROM classes */
    CAROM::Options svd_options(dim, totalSamples, 1);
    CAROM::BasisGenerator basis_generator(svd_options, false, basis_prefix);
    
    /* Collect the restart files */
    for (int k = minidx; k <= maxidx; k++)
    {
        filename = string_format(rom_options.restart_file_fmt, k);
        mgmol->loadRestartFile(filename);
        assert(dim == orbitals->getLocNumpt());
        assert(chrom_num == orbitals->chromatic_number());

        switch (rom_var)
        {
        case ROMVariable::ORBITALS:
            for (int i = 0; i < chrom_num; ++i)
                basis_generator.takeSample(orbitals->getPsi(i));
            break;

        case ROMVariable::POTENTIAL:
            basis_prefix += "_potential";
            /* we save total potential for now */
            basis_generator.takeSample(pot.vtot());
            break;
        }
    }
    basis_generator.writeSnapshot();
    basis_generator.endSamples();
}

template void readRestartFiles<LocGridOrbitals>(MGmolInterface *mgmol_);
template void readRestartFiles<ExtendedGridOrbitals>(MGmolInterface *mgmol_);