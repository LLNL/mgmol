// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_GRASSMANLINEMINIMIZATION_H
#define MGMOL_GRASSMANLINEMINIMIZATION_H

#include "Hamiltonian.h"
#include "Ions.h"
#include "OrbitalsStepper.h"
#include "ProjectedMatricesInterface.h"

#include <iostream>

template <class T>
class MGmol;

template <class T>
class GrassmanLineMinimization : public OrbitalsStepper<T>
{
protected:
    Hamiltonian<T>* hamiltonian_;
    ProjectedMatricesInterface* proj_matrices_;
    MGmol<T>* mgmol_strategy_;
    std::ostream& os_;

    static bool accelerate_;
    static bool pbset_;
    static bool conjugate_;

    Ions* ptr2ions_;

    // Pointers to residual and search direction
    T* new_grad_;
    T* new_pcgrad_;

    T* grad_;
    T* pcgrad_;
    T* sdir_;

    static Timer line_min_tm_;
    static Timer nl_update_tm_;
    static Timer comp_res_tm_;
    static Timer update_states_tm_;

    void update_states(T& orbitals, const double precond_factor);

public:
    GrassmanLineMinimization(Hamiltonian<T>* hamiltonian,
        ProjectedMatricesInterface* proj_matrices, MGmol<T>* mgmol_strategy,
        Ions& ions, std::ostream& os)
        : hamiltonian_(hamiltonian),
          proj_matrices_(proj_matrices),
          mgmol_strategy_(mgmol_strategy),
          os_(os)
    {
        ptr2ions_  = &ions;
        conjugate_ = false;
    }

    ~GrassmanLineMinimization() override
    {
        delete new_grad_;
        delete new_pcgrad_;
        delete grad_;
        delete pcgrad_;
        delete sdir_;
    }

    void setup(T&) override {};

    int updateWF(T& orbitals, Ions& ions, const double precond_factor,
        const bool orthof, T& work_orbitals, const bool accelerate,
        const bool print_res, const double atol) override;

    virtual void conjugate()                                               = 0;
    virtual double computeStepSize(T& orbitals)                            = 0;
    virtual void parallelTransportUpdate(const double lambda, T& orbitals) = 0;
    static void printTimers(std::ostream& os);
};
// Instantiate static variables here to avoid clang warnings
template <class T>
bool GrassmanLineMinimization<T>::pbset_ = false;
template <class T>
bool GrassmanLineMinimization<T>::accelerate_ = false;
template <class T>
bool GrassmanLineMinimization<T>::conjugate_ = false;
#endif
