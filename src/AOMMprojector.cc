// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "AOMMprojector.h"
#include "LocGridOrbitals.h"
#include "MasksSet.h"
#include "Mesh.h"
#include "ProjectedMatrices.h"
#include "ProjectedMatricesSparse.h"
#include "SubspaceProjector.h"

AOMMprojector::AOMMprojector(LocGridOrbitals& phi, LocalizationRegions& lrs)
{
    Control& ct     = *(Control::instance());
    Mesh* mymesh    = Mesh::instance();
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    bool with_spin      = (mmpi.nspin() > 1);
    const short subdivx = mymesh->subdivx();

    // radius of kernel functions
    const double radius = ct.AOMMradius();

    counter_ = 1;

    // threshold on distance between centers below which projectors are dropped
    const double threshold = ct.AOMMthresholdFactor() * radius;

    if (onpe0)
        cout << "AOMM: setup with radius " << radius << " and threshold "
             << threshold << endl;

    kernelMasks_ = new MasksSet(false);
    kernelMasks_->setup(lrs, radius);

    if (ct.short_sighted)
        kernel_proj_matrices_ = new ProjectedMatricesSparse(ct.numst, &lrs);
    else
        kernel_proj_matrices_ = new ProjectedMatrices(ct.numst, with_spin);

    // kernel functions use their own projected matrices and masks
    kernel_phi_
        = new LocGridOrbitals("AOMM", phi, kernel_proj_matrices_, kernelMasks_, 0);
    kernel_phi_->initGauss(0.5 * radius, lrs);

    kernel_phi_->setIterativeIndex(phi.getIterativeIndex() * 10);

    // truncate kernel functions
    kernel_phi_->applyMask();

    kernel_phi_->computeGramAndInvS(ct.verbose);

    kernelprojector_ = new SubspaceProjector(*kernel_phi_);

    matrix_mask_ = new SquareLocalMatrices<MATDTYPE>(
        subdivx, kernel_phi_->chromatic_number());
    lrs.getMatrixDistances(*matrix_mask_, phi.getOverlappingGids());

#if 0
    if( onpe0 )
        (*MPIdata::sout)<<"Matrix of distances"<<endl;
    matrix_mask_->print((*MPIdata::sout));
#endif

    // don't project out component corresponding to very close and very far
    // centers
    matrix_mask_->setMaskThreshold(threshold, ct.cut_radius - radius);
    //    matrix_mask_->setMaskThreshold(threshold, 10000.);
}

void AOMMprojector::resetProjectors(LocGridOrbitals& phi)
{
    if (onpe0) cout << "AOMM: reset projectors..." << endl;

    kernel_phi_->assign(phi);

    // truncate kernel functions
    kernel_phi_->applyMask();

    kernel_phi_->computeGramAndInvS(0);

    delete kernelprojector_;
    kernelprojector_ = new SubspaceProjector(*kernel_phi_);
}

void AOMMprojector::projectOut(LocGridOrbitals& phi)
{
    assert(kernelprojector_ != 0);
    assert(matrix_mask_ != 0);

    // if( counter_%10==0 )resetProjectors(phi);

    kernelprojector_->projectOut(phi, matrix_mask_);

    counter_++;
}

AOMMprojector::~AOMMprojector()
{
    assert(kernel_phi_ != 0);

    delete kernelprojector_;
    kernelprojector_ = 0;
    delete kernel_phi_;
    kernel_phi_ = 0;
    delete kernel_proj_matrices_;
    kernel_proj_matrices_ = 0;
    delete kernelMasks_;
    kernelMasks_ = 0;
    delete matrix_mask_;
    matrix_mask_ = 0;
}
