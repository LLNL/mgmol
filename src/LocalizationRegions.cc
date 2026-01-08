// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "global.h"

#include "LocGridOrbitals.h"
#include "LocalizationRegions.h"
#include "Mesh.h"
#include "SpreadsAndCenters.h"
#include "SquareLocalMatrices.h"
#include "SymmetricPair.h"
#include "hdf_tools.h"
#include "tools.h"

#ifdef MGMOL_USE_SCALAPACK
#include "OrbitalsTransform.h"
#endif

#include <algorithm>
#include <iomanip>
#include <sstream>

using namespace std;

double LocalizationRegions::fourthirdpi_ = 4. * M_PI / 3.;
Timer LocalizationRegions::setupLocalRegionsFromOverlapRegions_tm_(
    "LocalizationRegions::setupLocRegFromOvlpRegs");
Timer LocalizationRegions::syncCenter_tm_("LocalizationRegions::syncCenters");

Timer LR_moveTo_tm("LR_moveTo");

void LocalizationRegions::clear()
{
    // if( onpe0 )
    //    (*MPIdata::sout)<<"LocalizationRegions::clear()"<<endl;

    clearOldCenters();
    all_regions_.clear();
    local_regions_.clear();
    volume_ = 0.;
}

void LocalizationRegions::clearOldCenters()
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 3)
        (*MPIdata::sout) << "LocalizationRegions::clearOldCenters()" << endl;

    for (unsigned int i = 0; i < old_centers_.size(); i++)
    {
        old_centers_[i].clear();
    }

    old_centers_set_    = false;
    num_extrapolations_ = 0;
}

float LocalizationRegions::max_radii() const
{
    float maxRadii                    = 0.;
    vector<LRData>::const_iterator pr = local_regions_.begin();
    while (pr != local_regions_.end())
    {
        maxRadii = max(maxRadii, pr->radius);
        pr++;
    }

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allreduce(&maxRadii, 1, MPI_MAX);

    return maxRadii;
}

float LocalizationRegions::min_radii() const
{
    float minRadii                    = 100000.;
    vector<LRData>::const_iterator ir = local_regions_.begin();
    while (ir != local_regions_.end())
    {
        minRadii = min(minRadii, ir->radius);
        ir++;
    }

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allreduce(&minRadii, 1, MPI_MIN);

    return minRadii;
}

void LocalizationRegions::extrapolateCentersLinear(
    const bool force, const bool move_centers)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 1)
        (*MPIdata::sout) << "LocalizationRegions::extrapolateCentersLinear()..."
                         << endl;

    if (!old_centers_set_) resetOldCenters();

    assert((int)overlap_regions_.size() == (int)old_centers_[1].size());

    vector<vector<Vector3D>::iterator> it_old;

    const int old_centers_size = old_centers_.size();
    it_old.reserve(old_centers_size);
    for (int i = 0; i < old_centers_size; i++)
    {
        it_old.push_back(old_centers_[i].begin());
    }

    for (vector<LRData>::iterator it = overlap_regions_.begin();
         it != overlap_regions_.end(); ++it)
    {
        assert(it_old[1] != old_centers_[1].end());
        assert(it_old[2] != old_centers_[2].end());

        Vector3D d1 = it->center.vminimage(*it_old[1], cell_, ct.bcPoisson);
        Vector3D d2 = it->center.vminimage(*it_old[2], cell_, ct.bcPoisson);

        // compare last move with previous one
        d2 *= 0.5;
        d2 -= d1;
        double dev = length(d2);

        for (unsigned int i = 0; i < old_centers_.size() - 2; i++)
        {
            *it_old[old_centers_.size() - 1 - i]
                = *it_old[old_centers_.size() - 2 - i];
        }
        *it_old[1] = it->center;

        // if the last two moves were similar, do interpolation
        if ((dev < tol_displ_ || force) && ct.dt > 0)
        {
            *it_old[0] = it->center + d1;

            if (move_centers)
            {
                (it->center) = *it_old[0];
            }
        }
        else
        {
            *it_old[0] = it->center;
        }

        for (unsigned int i = 0; i < old_centers_.size(); i++)
        {
            ++it_old[i];
        }
    }

    if (num_extrapolations_ < static_cast<int>(old_centers_.size()))
        num_extrapolations_++;

    if (move_centers)
    {
        moveToExtrapolatedCenters();
    }
}

void LocalizationRegions::moveToExtrapolatedCenters()
{
    vector<Vector3D>::const_iterator it_extrapolated = old_centers_[0].begin();

    for (vector<LRData>::iterator it = overlap_regions_.begin();
         it != overlap_regions_.end(); ++it)
    {
        it->center = *it_extrapolated;
        ++it_extrapolated;
    }

    updateOverlapRegions();
    setupLocalRegionsFromOverlapRegions();
}

// "reversible" Verlet scheme
void LocalizationRegions::extrapolateCentersVerlet(
    const bool force, const bool move_centers)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 1)
        (*MPIdata::sout) << "LocalizationRegions::extrapolateCentersVerlet()..."
                         << endl;

    if (!old_centers_set_) resetOldCenters();

    assert((int)overlap_regions_.size() == (int)old_centers_[1].size());

    vector<vector<Vector3D>::iterator> it_old;

    const int old_centers_size = old_centers_.size();
    it_old.reserve(old_centers_size);
    for (int i = 0; i < old_centers_size; i++)
    {
        it_old.push_back(old_centers_[i].begin());
    }

    for (vector<LRData>::iterator it = overlap_regions_.begin();
         it != overlap_regions_.end(); ++it)
    {
        assert(it_old[1] != old_centers_[1].end());
        assert(it_old[2] != old_centers_[2].end());

        Vector3D d = it_old[1]->vminimage(*it_old[2], cell_, ct.bcPoisson);

        // extrapolation
        if (force)
        {
            *it_old[0] = it->center + d;

            if (move_centers)
            {
                it->center = *it_old[0];
            }
        }
        else
        {
            Vector3D d2 = it->center.vminimage(*it_old[1], cell_, ct.bcPoisson);
            double dev  = length(d2);
            if (dev < tol_displ_)
            {
                *it_old[0] = it->center + d;

                if (move_centers)
                {
                    it->center = *it_old[0];
                }
            }
            else
            {
                *it_old[0] = it->center;
            }
        }

        *it_old[2] = *it_old[1];
        // save extrapolated center into old1 (to be used for next
        // extrapolation)
        *it_old[1] = it->center;

        for (unsigned int i = 0; i < old_centers_.size(); i++)
        {
            ++it_old[i];
        }
    }

    updateOverlapRegions();

    setupLocalRegionsFromOverlapRegions();
}

void LocalizationRegions::extrapolateCentersQuadratic()
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 1)
        (*MPIdata::sout)
            << "LocalizationRegions::extrapolateCentersQuadratic()..." << endl;

    if (old_centers_[2].empty())
    {
        extrapolateCentersLinear(true, true);
        return;
    }

    vector<Vector3D>::iterator it_old1 = old_centers_[1].begin();
    vector<Vector3D>::iterator it_old2 = old_centers_[2].begin();

    for (vector<LRData>::iterator it = overlap_regions_.begin();
         it != overlap_regions_.end(); ++it)
    {
        assert(it_old1 != old_centers_[1].end());
        assert(it_old2 != old_centers_[2].end());

        Vector3D d1 = it->center.vminimage(*it_old1, cell_, ct.bcPoisson);
        Vector3D y2(*it_old2);

        *it_old2 = *it_old1;
        *it_old1 = it->center;

        (it->center) = y2;
        (it->center) += 3. * d1;

        ++it_old1;
        ++it_old2;
    }

    updateOverlapRegions();

    setupLocalRegionsFromOverlapRegions();
}

void LocalizationRegions::resetOldCenters()
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 1)
        (*MPIdata::sout) << "LocalizationRegions::resetOldCenters()..." << endl;

    for (unsigned int i = 0; i < old_centers_.size(); i++)
    {
        old_centers_[i].clear();
    }

    for (vector<LRData>::const_iterator it = overlap_regions_.begin();
         it != overlap_regions_.end(); ++it)
    {
        for (unsigned int i = 0; i < old_centers_.size(); i++)
        {
            old_centers_[i].push_back(it->center);
        }
    }

    old_centers_set_    = true;
    num_extrapolations_ = 0;
}

template <class T>
float LocalizationRegions::move(const SpreadsAndCenters<T>& sc, const bool flag)
{
    vector<Vector3D> centers;
    sc.computeCenters(centers);

    max_displ_ = moveTo(centers);

    if (flag) resetOldCenters();

    return max_displ_;
}

float LocalizationRegions::getMeanRadius()
{
    float rcube = 0.;
    int nlr     = 0;
    for (std::vector<LRData>::iterator it = local_regions_.begin();
         it != local_regions_.end(); it++)
    {
        const float alpha = it->radius;
        rcube += (alpha * alpha * alpha);
        nlr++;
    }

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allreduce(&rcube, 1, MPI_SUM);
    mmpi.allreduce(&nlr, 1, MPI_SUM);

    assert(rcube > 0.);

    return cbrt(rcube / (float)nlr);
}

template <class T>
float LocalizationRegions::updateRadii(
    const SpreadsAndCenters<T>& sc, const float ratio)
{
    for (vector<LRData>::iterator it = overlap_regions_.begin();
         it != overlap_regions_.end(); ++it)
    {
        const double spread = sc.computeSpread(it->gid[0]);
        it->radius          = ratio * spread;
    }

    return getMeanRadius();
}

template <class T>
float LocalizationRegions::updateRadiiConstVol(const SpreadsAndCenters<T>& sc)
{
    assert(volume_ > 0.);

    // first assign new local radii to compute new volume
    for (vector<LRData>::iterator it = local_regions_.begin();
         it != local_regions_.end(); ++it)
    {
        const float spread = sc.computeSpread(it->gid[0]);
        assert(spread > 0.);
        it->radius = spread;
    }
    // compute resulting new volume
    float new_vol = computeVolumePrivate();
    // get ratio = old volume/ new volume
    const float ratio = cbrt(volume_ / new_vol);
    if (onpe0) (*MPIdata::sout) << "Ratio Rc/spread: " << ratio << endl;

    assert(ratio > 0.);

    // reset new radii to conserve volume
    return updateRadii(sc, ratio);
}

#ifdef MGMOL_USE_SCALAPACK
float LocalizationRegions::updateRadii(
    const OrbitalsTransform* ot, const float ratio)
{
    for (vector<LRData>::iterator it = overlap_regions_.begin();
         it != overlap_regions_.end(); ++it)
    {
        const double spread = ot->spread(it->gid[0]);
        it->radius          = ratio * spread;
    }

    return getMeanRadius();
}
#endif

float LocalizationRegions::moveTo(const vector<Vector3D>& target_centers)
{
    assert(overlap_regions_.size() == target_centers.size());

    float max_displ = 0.;
    assert(nglobal_ > 0);
    Control& ct = *(Control::instance());

    // existing centers...
    LR_moveTo_tm.start();

    int i                       = 0;
    vector<LRData>::iterator it = overlap_regions_.begin();
    while (it != overlap_regions_.end())
    {
        Vector3D& region_center(it->center);
        Vector3D delta
            = target_centers[i].vminimage(region_center, cell_, ct.bcPoisson);

        float nrm2 = length(delta);
        if (nrm2 > max_displ) max_displ = nrm2;
        if (nrm2 > 1.) delta /= nrm2;
        region_center.axpy(1., delta);

        it++;
        i++;
    }

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allreduce(&max_displ, 1, MPI_MAX);

    if (!old_centers_set_) resetOldCenters();

    updateOverlapRegions();
    LR_moveTo_tm.stop();

    setupLocalRegionsFromOverlapRegions();

    if (onpe0)
        (*MPIdata::sout) << setprecision(4) << fixed
                         << "Max. Orbital center move: " << max_displ << endl;

    iterative_index_++;

    return max_displ;
}

void LocalizationRegions::setUniformRadii(const float radius)
{
    const int n = (int)overlap_regions_.size();
    for (int i = 0; i < n; i++)
    {
        overlap_regions_[i].radius = radius;
    }
}

float LocalizationRegions::computeVolume()
{
    volume_ = 0.;
    for (vector<LRData>::const_iterator it = local_regions_.begin();
         it != local_regions_.end(); ++it)
    {
        const float alpha = it->radius;
        volume_ += (alpha * alpha * alpha);
    }
    volume_ *= fourthirdpi_;

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allreduce(&volume_, 1, MPI_SUM);

    return volume_;
}

void LocalizationRegions::bcastLRs()
{
    Control& ct = *(Control::instance());

    if (ct.verbose > 0)
        printWithTimeStamp("LocalizationRegions::bcast()...", cout);

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    MPI_Comm comm = mmpi.commSameSpin();

    nglobal_     = (int)all_regions_.size();
    int retbcast = MPI_Bcast(&nglobal_, 1, MPI_INT, 0, comm);
    if (retbcast != MPI_SUCCESS)
    {
        cerr << "ERROR!!!! LocalizationRegions::bcast(), Failure in MPI_Bcast "
                "of 'nglobal_'!!!"
             << endl;
        MPI_Abort(comm, EXIT_FAILURE);
    }

    if (nglobal_ == 0) return;
    if (!onpe0) all_regions_.resize(nglobal_);

    // Bcast centers
    {
        vector<double> buffer_centers(3 * nglobal_);
        if (mmpi.instancePE0())
        {
            int i                       = 0;
            vector<LRData>::iterator it = all_regions_.begin();
            while (it != all_regions_.end())
            {
                buffer_centers[i++] = it->center[0];
                buffer_centers[i++] = it->center[1];
                buffer_centers[i++] = it->center[2];
                it++;
            }
            assert(i == 3 * nglobal_);
        }

        retbcast
            = MPI_Bcast(&buffer_centers[0], 3 * nglobal_, MPI_DOUBLE, 0, comm);
        if (retbcast != MPI_SUCCESS)
        {
            cerr << "ERROR!!!! LocalizationRegions::bcast(), Failure in "
                    "MPI_Bcast of 'centers'!!!"
                 << endl;
            MPI_Abort(comm, EXIT_FAILURE);
        }
        vector<LRData>::iterator it = all_regions_.begin();
        int i                       = 0;
        while (it != all_regions_.end())
        {
            it->center[0] = buffer_centers[i++];
            it->center[1] = buffer_centers[i++];
            it->center[2] = buffer_centers[i++];
            it++;
        }
        assert(i == 3 * nglobal_);
    }

    // Bcast radius
    {
        vector<float> buffer_radius(nglobal_);
        int i                       = 0;
        vector<LRData>::iterator it = all_regions_.begin();
        while (it != all_regions_.end())
        {
            buffer_radius[i++] = it->radius;
            it++;
        }
        assert(i == nglobal_);

        retbcast = MPI_Bcast(&buffer_radius[0], nglobal_, MPI_FLOAT, 0, comm);
        if (retbcast != MPI_SUCCESS)
        {
            cerr << "ERROR!!!! LocalizationRegions::bcast(), Failure in "
                    "MPI_Bcast of 'radius'!!!"
                 << endl;
            MPI_Abort(comm, EXIT_FAILURE);
        }
        it = all_regions_.begin();
        i  = 0;
        while (it != all_regions_.end())
        {
            it->radius = buffer_radius[i++];
            it++;
        }
        assert(i == nglobal_);
    }

    // bcast gid
    {
        vector<short> bufsizes;
        bufsizes.reserve(nglobal_);
        vector<LRData>::iterator it = all_regions_.begin();
        while (it != all_regions_.end())
        {
            bufsizes.push_back((short)it->gid.size());
            it++;
        }
        MPI_Bcast(&bufsizes[0], nglobal_, MPI_SHORT, 0, comm);
        it    = all_regions_.begin();
        int i = 0;
        if (!mmpi.instancePE0())
        {
            while (it != all_regions_.end())
            {
                assert(bufsizes[i] > 0);
                it->gid.resize(bufsizes[i++]);
                it++;
            }
            assert(i == nglobal_);
        }

        int sumsizes = 0;
        for (int j = 0; j < nglobal_; j++)
        {
            sumsizes += bufsizes[j];
        }
        assert(sumsizes == nglobal_);
        vector<int> bufgids(sumsizes);

        // pack
        if (mmpi.instancePE0())
        {
            it = all_regions_.begin();
            i  = 0;
            while (it != all_regions_.end())
            {
                vector<int>::iterator git = it->gid.begin();
                while (git != it->gid.end())
                {
                    bufgids[i++] = *git;
                    git++;
                }
                it++;
            }
            assert(i == sumsizes);
        }

        MPI_Bcast(&bufgids[0], sumsizes, MPI_INT, 0, comm);

        // unpack
        assert(nglobal_ > 0);
        it = all_regions_.begin();
        i  = 0;
        while (it != all_regions_.end())
        {
            vector<int>::iterator git = it->gid.begin();
            assert(it->gid.size() > 0);
            while (git != it->gid.end())
            {
                *git = bufgids[i++];
                git++;
            }
            it++;
        }
        assert(i == sumsizes);
    }

    retbcast = MPI_Bcast(&volume_, 1, MPI_FLOAT, 0, comm);
    if (retbcast != MPI_SUCCESS)
    {
        cerr << "ERROR!!!! LocalizationRegions::bcast(), Failure in MPI_Bcast "
                "of 'volume'!!!"
             << endl;
        MPI_Abort(comm, EXIT_FAILURE);
    }
    if (ct.verbose > 0)
        printWithTimeStamp("LocalizationRegions::bcast() done...", cout);
}

void LocalizationRegions::printAllRegions(ostream& os)
{
    // sync centers to update regions_ data
    syncCenters();

    const float minr = min_radii();
    const float maxr = max_radii();

    if (onpe0)
    {
        os << "Localization Regions:" << endl;
        os << "Min. Radius: " << minr << endl;
        os << "Max. Radius: " << maxr << endl;
        for (int i = 0; i < nglobal_; i++)
        {
            os << "@@ " << setw(4) << i + 1 << " Center ";
            os.setf(ios::fixed, ios::floatfield);
            os.setf(ios::right, ios::adjustfield);
            os << all_regions_[i].center;
            os << ", radius " << all_regions_[i].radius << ", multiplicity "
               << all_regions_[i].gid.size() << endl;
        }
    }
}

/* This function needs updating to use local info only */
float LocalizationRegions::imposeMinDistanceBetweenCentersGlobal(
    const float drmin)
{

    /* First call syncCenters() to update regions_ */
    //    assert(!local_regions_.empty());
    syncCenters();

    float maxdispl = 0.;
    Control& ct    = *(Control::instance());
    int step       = 1;
    int npairs     = 1000;
    while (npairs > 0)
    {
        npairs = 0;
        vector<float> forces(3 * nglobal_, 0.);

        for (int index0 = 0; index0 < nglobal_; index0++)
        {
            float* pforce0 = &forces[3 * index0];
            for (int index1 = 0; index1 < index0; index1++)
            {
                float* pforce1 = &forces[3 * index1];
                Vector3D dr    = all_regions_[index1].center.vminimage(
                    all_regions_[index0].center, cell_, ct.bcPoisson);
                const float dr0 = length(dr);

                if (dr0 < drmin)
                {
                    const float factor = 0.5;
                    float fx, fy, fz;
                    if (dr0 > 1.e-6)
                    {
                        float af     = factor * (drmin - dr0 + 1.e-4);
                        float invdr0 = 1. / dr0;
                        fx           = af * dr[0] * invdr0;
                        fy           = af * dr[1] * invdr0;
                        fz           = af * dr[2] * invdr0;
                    }
                    else
                    {
                        fx = factor * drmin * ((float)rand() / (float)RAND_MAX);
                        fy = factor * drmin * ((float)rand() / (float)RAND_MAX);
                        fz = factor * drmin * ((float)rand() / (float)RAND_MAX);
                    }
                    pforce0[0] = pforce0[0] - fx;
                    pforce0[1] = pforce0[1] - fy;
                    pforce0[2] = pforce0[2] - fz;
                    pforce1[0] = pforce1[0] + fx;
                    pforce1[1] = pforce1[1] + fy;
                    pforce1[2] = pforce1[2] + fz;

                    npairs++;
                }
            }
        }
        step++;

        // move centers
        if (npairs > 0)
        {
            for (int index = 0; index < nglobal_; index++)
            {
                all_regions_[index].center[0]
                    = all_regions_[index].center[0] + forces[3 * index + 0];
                all_regions_[index].center[1]
                    = all_regions_[index].center[1] + forces[3 * index + 1];
                all_regions_[index].center[2]
                    = all_regions_[index].center[2] + forces[3 * index + 2];
                float norm2f = forces[3 * index + 0] * forces[3 * index + 0];
                norm2f = norm2f + forces[3 * index + 1] * forces[3 * index + 1];
                norm2f = norm2f + forces[3 * index + 2] * forces[3 * index + 2];
                if (norm2f > maxdispl) maxdispl = norm2f;
            }
        }
    }
    /* now we need to update overlapping regions data members */
    if (!overlap_regions_.empty())
    {
        computeOverlapRegionsGlobal();
        computeOverlapGids();
        computeSubdivOverlapGids();
    }
    // update local regions
    setupLocalRegionsFromOverlapRegions();

    return sqrt(maxdispl);
}

float LocalizationRegions::getStatesWithClosestCenters(
    int* const st1, int* const st2, const set<SymmetricPair>& exclude_set) const
{
    float drmin2 = 1.e12;
    Control& ct  = *(Control::instance());

    for (int index0 = 0; index0 < (int)overlap_regions_.size(); index0++)
    {
        for (int index1 = 0; index1 < index0; index1++)
        {

            Vector3D dr = overlap_regions_[index1].center.vminimage(
                overlap_regions_[index0].center, cell_, ct.bcPoisson);
            const float dr2 = norm2(dr);
            // if(onpe0)(*MPIdata::sout)<<"dr2="<<dr2<<endl;
            if (dr2 < drmin2)
            {
                const int sti0 = overlap_regions_[index0].gid[0];
                const int sti1 = overlap_regions_[index1].gid[0];
                SymmetricPair pair(sti0, sti1);
                if (exclude_set.find(pair) == exclude_set.end())
                {
                    drmin2 = dr2;
                    *st1   = sti0;
                    *st2   = sti1;
                }
            }
        }
    }

    /* communicate info */
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    MPI_Comm comm = mmpi.commSameSpin();
    int npes;
    MPI_Comm_size(comm, &npes);

    int rank;
    mmpi.getRankMinVal(drmin2, &rank);
    int st[2] = { *st1, *st2 };
    mmpi.bcast(st, 2, rank);
    *st1 = st[0];
    *st2 = st[1];
    mmpi.bcast(&drmin2, 1, rank);

    if (onpe0) (*MPIdata::sout) << "drmin2=" << drmin2 << endl;
    return sqrt(drmin2);
}

// out: global indices of functions centered in local subdomain
// returned in vector<int>& indices
void LocalizationRegions::getLocalSubdomainIndices(vector<int>& indices)
{
    if (!gid_local_regions_set_)
    {
        setupLocalRegionsFromOverlapRegions();
    }

    indices.clear();
    for (vector<LRData>::const_iterator itg = local_regions_.begin();
         itg != local_regions_.end(); itg++)
    {
        for (vector<int>::const_iterator iitg = itg->gid.begin();
             iitg != itg->gid.end(); iitg++)
            indices.push_back(*iitg);
    }
}

bool LocalizationRegions::isLRcenterLocal(const LRData& lr) const
{
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    const Vector3D& center(lr.center);

    double t[3] = { center[0], center[1], center[2] };
    for (short i = 0; i < 3; i++)
    {
        while (t[i] >= subdom_end_[i])
            t[i] -= mygrid.ll(i);
        while (t[i] < subdom_origin_[i])
            t[i] += mygrid.ll(i);
        t[i] -= subdom_lower_left_[i];
    }
    if (t[0] >= 0. && t[0] < (subdom_lattice_[0]))
        if (t[1] >= 0. && t[1] < (subdom_lattice_[1]))
            if (t[2] >= 0. && t[2] < (subdom_lattice_[2]))
            {
                return true;
            }

    return false;
}

void LocalizationRegions::setupLocalRegionsFromOverlapRegions()
{

    setupLocalRegionsFromOverlapRegions_tm_.start();

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));

    Control& ct = *(Control::instance());

    if (ct.verbose > 0)
    {
        printWithTimeStamp(
            "LocalizationRegions::setupLocalRegionsFromOverlapRegions()...",
            (*MPIdata::sout));
        //(*MPIdata::sout)<<"Number of overlap regions: "
        //                <<overlap_regions_.size()<<endl;
    }

    local_regions_.clear();

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    int count              = 0;

    vector<LRData>::const_iterator itc = overlap_regions_.begin();
    while (itc != overlap_regions_.end())
    {
        const Vector3D& center(itc->center);
        double t[3] = { center[0], center[1], center[2] };
        for (short i = 0; i < 3; i++)
        {
            while (t[i] >= subdom_end_[i])
                t[i] -= mygrid.ll(i);
            while (t[i] < subdom_origin_[i])
                t[i] += mygrid.ll(i);
            t[i] -= subdom_lower_left_[i];
        }
        if (t[0] >= 0. && t[0] < (subdom_lattice_[0]))
            if (t[1] >= 0. && t[1] < (subdom_lattice_[1]))
                if (t[2] >= 0. && t[2] < (subdom_lattice_[2]))
                {
                    LRData data;
                    data.center = center;
                    data.radius = itc->radius;
                    data.gid    = itc->gid;
                    local_regions_.push_back(data);
                    count++;
                }

        itc++;
    }

    // cout<<"count="<<count<<", n. regions="<<nglobal_<<endl;
    mmpi.allreduce(&count, 1, MPI_SUM);
    if (count == 0)
    {
        cerr << "ERROR in distribution of localization centers: count=" << count
             << endl;
        cerr << "global number of regions=" << nglobal_ << endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    gid_local_regions_set_ = true;
    setupLocalRegionsFromOverlapRegions_tm_.stop();
}

void LocalizationRegions::setupSubdomainGeometry()
{

    Mesh* mymesh             = Mesh::instance();
    const pb::Grid& mygrid   = mymesh->grid();
    const pb::PEenv& myPEenv = mymesh->peenv();

    for (short i = 0; i < 3; i++)
    {
        subdom_origin_[i]  = mygrid.origin(i);
        subdom_lattice_[i] = mygrid.ll(i) / (double)(myPEenv.n_mpi_task(i));
    }
    for (short i = 0; i < 3; i++)
        subdom_end_[i] = subdom_origin_[i] + mygrid.ll(i);
    for (short i = 0; i < 3; i++)
        subdom_lower_left_[i] = (double)myPEenv.my_mpi(i) * subdom_lattice_[i]
                                + subdom_origin_[i];
    for (short i = 0; i < 3; i++)
        subdom_upper_right_[i] = subdom_lower_left_[i] + subdom_lattice_[i];
    for (short i = 0; i < 3; i++)
        subdom_center_[i]
            = 0.5 * (subdom_lower_left_[i] + subdom_upper_right_[i]);

    const double h[3] = { mygrid.hgrid(0), mygrid.hgrid(1), mygrid.hgrid(2) };
    subdom_delta_     = sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]) + 1.e-8;
    // Take into account Differential operator
    subdom_delta_ *= (double)mygrid.ghost_pt();

    // correct div_lattice to include only mesh points, not domain
    subdom_div_lattice_[0] = subdom_lattice_[0] - h[0];
    subdom_div_lattice_[1] = subdom_lattice_[1] - h[1];
    subdom_div_lattice_[2] = subdom_lattice_[2] - h[2];
}

bool LocalizationRegions::addOverlappingRegion(const LRData& region)
{
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    bool overlaps = false;

    const Vector3D& center(region.center);
    double t[3] = { center[0], center[1], center[2] };
    for (short i = 0; i < 3; i++)
    {
        while (t[i] >= subdom_center_[i] + 0.5 * mygrid.ll(i))
            t[i] -= mygrid.ll(i);
        while (t[i] < subdom_center_[i] - 0.5 * mygrid.ll(i))
            t[i] += mygrid.ll(i);

        t[i] -= subdom_lower_left_[i];
    }

    // find "p" closest point within subdomain
    double p[3] = { t[0], t[1], t[2] };
    if (p[0] < 0.) p[0] = 0.;
    if (p[1] < 0.) p[1] = 0.;
    if (p[2] < 0.) p[2] = 0.;
    if (p[0] > subdom_div_lattice_[0]) p[0] = subdom_div_lattice_[0];
    if (p[1] > subdom_div_lattice_[1]) p[1] = subdom_div_lattice_[1];
    if (p[2] > subdom_div_lattice_[2]) p[2] = subdom_div_lattice_[2];

    const float radius = region.radius + subdom_delta_;
    double d2 = (p[0] - t[0]) * (p[0] - t[0]) + (p[1] - t[1]) * (p[1] - t[1])
                + (p[2] - t[2]) * (p[2] - t[2]);
    if (d2 < radius * radius)
    {
        assert(!region.gid.empty());
        overlap_regions_.push_back(region);
        overlaps = true;
    }
    return overlaps;
}

void LocalizationRegions::computeOverlapRegionsGlobal()
{
    Control& ct = *(Control::instance());

    if (ct.verbose > 0)
        printWithTimeStamp(
            "LocalizationRegions::computeOverlapRegionsGlobal()...",
            (*MPIdata::sout));

    overlap_regions_.clear();

    vector<LRData>::const_iterator itc = all_regions_.begin();
    while (itc != all_regions_.end())
    {
        assert(!itc->gid.empty());
        const LRData& region(*itc);
        addOverlappingRegion(region);
        itc++;
    }

    // cout<<"Num. overalapping regions: "<<overlap_regions_.size()<<endl;
    if (ct.verbose > 0)
        printWithTimeStamp(
            "LocalizationRegions::computeOverlapRegionsGlobal() done...",
            (*MPIdata::sout));

    // initialize old centers
    resetOldCenters();
}

void LocalizationRegions::computeSubdivOverlapGids()
{
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    const short subdivx = mymesh->subdivx();
    subdiv_overlap_gids_.resize(subdivx);

    const double loc_length = subdom_lattice_[0] / subdivx;
    const double hgrid0     = mygrid.hgrid(0);

    for (short iloc = 0; iloc < subdivx; iloc++)
        subdiv_overlap_gids_[iloc].clear();

    for (vector<LRData>::const_iterator it = overlap_regions_.begin();
         it != overlap_regions_.end(); it++)
    {
        const Vector3D& center(it->center);

        const float radius = it->radius + subdom_delta_;
        assert(radius > 0.1);
        for (short iloc = 0; iloc < subdivx; iloc++)
        {
            // subdiv center
            double cc[3] = { 0.5
                                 * (2. * subdom_lower_left_[0]
                                     + (2 * iloc + 1) * loc_length),
                0.5 * (2. * subdom_lower_left_[1] + subdom_div_lattice_[1]),
                0.5 * (2. * subdom_lower_left_[2] + subdom_div_lattice_[2]) };
            double t[3]  = { center[0], center[1], center[2] };
            for (short i = 0; i < 3; i++)
            {
                while (t[i] >= cc[i] + 0.5 * mygrid.ll(i))
                    t[i] -= mygrid.ll(i);
                while (t[i] < cc[i] - 0.5 * mygrid.ll(i))
                    t[i] += mygrid.ll(i);

                t[i] -= subdom_lower_left_[i];
            }

            // find "p" closest point within subdomain (projection)
            double p[3] = { t[0], t[1], t[2] };
            if (p[0] < iloc * loc_length) p[0] = iloc * loc_length;
            if (p[1] < 0.) p[1] = 0.;
            if (p[2] < 0.) p[2] = 0.;
            if (p[0] > (iloc + 1) * loc_length - hgrid0)
                p[0] = (iloc + 1) * loc_length - hgrid0;
            if (p[1] > subdom_div_lattice_[1]) p[1] = subdom_div_lattice_[1];
            if (p[2] > subdom_div_lattice_[2]) p[2] = subdom_div_lattice_[2];
            double d2 = (p[0] - t[0]) * (p[0] - t[0])
                        + (p[1] - t[1]) * (p[1] - t[1])
                        + (p[2] - t[2]) * (p[2] - t[2]);
            if (d2 < radius * radius)
            {
                for (vector<int>::const_iterator ii = it->gid.begin();
                     ii != it->gid.end(); ii++)
                {
                    subdiv_overlap_gids_[iloc].push_back(*ii);
                }
            }
        }
    }
}

bool LocalizationRegions::overlap(const int gid1, const int gid2)
{
    assert(subdiv_overlap_gids_.size() > 0);

    Mesh* mymesh        = Mesh::instance();
    const short subdivx = mymesh->subdivx();

    for (short iloc = 0; iloc < subdivx; iloc++)
    {
        const vector<int>& ov(subdiv_overlap_gids_[iloc]);
        // assert( ov.size()>0 );
        bool f1 = (find(ov.begin(), ov.end(), gid1) != ov.end());
        bool f2 = (find(ov.begin(), ov.end(), gid2) != ov.end());
        if (f1 && f2) return true;
    }
    return false;
}

bool LocalizationRegions::overlapSubdiv(const int gid, const short iloc) const
{
    assert(subdiv_overlap_gids_.size() > 0);

    const vector<int>& ov(subdiv_overlap_gids_[iloc]);
    return (find(ov.begin(), ov.end(), gid) != ov.end());
}

bool LocalizationRegions::overlapSubdiv(const int gid) const
{
    assert(subdiv_overlap_gids_.size() > 0);

    Mesh* mymesh        = Mesh::instance();
    const short subdivx = mymesh->subdivx();

    short count = 0;
    for (short iloc = 0; iloc < subdivx; iloc++)
    {
        count += overlapSubdiv(gid, iloc);
    }
    return (count > 0);
}

void LocalizationRegions::computeOverlapGids()
{
    overlap_gids_.clear();
    vector<LRData>::const_iterator it = overlap_regions_.begin();
    while (it != overlap_regions_.end())
    {
        assert(!it->gid.empty());
        vector<int>::const_iterator git = it->gid.begin();
        while (git != it->gid.end())
        {
            overlap_gids_.push_back(*git);
            git++;
        }
        it++;
    }
}

void LocalizationRegions::getGidsGlobal(std::vector<int>& gids)
{
    syncCenters();

    gids.clear();
    gids.reserve(128);

    vector<LRData>::const_iterator it = all_regions_.begin();
    while (it != all_regions_.end())
    {
        vector<int>::const_iterator git = it->gid.begin();
        while (git != it->gid.end())
        {
            gids.push_back(*git);
            git++;
        }
        it++;
    }
}

void LocalizationRegions::syncCenters()
{
    syncCenter_tm_.start();
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    Control& ct = *(Control::instance());
    if (ct.verbose > 0)
        printWithTimeStamp(
            "LocalizationRegions::syncCenters()...", (*MPIdata::sout));

    vector<int> local_indexes;
    vector<double> local_centers;
    vector<float> local_radius;

    int num_gid = 0;
    for (vector<LRData>::const_iterator it = local_regions_.begin();
         it != local_regions_.end(); ++it)
    {
        assert(!it->gid.empty());
        assert(it->gid[0] >= 0);

        local_centers.push_back(it->center[0]);
        local_centers.push_back(it->center[1]);
        local_centers.push_back(it->center[2]);
        local_indexes.push_back(it->gid[0]); /// assuming only one gid/center
        num_gid++;
        local_radius.push_back(it->radius);
    }
    // cout<<"LocalizationRegions::syncCenters() -- Nb. gids = "<<num_gid
    //    <<" on PE "<<mmpi.mype()<<endl;

    mmpi.allreduce(&num_gid, 1, MPI_SUM);
    assert(num_gid > 0);
    if (onpe0 && ct.verbose > 0)
        cout << "LocalizationRegions::syncCenters() -- Total Nb. gids = "
             << num_gid << endl;

    vector<double> centers(3 * num_gid);
    vector<int> indexes(num_gid);
    vector<float> radius(num_gid);

    mmpi.allGatherV(local_centers, centers);
    mmpi.allGatherV(local_indexes, indexes);
    mmpi.allGatherV(local_radius, radius);

    all_regions_.resize(num_gid);
    nglobal_ = (int)all_regions_.size();

    int j = 0;
    for (int i = 0; i < nglobal_; i++)
    {
        assert(radius[i] > 0.001);
        assert(indexes[i] < num_gid);

        LRData& region(all_regions_[indexes[i]]);

        region.center.assign(centers[j], centers[j + 1], centers[j + 2]);
        region.radius = radius[i];
        region.gid.clear();
        region.gid.push_back(indexes[i]);

        j += 3;
    }

    for (vector<LRData>::const_iterator it = all_regions_.begin();
         it != all_regions_.end(); ++it)
    {
        assert(!it->gid.empty());
    }

    /* compute volume */
    computeVolume();

    if (ct.verbose > 0)
        printWithTimeStamp(
            "LocalizationRegions::syncCenters() done...", (*MPIdata::sout));
    syncCenter_tm_.stop();
}

/* update local overlap regions by considering neighboring information */
/* NOTE!!!! We assume function multiplicity is 1 !!!! */
void LocalizationRegions::updateOverlapRegions()
{
    /* first gather list of existing local data */
    const int init_size = overlap_regions_.size();
    set<int> overlap_gid_set;
    vector<LRData>::iterator it = overlap_regions_.begin();
    while (it != overlap_regions_.end())
    {
        vector<int>::iterator git = it->gid.begin();
        while (git != it->gid.end())
        {
            overlap_gid_set.insert(*git);
            assert(*git >= 0);
            assert(*git <= 1e6);
            git++;
        }
        it++;
    }

    /* loop over the directions to gather data */
    for (short dir = 0; dir < 3; dir++)
    {
        /* send data to the left and recv from right */
        /* then, send data to the right and recv from left */
        for (short disp = -1; disp < 2; disp = disp + 2)
        {
            /* copy current overlap_regions data */
            vector<LRData> ldata(overlap_regions_);

            /* Compute some data sizes */
            int ndata    = ldata.size();
            data_pos_[0] = 0;
            data_pos_[1] = (ndata + 1) * sizeof(int);
            int offset   = data_pos_[1] + ndata * sizeof(float);
            int padding  = (offset % sizeof(double));
            data_pos_[2] = offset + padding;
            const int bsiz
                = data_pos_[2]
                  + 3 * (old_centers_.size() + 2) * ndata * sizeof(double);
            augmentOverlapRegions(dir, disp, bsiz, ldata, overlap_gid_set);
        }
    }

    /* Now remove non-overlapping functions */
    vector<LRData> regions_buffer(overlap_regions_);
    overlap_regions_.clear();

    vector<vector<Vector3D>> oc;

    for (unsigned int i = 0; i < old_centers_.size(); i++)
    {
        oc.push_back(old_centers_[i]);
        old_centers_[i].clear();
    }

    int pos = 0;

    // add back "local" regions if they still overlap
    it                           = regions_buffer.begin();
    vector<LRData>::iterator it2 = regions_buffer.begin() + init_size;
    while (it != it2)
    {
        bool ovlps = addOverlappingRegion(*it);

        if (ovlps)
        {
            for (unsigned int i = 0; i < old_centers_.size(); i++)
            {
                old_centers_[i].push_back(oc[i][pos]);
            }
        }

        it++;
        pos++;
    }

    // add new regions
    for (it = it2; it != regions_buffer.end(); it++)
    {
        assert(!it->gid.empty());
        overlap_regions_.push_back(*it);
    }

    for (unsigned int i = 0; i < old_centers_.size(); i++)
    {
        for (vector<Vector3D>::const_iterator itc = oc[i].begin() + init_size;
             itc != oc[i].end(); ++itc)
        {
            old_centers_[i].push_back(*itc);
        }
    }

    /* Now update overlap_gids_ and subdiv_overlap_gids_ members */
    computeOverlapGids();

    computeSubdivOverlapGids();
}

void LocalizationRegions::augmentOverlapRegions(const short dir,
    const short disp, const int buffer_size, const vector<LRData>& ldata,
    set<int>& overlap_gid_set)
{
    int lsize                = ldata.size();
    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();
    MPI_Comm cart_comm       = myPEenv.cart_comm();

    /*setup communication buffers */
    // int bsiz = 2*buffer_size;
    int bsiz  = buffer_size;
    char* buf = new char[bsiz];

    /* pack overlap_regions_ data */
    /* size of overlap_regions */
    assert(data_pos_[0] >= 0);
    assert(data_pos_[0] < bsiz);
    int* iptr = (int*)&buf[data_pos_[0]];
    *(iptr++) = lsize;

    /* pack overlap gids */
    for (vector<LRData>::const_iterator it = ldata.begin(); it != ldata.end();
         it++)
        *(iptr++) = it->gid[0];

    /* pack overlap radii */
    assert(data_pos_[1] < bsiz);
    float* fptr = (float*)&buf[data_pos_[1]];
    for (vector<LRData>::const_iterator it = ldata.begin(); it != ldata.end();
         it++)
        *(fptr++) = it->radius;

    /* pack overlap centers */
    assert(data_pos_[2] <= bsiz); // can be equal for n==0
    double* dptr = (double*)&buf[data_pos_[2]];
    int ndoubles = 0;
    for (vector<LRData>::const_iterator it = ldata.begin(); it != ldata.end();
         it++)
    {
        for (short i = 0; i < 3; i++)
            *(dptr++) = it->center[i];
        ndoubles += 3;
    }

    vector<Vector3D>::iterator it;

    /* pack old centers */
    for (unsigned int i = 0; i < old_centers_.size(); i++)
    {
        for (it = old_centers_[i].begin(); it != old_centers_[i].end(); it++)
        {
            for (short j = 0; j < 3; j++)
                *(dptr++) = (*it)[j];
            ndoubles += 3;
        }
    }

    /* Done packing data */
    int siz = data_pos_[2] + ndoubles * sizeof(double);
    assert(siz <= bsiz);

    /* Begin data distribution */
    MPI_Request request[2] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL };

    /* Get source and destination ID to send and recv data */
    int source, dest;
    MPI_Cart_shift(cart_comm, dir, disp, &source, &dest);

    /* step one time to get data from nearest neighbor in dir dirction */
    /* Send and receive size of data */
    int rsiz = 0;
    MPI_Irecv(&rsiz, 1, MPI_INT, source, 0, cart_comm, &request[0]);
    MPI_Isend(&siz, 1, MPI_INT, dest, 0, cart_comm, &request[1]);
    MPI_Waitall(2, request, MPI_STATUS_IGNORE);

    /* Send and receive data */
    char* rbuf              = new char[rsiz];
    MPI_Request request2[2] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL };
    MPI_Irecv(rbuf, rsiz, MPI_CHAR, source, 0, cart_comm, &request2[0]);
    MPI_Isend(buf, siz, MPI_CHAR, dest, 0, cart_comm, &request2[1]);
    MPI_Waitall(2, request2, MPI_STATUS_IGNORE);

    /* set pointer positions */
    int* const recv_size = (int*)rbuf;
    const int nrecv      = *recv_size;
    assert(nrecv >= 0);
    assert(nrecv < 10000);
    int* const recv_gids    = recv_size + 1;
    float* const recv_radii = (float*)(rbuf + (nrecv + 1) * sizeof(int));
    const int offset        = (nrecv + 1) * sizeof(int) + nrecv * sizeof(float);
    const int padding       = offset % sizeof(double);
    double* const recv_centers = (double*)(rbuf + offset + padding);
    assert(offset + padding + 6 * nrecv * sizeof(double)
           <= static_cast<unsigned int>(rsiz));
    std::vector<double*> old_centers(old_centers_.size());
    for (unsigned int i = 0; i < old_centers_.size(); i++)
    {
        old_centers[i]
            = (double*)(rbuf + offset + padding
                        + 3 * (i + 1) * nrecv
                              * sizeof(double)); // recv_centers + 3*nrecv;
    }

    /* loop over recv'd data and augment overlap_regions_ */
    for (int i = 0; i < nrecv; i++)
    {
        // if function already overlaps on local subdomain, ignore
        const int rgid = recv_gids[i];
        assert(rgid >= 0);
        assert(rgid <= 1e6);
        if (overlap_gid_set.find(rgid) != overlap_gid_set.end()) continue;

        /* update map */
        overlap_gid_set.insert(rgid);

        // include potentially new overlap
        const int idx = 3 * i;
        assert(offset + padding + (idx + 2) * sizeof(double)
               < static_cast<unsigned int>(rsiz));
        Vector3D region_center(
            recv_centers[idx], recv_centers[idx + 1], recv_centers[idx + 2]);
        LRData region;
        region.center = region_center;
        region.radius = recv_radii[i];
        vector<int> new_func(1, recv_gids[i]);
        region.gid    = new_func;
        bool overlaps = addOverlappingRegion(region);

        /* update old centers */
        vector<Vector3D> oc;
        for (unsigned int i = 0; i < old_centers_.size(); i++)
        {
            assert(offset + padding + 3 * (i + 1) * nrecv * sizeof(double)
                       + (idx + 2) * sizeof(double)
                   < static_cast<unsigned int>(rsiz));
            oc.push_back(Vector3D(old_centers[i][idx], old_centers[i][idx + 1],
                old_centers[i][idx + 2]));
            if (overlaps) old_centers_[i].push_back(oc[i]);
        }
    }

    delete[] buf;
    delete[] rbuf;
}

void LocalizationRegions::printTimers(ostream& os)
{
    setupLocalRegionsFromOverlapRegions_tm_.print(os);
    syncCenter_tm_.print(os);
}

void LocalizationRegions::getMatrixDistances(
    SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& mat,
    const vector<vector<int>>& gids)
{
    Mesh* mymesh        = Mesh::instance();
    const short subdivx = mymesh->subdivx();
    mat.reset();

    for (short iloc = 0; iloc < subdivx; iloc++)
    {
        const vector<int>& gids_iloc(gids[iloc]);
        const int n = gids_iloc.size();
        std::vector<double> tmp(n * n);
        memset(tmp.data(), 0, n * n * sizeof(double));
        short i = 0;
        for (vector<int>::const_iterator it1 = gids_iloc.begin();
             it1 != gids_iloc.end(); ++it1)
        {
            if ((*it1) >= 0)
            {
                short j = 0;
                for (vector<int>::const_iterator it2 = gids_iloc.begin();
                     it2 != gids_iloc.end(); ++it2)
                {
                    if ((*it2) >= 0)
                        tmp[i + j * n] = sqrt(getDistance2(*it1, *it2));
                    j++;
                }
            }
            i++;
        }
        mat.setValues(tmp.data(), n, iloc);
    }
}

void LocalizationRegions::push_back_local(
    const Vector3D& center, const float radius, const std::vector<int>& gid)
{
    assert(radius > 0.);
    assert(gid.size() > 0);

    LRData lr;
    lr.center = center;
    lr.radius = radius;
    lr.gid    = gid;

    if (isLRcenterLocal(lr))
    {
        // cout<<"Added region with gid "<<gid[0]<<std::endl;
        local_regions_.push_back(lr);

        volume_
            += fourthirdpi_ * (radius * radius * radius) * (float)gid.size();

        assert(volume_ > 0.);
    }
}

void LocalizationRegions::push_back_global(
    const Vector3D& center, const float radius, const std::vector<int>& gid)
{
    assert(radius > 0.1);
    assert(gid.size() > 0);

    LRData lr;
    lr.center = center;
    lr.radius = radius;
    lr.gid    = gid;

    all_regions_.push_back(lr);
    nglobal_ = all_regions_.size();

    volume_ += fourthirdpi_ * (radius * radius * radius) * (float)gid.size();

    assert(volume_ > 0.);
}

// random generator function:
// This function shall return a non-negative value less than its argument.
// rand() returns a pseudo-random integral number in the range between 0
// and RAND_MAX.
static int myrandom(int i) { return std::rand() % i; }

void LocalizationRegions::randomizeGids()
{
    /* initialize random seed: */
    srand(1234);

    random_shuffle(all_regions_.begin(), all_regions_.end(), myrandom);

    int new_gid = 0;
    for (std::vector<LRData>::iterator it = all_regions_.begin();
         it != all_regions_.end(); ++it)
    {
        vector<int>& gid = it->gid;
        for (vector<int>::iterator git = gid.begin(); git != gid.end(); ++git)
        {
            *git = new_gid;
            new_gid++;
        }
    }
}

// get distance between 2 orbitals centers
double LocalizationRegions::distanceBetweenCenters(
    const int gid1, const int gid2) const
{
    Control& ct            = *(Control::instance());
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    Vector3D ll(mygrid.ll(0), mygrid.ll(1), mygrid.ll(2));

    const Vector3D& center1(getCenter(gid1));
    const Vector3D& center2(getCenter(gid2));

    double d = center1.minimage(center2, ll, ct.bcPoisson);

    return d;
}

void LocalizationRegions::writeOldCenters(HDFrestart& h5f_file)
{
    Control& ct(*(Control::instance()));

    if (onpe0 && ct.verbose > 1)
    {
        (*MPIdata::sout) << "LocalizationRegions::writeOldCenters" << endl;
    }

    vector<int> gids;
    for (vector<LRData>::const_iterator it = overlap_regions_.begin();
         it != overlap_regions_.end(); ++it)
    {
        gids.push_back(it->gid[0]);
    }

    hid_t file_id = h5f_file.file_id();
    if (file_id >= 0)
    {
        string datasetname = "GidsList";

        // fill up data array to dimension common to all tasks
#ifdef MGMOL_USE_HDF5P
        if (h5f_file.useHdf5p())
        {
            short s = gids.size();
            short ms;
            mgmol_tools::allreduce(&s, &ms, 1, MPI_MAX, h5f_file.comm_active());
            for (short i = s; i < ms; i++)
                gids.push_back(-1);

            cerr << "writeOldCenters(), parallelWrite1d not yet implemented."
                 << endl;
        }
        else
#endif
        {
            mgmol_tools::write1d(file_id, datasetname, gids, gids.size());
        }
    }

    h5f_file.add2File(num_extrapolations_, "NumOldCenters");

    for (int i = 0; i < num_extrapolations_; i++)
    {
        writeOldCenter(h5f_file, i);
    }
}

void LocalizationRegions::writeOldCenter(HDFrestart& h5f_file, int i)
{
    assert(overlap_regions_.size() == old_centers_[i].size());

    vector<double> data;
    vector<Vector3D>::const_iterator center = old_centers_[i].begin();

    while (center != old_centers_[i].end())
    {
        data.push_back((*center)[0]);
        data.push_back((*center)[1]);
        data.push_back((*center)[2]);
        center++;
    }

    hid_t file_id = h5f_file.file_id();
    if (file_id >= 0)
    {
        stringstream datasetname;
        datasetname << "OldCenter_" << i;

#ifdef MGMOL_USE_HDF5P
        // fill up data array to dimension common to all tasks
        if (h5f_file.useHdf5p())
        {
            short s = data.size();
            short ms;
            mgmol_tools::allreduce(&s, &ms, 1, MPI_MAX, h5f_file.comm_active());
            for (short j = s; j < ms; j++)
                data.push_back(1.e32);

            size_t dims[2] = { data.size() / 3, 3 };

            mgmol_tools::parallelWrite2d(
                file_id, datasetname.str(), data, dims, h5f_file.comm_active());
        }
        else
#endif
        {
            size_t dims[2] = { data.size() / 3, 3 };
            mgmol_tools::write2d(file_id, datasetname.str(), data, dims);
        }
    }
}

void LocalizationRegions::setupOldCenters(HDFrestart& h5_file)
{
    vector<int> gids;
    std::string datasetname("GidsList");
    h5_file.readAtomicData(datasetname, gids);

    map<int, int> gids_map;
    for (unsigned int i = 0; i < gids.size(); i++)
    {
        gids_map[gids[i]] = i;
    }

    vector<double> data;
    vector<Vector3D>::iterator center;

    num_extrapolations_ = h5_file.getFromFile("NumOldCenters");

    // cout<<num_extrapolations_<<endl;
    // cout<<old_centers_set_<<endl;

    for (int i = 0; i < num_extrapolations_; i++)
    {
        data.clear();
        h5_file.readOldCenter(data, i);

        center = old_centers_[i].begin();

        int idx;
        int gid;
        for (unsigned int j = 0; j < overlap_regions_.size(); j++)
        {
            gid = overlap_regions_[j].gid[0];
            idx = gids_map[gid];
            center->assign(data[3 * idx], data[3 * idx + 1], data[3 * idx + 2]);
            center++;
        }
    }
}

short LocalizationRegions::getNumOldCenters() const
{

    return num_extrapolations_;
}

double LocalizationRegions::computeMinDistBetweenLocalPairs(
    std::ostream& os, const bool print) const
{
    Control& ct            = *(Control::instance());
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    Vector3D ll(mygrid.ll(0), mygrid.ll(1), mygrid.ll(2));

    double distance = 1.e18;

    for (vector<LRData>::const_iterator it1 = local_regions_.begin();
         it1 != local_regions_.end(); it1++)
        for (vector<LRData>::const_iterator it2 = local_regions_.begin();
             it2 != local_regions_.end(); it2++)
        {
            if (it1 != it2)
            {
                double d
                    = (it1->center).minimage(it2->center, ll, ct.bcPoisson);
                if (d < distance) distance = d;
            }
        }

    if (distance < 1.e-3 && print)
    {
        for (vector<LRData>::const_iterator it1 = local_regions_.begin();
             it1 != local_regions_.end(); it1++)
            os << "center at " << it1->center << endl;
    }

    return distance;
}

// void LocalizationRegions::setOldCenters(const vector<double>& tau)
//{
//    assert(tau.size() == 3*local_ions_.size());
//
//    vector<Ion*>::iterator ion=local_ions_.begin();
//    int ia=0;
//    while(ion!=local_ions_.end()){
//        (*ion)->setPosition(tau[3*ia+0],
//                          tau[3*ia+1],
//                          tau[3*ia+2]);
//
//        ion++;
//        ia++;
//    }
//
//    setup_=false;
//}
//
// void setOldCenter(const double x, const double y, const double z)
//{
//    old_position_[0]=position_[0];
//    old_position_[1]=position_[1];
//    old_position_[2]=position_[2];
//
//    position_[0]=x;
//    position_[1]=y;
//    position_[2]=z;
//
//    kbproj_.clear();
//}

template float LocalizationRegions::move(
    const SpreadsAndCenters<LocGridOrbitals<ORBDTYPE>>& sc, const bool flag);
template float LocalizationRegions::updateRadiiConstVol(
    const SpreadsAndCenters<LocGridOrbitals<ORBDTYPE>>& sc);
template float LocalizationRegions::updateRadii(
    const SpreadsAndCenters<LocGridOrbitals<ORBDTYPE>>& sc, const float ratio);
