// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_LOCALIZATIONREGIONS_H
#define MGMOL_LOCALIZATIONREGIONS_H

#include "Control.h"
#include "HDFrestart.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"
#include "SpreadsAndCenters.h"
#include "SquareLocalMatrices.h"
#include "Timer.h"
#include "Vector3D.h"
#include "global.h"
#include "mgmol_mpi_tools.h"
#include "tools.h"

#ifdef MGMOL_USE_SCALAPACK
#include "OrbitalsTransform.h"
#endif

#include <iostream>
#include <set>
#include <vector>

class SymmetricPair;

typedef struct LRData
{
    Vector3D center;
    float radius;
    std::vector<int> gid;
} LRData;

// container class for localization regions (centers and radii)

class LocalizationRegions
{

private:
    static Timer setupLocalRegionsFromOverlapRegions_tm_;
    static Timer syncCenter_tm_;

    static double fourthirdpi_;

    bool gid_local_regions_set_;

    std::vector<int> overlap_gids_;
    std::vector<std::vector<int>> subdiv_overlap_gids_;

    int nglobal_;

    float max_displ_;

    void bcastLRs();
    void computeOverlapGids();
    void computeOverlapRegionsGlobal();
    void computeSubdivOverlapGids();
    bool isLRcenterLocal(const LRData& lr) const;
    void syncCenters();

    float getMeanRadius();

    float computeVolumePrivate()
    {
        float volume = 0.;
        for (std::vector<LRData>::iterator it = local_regions_.begin();
             it != local_regions_.end(); it++)
        {
            const float alpha = it->radius;
            volume += (alpha * alpha * alpha);
        }

        volume *= fourthirdpi_;

        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        mmpi.allreduce(&volume, 1, MPI_SUM);

        return volume;
    }
    void augmentOverlapRegions(const short dir, const short disp,
        const int buffer_size, const std::vector<LRData>& aug_data,
        std::set<int>& overlap_gid_map);
    bool addOverlappingRegion(const LRData& region);

    // local subdomain geometry
    double subdom_origin_[3];
    double subdom_end_[3];
    double subdom_lower_left_[3];
    double subdom_upper_right_[3];
    double subdom_center_[3];
    double subdom_div_lattice_[3];
    double subdom_lattice_[3];
    double subdom_delta_;
    void setupSubdomainGeometry();

    //  private members for data distribution of LRData
    int data_pos_[3]; // starting position of ints (gids), floats (radii) and
                      // doubles (centers)

    void push_back_global(const Vector3D& center, const float radius,
        const std::vector<int>& gid);

    void push_back_local(const Vector3D& center, const float radius,
        const std::vector<int>& gid);

    void writeOldCenter(HDFrestart& h5f_file, int i);
    void writeOldCenterOnMesh(HDFrestart& h5f_file, int i);

    void shiftCenters(Vector3D dir);

    void shiftCenters(std::vector<Vector3D> shifts);

    unsigned int iterative_index_;

protected:
    std::vector<LRData> all_regions_; // global

    float volume_; // total volume of LRs

    std::vector<LRData>
        local_regions_; // center of LR located in local subdomain

    bool old_centers_set_;

    std::vector<LRData> overlap_regions_; // LR overlaps with local subdomain

    Vector3D cell_;

    float tol_displ_;

    std::vector<std::vector<Vector3D>> old_centers_; // old centers of LRs

    short K_; // number of old centers to store

    int num_extrapolations_;

    void setupLocalRegionsFromOverlapRegions();

public:
    LocalizationRegions(const Vector3D& cell, const float tol_displ = 1.e10)
        : gid_local_regions_set_(false),
          nglobal_(-1),
          volume_(0.),
          cell_(cell),
          tol_displ_(tol_displ)
    {
        assert(cell_[0] > 0.01);
        assert(cell_[1] > 0.01);
        assert(cell_[2] > 0.01);

        setupSubdomainGeometry();

        max_displ_ = 0.;

        K_ = 1;
        if (K_ == 1) K_++;

        int len_old_centers = (K_ < 3) ? 3 : K_;

        old_centers_.resize(len_old_centers);

        old_centers_set_    = false;
        num_extrapolations_ = 0;

        iterative_index_ = 0;
    }
    void printTimers(std::ostream& os);

    void writeOldCenters(HDFrestart& h5f_file);
    void setupOldCenters(HDFrestart& h5_file);

    int globalNumLRs() const { return nglobal_; }

    virtual ~LocalizationRegions() { }

    virtual void setup()
    {
        Control& ct = *(Control::instance());

        if (ct.verbose > 0)
            printWithTimeStamp(
                "LocalizationRegions::setup()...", (*MPIdata::sout));

        clearOldCenters();

        if (ct.restart_info <= 2 || !ct.isLocMode())
        {
            // bcast data regions_ read from PE0
            bcastLRs();
            computeOverlapRegionsGlobal(); // compute overlaps before setting up
                                           // local regions
            setupLocalRegionsFromOverlapRegions();
        }
        else
        {
            // Use local_regions_ to setup regions_
            syncCenters();
            computeOverlapRegionsGlobal();
        }

        computeOverlapGids();

        computeSubdivOverlapGids();
    }

    void printInfo(std::ostream& os)
    {
        Control& ct = *(Control::instance());
        if (onpe0 && ct.verbose > 0)
        {
            os << "LocalizationRegions: Number of localization regions = "
               << nglobal_ << std::endl;
        }
    }

    float moveTo(const std::vector<Vector3D>& target_centers);

    bool moveIsSmall() { return (max_displ_ < tol_displ_); }

    bool hasNcentersReachedNumber(const int max_nb_lrs)
    {
        return (nglobal_ >= max_nb_lrs);
    }

    const Vector3D& getCenter(const int gid) const
    {
        assert(0 <= gid);
        assert(gid < nglobal_);

        for (std::vector<LRData>::const_iterator it = overlap_regions_.begin();
             it != overlap_regions_.end(); ++it)
        {
            const std::vector<int>& gids(it->gid);
            if (gid == gids[0]) return it->center;
        }

        std::cout << "mype:" << mype
                  << ", WARNING: didn't find a center for gid=" << gid
                  << std::endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        return (overlap_regions_.end())->center;
    }

    int getLRsColor(const int gid) const
    {
        assert(0 <= gid);
        assert(gid < nglobal_);

        const short iloc = 0;

        for (unsigned int i = 0; i < subdiv_overlap_gids_[iloc].size(); i++)
        {
            if (subdiv_overlap_gids_[iloc][i] == gid) return i;
        }

        std::cout << "mype:" << mype
                  << ", WARNING: didn't find a color for gid=" << gid
                  << std::endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        return -1;
    }

    double getCenterAndRadius(const int gid, Vector3D& center) const
    {
        assert(0 <= gid);
        assert(gid < nglobal_);

        for (std::vector<LRData>::const_iterator it = overlap_regions_.begin();
             it != overlap_regions_.end(); ++it)
        {
            const std::vector<int>& gids(it->gid);
            assert(!gids.empty());
            if (gid == gids[0])
            {
                center = it->center;
                return it->radius;
            }
        }
        std::cout << "mype:" << mype
                  << ", WARNING: didn't find a center and radius for gid="
                  << gid << std::endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

        center = (overlap_regions_.end())->center;
        return (overlap_regions_.end())->radius;
    }

    bool checkCenter(const Vector3D& center) const
    {
        assert(fabs(center[0]) < 1.e30);
        assert(fabs(center[1]) < 1.e30);
        assert(fabs(center[2]) < 1.e30);
        assert(overlap_regions_.size() > 0);

        bool check = false;
        for (std::vector<LRData>::const_iterator it = overlap_regions_.begin();
             it != overlap_regions_.end(); ++it)
        {
            check = (check || center == it->center);
        }

        return check;
    }

    void setCenters(const std::vector<double>& tau)
    {
        const int nc = (int)overlap_regions_.size();
        for (int j = 0; j < nc; j++)
        {
            overlap_regions_[j].center.assign(
                tau[3 * j], tau[3 * j + 1], tau[3 * j + 2]);
        }

        setupLocalRegionsFromOverlapRegions();
    }

    void fillWithZeroCenters(const int n)
    {
        Vector3D tmp_center(0., 0., 0.);
        const int istart = (int)all_regions_.size();
        for (int i = istart; i < n; i++)
            push_back_global(tmp_center, 10000.);

        nglobal_ = (int)all_regions_.size();
    }

    void setCenter(const int gid, const Vector3D& center)
    {
        assert(gid < (int)all_regions_.size());

        for (std::vector<LRData>::iterator it = overlap_regions_.begin();
             it != overlap_regions_.end(); ++it)
        {
            const std::vector<int>& gids(it->gid);
            if (gid == gids[0]) it->center = center;
        }
    }

    // return radius of regions gid
    // return -1 if gid not known locally
    float radius(const int gid) const
    {
        assert(gid < nglobal_);
        for (std::vector<LRData>::const_iterator it = overlap_regions_.begin();
             it != overlap_regions_.end(); ++it)
        {
            const std::vector<int>& gids(it->gid);
            if (gid == gids[0]) return it->radius;
        }

        return -1.;
    }

    float volume() const { return volume_; }

    // create a new function with a new gid
    void push_back_global(
        const Vector3D& center, const float radius, const int gid = -1)
    {
        assert(radius > 0.);

        const int new_gid = (gid == -1) ? (int)all_regions_.size() : gid;
        std::vector<int> gids(1);
        gids[0] = new_gid;

        push_back_global(center, radius, gids);
    }

    void push_back_local(
        const Vector3D& center, const float radius, const int gid = -1)
    {
        assert(radius > 0.);

        const int new_gid = (gid == -1) ? (int)local_regions_.size() : gid;
        std::vector<int> gids(1);
        gids[0] = new_gid;

        push_back_local(center, radius, gids);
    }

    // compute distance square between two centers of global indexes gid0 and
    // gid1
    double getDistance2(const int gid0, const int gid1)
    {
        Control& ct = *(Control::instance());

        const Vector3D& center0(getCenter(gid0));
        const Vector3D& center1(getCenter(gid1));

        Vector3D dr = center1.vminimage(center0, cell_, ct.bcPoisson);
        return norm2(dr);
    }

    template <class T>
    float move(const SpreadsAndCenters<T>& sc, const bool flag = false);
#ifdef MGMOL_USE_SCALAPACK
    float updateRadii(const OrbitalsTransform* ot, const float ratio);
#endif
    template <class T>
    float updateRadii(const SpreadsAndCenters<T>& sc, const float ratio);
    template <class T>
    float updateRadiiConstVol(const SpreadsAndCenters<T>& sc);

    float resetRadii(const std::vector<double>& new_radii, const float ratio);
    float resetRadiiConstVol(const std::vector<double>& new_radii);
    void setUniformRadii(const float radius);
    float computeVolume();

    void extrapolateCentersLinear(
        const bool force, const bool move_centers = true);
    void extrapolateCentersQuadratic();
    void extrapolateCentersVerlet(
        const bool force, const bool move_centers = true);
    void moveToExtrapolatedCenters();

    void resetOldCenters();

    void clear();
    void clearOldCenters();

    void printAllRegions(std::ostream&);

    float max_radii() const;
    float min_radii() const;
    float imposeMinDistanceBetweenCentersGlobal(const float);
    float getStatesWithClosestCenters(int* const st1, int* const st2,
        const std::set<SymmetricPair>& exclude_set) const;

    void getLocalSubdomainIndices(std::vector<int>& indices);
    const std::vector<int>& getOverlapGids() const { return overlap_gids_; }
    const std::vector<std::vector<int>>& getSubdivOverlapGids() const
    {
        return subdiv_overlap_gids_;
    }
    int getNumOverlapGids() const { return overlap_gids_.size(); }
    // bool gidIncluded(const int gid)
    //{
    //    return(
    //    find(overlap_gids_.begin(),overlap_gids_.end(),gid)!=overlap_gids_.end()
    //    );
    //}

    bool overlap(const int gid1, const int gid2);
    bool overlapSubdiv(const int gid, const short iloc) const;
    bool overlapSubdiv(const int gid) const;

    void getGidsGlobal(std::vector<int>& gids);
    void updateOverlapRegions();
    void getMatrixDistances(
        SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& mat,
        const std::vector<std::vector<int>>& gids);

    void printMinMaxRadius(std::ostream& os)
    {
        const float minr = min_radii();
        const float maxr = max_radii();
        if (onpe0)
        {
            os << "Radius LRs: min = " << minr << ", max = " << maxr
               << std::endl;
        }
    }

    double distanceBetweenCenters(const int gid1, const int gid2) const;

    void randomizeGids();
    void subdomainCenter(double center[3])
    {
        for (short i = 0; i < 3; i++)
            center[i] = subdom_center_[i];
    }
    void subdomainCenter(Vector3D& center)
    {
        center.assign(subdom_center_[0], subdom_center_[1], subdom_center_[2]);
    }

    short getNumOldCenters() const;

    unsigned int getIterativeIndex() const { return iterative_index_; }

    double computeMinDistBetweenLocalPairs(
        std::ostream& os, const bool print) const;
};

#endif
