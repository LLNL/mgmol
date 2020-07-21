// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

/*!
 * Main header file for clustering or distributing orbitals across subdomains
 * for load balancing.
 */
#ifndef CLUSTERORBITALS_H
#define CLUSTERORBITALS_H

#include "DataDistribution.h"
#include "LocalizationRegions.h"
#include "Timer.h"
#include "VariableSizeMatrix.h"
#include "Vector3D.h"
#include <vector>

class ClusterOrbitals
{

private:
    static Timer computeClusters_tm_;
    static Timer setupClusters_tm_;

    std::shared_ptr<LocalizationRegions>
        lrs_; // pointer to localizationRegions object

    bool isswitched_; // check to see if local cluster size crosses average size
    float alpha_; // tuning paramemter for computing bias
    float damping_tol_; // damping parameter for alpha_
    short iters_; // number of iterations
    std::vector<int>
        locfcns_; // cluster of variables assigned to this proc. Initially local
                  // variables centered on this proc.
    double cluster_range_radius_; // radius of how far to gather data to
                                  // determine cluster for each proc.
    int subdom_steps_; // maximum number of steps to gather regions information
                       // from subdomains
    double subdom_bias_; // subdomain bias
    double avg_locfcns_global_; // global average number of centered localized
                                // functions
    double
        avg_locfcns_; // "local" average number of centered localized functions
    double old_avg_locfcns_; // previous "local" average number of centered
                             // localized functions
    int max_cluster_size_; // max. size of cluster in neighborhood

    int subdom_id_; // subdomain id ( == pid)
    short subdom_id_pos_; // position of subdomain info in subdomains data
    Vector3D* subdom_ll_; // subdomain dimension
    Vector3D* subdom_center_; // subdomain center
    std::map<int, Vector3D>
        subdoms_data_; // indexes (pids) and centers of (neighboring) subdomains
    std::vector<int> subdoms_index_; // pids of subdomains
    std::vector<double>
        bias_data_; // vector of bias values subdomain and its neighbors

    std::map<int, Vector3D>
        regions_data_; // gids and centers of regions in clustering search space
    std::vector<int>
        cluster_indexes_; // indexes of regions owned by this subdomain cluster

    DataDistribution* distributor_; // data distribution object
    VariableSizeMatrix<SparseRow>* comm_data_; // communication container

    void packRegionsData();
    void unpackAndSetRegionsData();
    void initializeRegionsData();
    void packSubdomData();
    void unpackAndSetSubdomData();
    void initializeSubdomainData();
    void computeLocalBias(); // compute/ update local bias. Reset initial bias
                             // to zero if needed.
    void updateBiasData();
    double computeSquaredDistanceBetweenCenters(
        const Vector3D& center1, const Vector3D& center2) const;
    bool checkConv();
    void reset();
    void computeLocalRegionOwnership();
    void updateCluster();
    void writeVTKHeader(
        const int npx, const int npy, const int npz, std::ofstream& os);
    void writeVTKDataset(
        const std::string& name, const pb::PEenv& myPEenv, std::ofstream& os);

public:
    ClusterOrbitals(std::shared_ptr<LocalizationRegions> lrs); // constructor
    void setup(); // setup some data structures and subdomain data communication
    int computeClusters(const short
            maxiters); // compute cluster of variables assigned to procs.
    ~ClusterOrbitals(); // destructor
    const std::vector<int>& getClusterIndices() const
    {
        return cluster_indexes_;
    }
    static void printTimers(std::ostream& os) // print timers
    {
        computeClusters_tm_.print(os);
        setupClusters_tm_.print(os);
    }
};
#endif
