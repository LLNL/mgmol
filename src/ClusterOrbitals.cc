// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "ClusterOrbitals.h"
#include "Control.h"
#include "DataDistribution.h"
#include "Mesh.h"

#include <sstream>
#include <vector>

using namespace std;

Timer ClusterOrbitals::computeClusters_tm_("ClusterOrbitals::ComputeClusters");
Timer ClusterOrbitals::setupClusters_tm_("ClusterOrbitals::Setup");

ClusterOrbitals::ClusterOrbitals(std::shared_ptr<LocalizationRegions> lrs)
{
    lrs_ = lrs;

    Control& ct = *(Control::instance());
    // tuning parameter for computing bias
    alpha_              = ct.load_balancing_alpha;
    damping_tol_        = ct.load_balancing_damping_tol;
    iters_              = 0;
    max_cluster_size_   = 0;
    avg_locfcns_        = 0.;
    old_avg_locfcns_    = 0.;
    avg_locfcns_global_ = 0.;

    // subdomain bias
    subdom_bias_ = 0.;

    // cluster subdomain steps (nearest neighbor by default)
    subdom_steps_ = 1;
}

ClusterOrbitals::~ClusterOrbitals()
{
    delete subdom_center_;
    delete subdom_ll_;
    delete distributor_;
    delete comm_data_;
}

void ClusterOrbitals::packRegionsData()
{
    // initialize communication matrix with locally centered data
    comm_data_->setupSparseRows(locfcns_);
    const int count = (int)locfcns_.size();
    // comm_data_ has been initialized with locfcns_ so loop index i
    // corresponds to local index of locfcns_[i]
    for (int i = 0; i < count; i++)
    {
        int cols[3]            = { 0, 1, 2 };
        Vector3D region_center = lrs_->getCenter(locfcns_[i]);
        double vals[3]
            = { (region_center)[0], (region_center)[1], (region_center)[2] };
        comm_data_->updateLocalRowInsert(3, i, &cols[0], &vals[0]);
    }
}

void ClusterOrbitals::unpackAndSetRegionsData()
{
    // unpack regions data from communication matrix (comm_data_)
    std::vector<int> regions_index = comm_data_->lvars();

    const int n = comm_data_->n();
    std::vector<double> center;

    for (int region = 0; region < n; region++)
    {
        comm_data_->getRowEntries(region, center);
        Vector3D region_center(center[0], center[1], center[2]);
        regions_data_[regions_index[region]] = region_center;
    }
}

void ClusterOrbitals::initializeRegionsData()
{
    // Get local functions centered on this processor
    lrs_->getLocalSubdomainIndices(locfcns_);
    int loc_size = (int)locfcns_.size();

    // compute average number of centered (localized) functions
    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();
    avg_locfcns_global_
        = (double)lrs_->globalNumLRs() / (double)myPEenv.n_mpi_tasks();

    for (std::vector<int>::iterator it = locfcns_.begin(); it != locfcns_.end();
         ++it)
    {
        Vector3D region_center = lrs_->getCenter(*it);
        regions_data_[*it]     = region_center;
    }

    // Communicate regions information to initialize some variables
    comm_data_->setupSparseRows(subdoms_index_);
    std::vector<double> vals(loc_size, 0.);
    comm_data_->initializeLocalRow(
        loc_size, subdom_id_pos_, &locfcns_[0], &vals[0]);
    // call updateCluster to communicate cluster info and initialize
    // cluster_indexes
    updateCluster();

    // compute "local average" cluster size
    //   assert(n > 0);
    //   avg_locfcns_ = (double)comm_data_->nnzmat()/ (double)comm_data_->n();
    // compute subdomain alpha
    //   double diff = (double)comm_data_->nzmin() - avg_locfcns_global_;
    //   alpha_ = fabs(diff);
    // initialize switch variable
    isswitched_ = ((double)loc_size < avg_locfcns_global_);

    // initialize cluster to locally centered data.
    // This need not be ordered to match ordering of regions data map object.
    //   cluster_indexes_ = locfcns_;

    // Pack communication matrix with regions data
    //   packRegionsData();

    // Distribute/ Gather regions data
    //   double domain[3] =
    //   {(*subdom_ll_)[0],(*subdom_ll_)[1],(*subdom_ll_)[2]}; const int
    //   max_steps[3] = {subdom_steps_, subdom_steps_, subdom_steps_};
    //   DataDistribution regions_distributor("Init_regions_data",max_steps,
    //   myPEenv, domain); regions_distributor.augmentLocalData(*comm_data_,
    //   true);

    // unpack and initialize
    //   unpackAndSetRegionsData();
}

// compute local bias
void ClusterOrbitals::computeLocalBias()
{
    int sz       = (int)cluster_indexes_.size();
    double denom = sz > 0 ? (double)sz : min(0.5, 0.5 * avg_locfcns_global_);
    double ratio = avg_locfcns_global_ / denom;

    // test switch
    bool test_switch = (denom < avg_locfcns_global_);
    if (isswitched_ != test_switch)
    {
        // modify alpha
        alpha_ = damping_tol_ * alpha_;
        // reset switch
        isswitched_ = test_switch;
    }

    // assume cluster radius to be distance between subdomain centers
    // which is approx. subdom. width. Here, we use cluster_range_radius_
    // variable, which is min(subdom_width) for the different directions.
    subdom_bias_ = subdom_bias_
                   + alpha_ * cluster_range_radius_ * cluster_range_radius_
                         * (pow(ratio, 0.67) - 1);
}

void ClusterOrbitals::packSubdomData()
{
    // reset/ prepare communication matrix
    comm_data_->reset();
    // pack data
    int cols[3]    = { 0, 1, 2 };
    double vals[3] = {};
    for (int k = 0; k < 3; k++)
    {
        vals[k] = (*subdom_center_)[k];
    }
    comm_data_->insertNewRow(3, subdom_id_, &cols[0], &vals[0], true);
}

void ClusterOrbitals::unpackAndSetSubdomData()
{
    // unpack subdomain data from communication matrix
    const int n = comm_data_->n();

    std::vector<int> indexes = comm_data_->lvars();
    std::vector<double> center;
    for (int subdom = 0; subdom < n; subdom++)
    {
        comm_data_->getRowEntries(subdom, center);
        Vector3D subdom_center(center[0], center[1], center[2]);
        subdoms_data_[indexes[subdom]] = subdom_center;
    }

    // save subdomain index ordering and position of local subdomain id
    for (std::map<int, Vector3D>::iterator it = subdoms_data_.begin();
         it != subdoms_data_.end(); ++it)
    {
        subdoms_index_.push_back(it->first);
    }
    subdom_id_pos_
        = std::distance(subdoms_data_.begin(), subdoms_data_.find(subdom_id_));
}

void ClusterOrbitals::initializeSubdomainData()
{
    // cleanup old data
    subdom_bias_ = 0.;
    bias_data_.clear();
    subdoms_index_.clear();
    subdoms_data_.clear();
    // pack communication matrix with subdomain data
    packSubdomData();
    // distribute/ gather subdomain data
    distributor_->augmentLocalData(*comm_data_, true);
    // updack and initialize
    unpackAndSetSubdomData();
}

// Update bias data vector. Get index information from previously
// populated subdomain index vector.
void ClusterOrbitals::updateBiasData()
{
    bias_data_.clear();
    // initialize communication matrix subdomain indexes
    comm_data_->setupSparseRows(subdoms_index_);
    const int count = (int)subdoms_index_.size();
    assert(count > 0);

    // pack communication matrix with bias data
    // Note: We only need to update row corresponding to
    // current subdomain prior to communication.
    // All other rows must be zero.

    const int col = 0;
    const int row = subdom_id_pos_;
    comm_data_->initializeLocalRow(1, row, &col, &subdom_bias_);

    // distribute/ gather bias data to fill in zero rows
    // Note: we do not expect any new rows during communication since number of
    // neighbors is fixed. each subdomain communicates their own uniquely
    // computed bias, hence we can use updateLocalRows.
    distributor_->updateLocalRows(*comm_data_);

    // unpack and update bias data
    std::vector<double> values;
    for (int i = 0; i < count; i++)
    {
        comm_data_->getRowEntries(i, values);
        bias_data_.push_back(values[0]);
    }
}

// get distance between 2 centers
double ClusterOrbitals::computeSquaredDistanceBetweenCenters(
    const Vector3D& center1, const Vector3D& center2) const
{
    Control& ct = *(Control::instance());
    double d    = center1.minimage(center2, *subdom_ll_, ct.bcPoisson);
    return (d * d);
}
// Check convergence
bool ClusterOrbitals::checkConv()
{
    // global convergence criterion
    // compute max cluster size
    //   int max_size = (int)cluster_indexes_.size();
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    //   mmpi.allreduce(&max_size, 1, MPI_MAX);

    // local convergence criterion

    // (strong convergence/ early termination)
    int sconv = (max_cluster_size_ == ceil(avg_locfcns_global_)) ? 0 : 1;
    // (weak convergence)
    double diff = old_avg_locfcns_ - avg_locfcns_;
    int wconv   = (fabs(diff) <= 1.0e-6) ? 0 : 1;
    int conv[2] = { sconv, wconv };

    mmpi.allreduce(&conv[0], 2, MPI_MAX);
    sconv = conv[0];
    wconv = conv[1];
    if (sconv == 0)
    {
        return true;
    }
    else
    {
        // weak convergence
        if ((iters_ > 10) && (wconv == 0))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}
// compute local cluster
int ClusterOrbitals::computeClusters(const short maxiters)
{
    Control& ct     = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    //   mmpi.barrier();
    // timer begin
    computeClusters_tm_.start();

    // reset data
    reset();
    // initialize locally centered regions data
    initializeRegionsData();

    ///*
    // write to vtk file
    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();
    std::ofstream vtkfile;
    std::stringstream ss;
    if (ct.write_clusters && ct.verbose > 1)
    {
        if (onpe0)
            vtkfile.open(ct.load_balancing_output_file.c_str(), std::ios::out);
        writeVTKHeader(myPEenv.n_mpi_task(0), myPEenv.n_mpi_task(1),
            myPEenv.n_mpi_task(2), vtkfile);

        ss << "cluster_data_iter_" << iters_;
        writeVTKDataset(ss.str(), myPEenv, vtkfile);
    }
    //   mmpi.barrier();
    //   writeVTKDataset("initial_clusters",myPEenv,vtkfile);
    //*/
    // Begin
    // Initial stats: max of cluster sizes

    int locmax = (int)locfcns_.size();
    mmpi.allreduce(&locmax, 1, MPI_MAX);
    if (onpe0) cout << "Initial statistics:: max cluster = " << locmax << endl;

    bool conv = false;
    // Do at least 1 iteration
    while (iters_ < maxiters)
    {
        // next iter
        iters_++;

        // setup/ compute subdomain bias data
        computeLocalBias();
        updateBiasData();
        // compute cluster membership for each locally centered region
        computeLocalRegionOwnership();
        // communicate cluster membership information
        updateCluster();

        if (ct.verbose > 1)
        {
            locmax = (int)cluster_indexes_.size();
            mmpi.allreduce(&locmax, 1, MPI_MAX);
            if ((onpe0))
                cout << "iter = " << iters_
                     << ", global max cluster size = " << locmax
                     << ", conv = " << conv << endl;
        }

        // write cluster data
        if (ct.write_clusters && ct.verbose > 1)
        {
            ss.str(std::string());
            ss << "cluster_data_iter_" << iters_;
            writeVTKDataset(ss.str(), myPEenv, vtkfile);
        }

        // check convergence
        conv = checkConv();
        if (conv) break;
    }
    // timer end
    //   mmpi.barrier();
    computeClusters_tm_.stop();

    // Final stats: get min and max of cluster sizes

    locmax = (int)cluster_indexes_.size();
    mmpi.allreduce(&locmax, 1, MPI_MAX);
    if (onpe0)
        cout << "Final Stats: num. iters = " << iters_
             << ", global max cluster = " << locmax << ", conv = " << conv
             << endl;

    ///*
    // write to vtk file
    //   writeVTKDataset("final_clusters",myPEenv,vtkfile);
    if (ct.write_clusters && ct.verbose > 1) vtkfile.close();
    //*/

    // MPI_Finalize();
    // exit(0);

    return iters_;
}

void ClusterOrbitals::setup()
{
    //   MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    //   mmpi.barrier();
    setupClusters_tm_.start();
    // setup some data structures
    // get subdomain center
    subdom_center_ = new Vector3D(0, 0, 0);
    lrs_->subdomainCenter(*subdom_center_);

    // build data distribution object
    Mesh* mymesh             = Mesh::instance();
    const pb::Grid& mygrid   = mymesh->grid();
    const pb::PEenv& myPEenv = mymesh->peenv();
    // assign domain extents
    double domain[3] = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };
    subdom_ll_       = new Vector3D(domain[0], domain[1], domain[2]);

    // compute cluster range radius. Use minimum subdomain width to
    // ensure communication takes only one step (adjacent neighbors only).
    double subdom_width[3] = {};
    for (int k = 0; k < 3; k++)
    {
        subdom_width[k] = domain[k] / myPEenv.n_mpi_task(k);
    }
    cluster_range_radius_  = *std::min_element(subdom_width, subdom_width + 3);
    const int max_steps[3] = { subdom_steps_, subdom_steps_, subdom_steps_ };
    //   distributor_ = new
    //   DataDistribution("computeClusters",cluster_range_radius_, myPEenv,
    //   domain);
    distributor_
        = new DataDistribution("computeClusters", max_steps, myPEenv, domain);

    // save my pid
    subdom_id_ = myPEenv.mytask();

    // construct matrix for data distribution
    comm_data_ = new VariableSizeMatrix<sparserow>("computeClusters", 512);

    // Initialize and communicate subdomain data
    initializeSubdomainData();

    //   mmpi.barrier();
    setupClusters_tm_.stop();
}

void ClusterOrbitals::reset()
{
    // cleanup/ reset old data
    locfcns_.clear();
    cluster_indexes_.clear();
    regions_data_.clear();

    subdom_bias_        = 0.;
    iters_              = 0;
    Control& ct         = *(Control::instance());
    alpha_              = ct.load_balancing_alpha;
    max_cluster_size_   = 0;
    avg_locfcns_        = 0.;
    old_avg_locfcns_    = 0.;
    avg_locfcns_global_ = 0.;
}

// determine which cluster owns which region
void ClusterOrbitals::computeLocalRegionOwnership()
{
    // reset and initialize communication matrix with subdomain info
    comm_data_->setupSparseRows(subdoms_index_);

    std::vector<int> marked_id;
    // clear locally centered list
    locfcns_.clear();
    // loop over locally centered regions
    for (std::map<int, Vector3D>::iterator region = regions_data_.begin();
         region != regions_data_.end(); ++region)
    {
        // compute initial minimum distance.
        std::map<int, Vector3D>::iterator subdom = subdoms_data_.begin();
        double mindist = computeSquaredDistanceBetweenCenters(
                             subdoms_data_[subdom_id_], region->second)
                         - bias_data_[subdom_id_pos_];
        int subdom_index_pos = 0;
        int min_subdom_pos   = subdom_id_pos_;
        // loop over remaining subdomains
        while (subdom != subdoms_data_.end())
        {
            double distance = computeSquaredDistanceBetweenCenters(
                                  subdom->second, region->second)
                              - bias_data_[subdom_index_pos];
            if (distance < mindist)
            {
                min_subdom_pos = subdom_index_pos;
                mindist        = distance;
            }
            subdom++;
            subdom_index_pos++;
        }
        // update subdomain cluster data with current localized region
        comm_data_->updateLocalRowInsert(min_subdom_pos, region->first, 0.);
    }
}

// communicate cluster information
void ClusterOrbitals::updateCluster()
{
    // communicate cluster info
    distributor_->augmentLocalData(*comm_data_, false);

    // unpack local data and update local cluster info
    cluster_indexes_.clear();
    const int lrindex = subdom_id_pos_;
    comm_data_->getColumnIndexes(lrindex, cluster_indexes_);

    max_cluster_size_ = comm_data_->nzmax();
    // compute "local average" cluster size
    //   assert(n > 0);
    old_avg_locfcns_ = avg_locfcns_;
    avg_locfcns_     = (double)comm_data_->nnzmat() / (double)comm_data_->n();
}
// Print results in vtk format for visualization in visIT
void ClusterOrbitals::writeVTKHeader(
    const int npx, const int npy, const int npz, std::ofstream& os)
{
    // write vtk header file
    if (onpe0)
    {
        os << "# vtk DataFile Version 3.1" << endl;
        os << "Subdomain cluster sizes" << endl;
        os << "ASCII" << endl;
        os << "DATASET UNSTRUCTURED_GRID" << endl;
        // define vertex points/ connectivity info
        const int npx1 = npx + 1;
        const int npy1 = npy + 1;
        const int npz1 = npz + 1;

        const int npts = (npx1) * (npy1) * (npz1);
        os << "POINTS " << npts << " FLOAT" << endl;
        for (int iz = 0; iz < npz1; iz++)
        {
            for (int iy = 0; iy < npy1; iy++)
            {
                for (int ix = 0; ix < npx1; ix++)
                {
                    os << ix << " " << iy << " " << iz << endl;
                }
            }
        }
        os << " " << endl;

        // define cells
        const int ncells          = npx * npy * npz;
        const int points_per_cell = 8;
        os << "CELLS " << ncells << " " << ncells * (points_per_cell + 1)
           << endl;
        for (int iz = 0; iz < npz; iz++)
        {
            for (int iy = 0; iy < npy; iy++)
            {
                for (int ix = 0; ix < npx; ix++)
                {
                    int p1 = ix + iy * npx1 + iz * npx1 * npy1;
                    int p5 = p1 + npx1 * npy1;
                    os << points_per_cell << " " << p1 << " " << p1 + 1 << " "
                       << p1 + npx1 << " " << p1 + npx1 + 1 << " " << p5 << " "
                       << p5 + 1 << " " << p5 + npx1 << " " << p5 + npx1 + 1
                       << endl;
                }
            }
        }

        // define cell types
        os << " " << endl;
        os << "CELL_TYPES " << ncells << endl;
        const int cell_type = 11; // use vtks voxel cell type
        for (int cell = 0; cell < ncells; cell++)
        {
            os << cell_type << " ";
        }
        os << endl;
        os << " " << endl;

        // prepare to receive cell data
        os << "CELL_DATA   " << ncells << endl;
    }
}

void ClusterOrbitals::writeVTKDataset(
    const std::string& name, const pb::PEenv& myPEenv, std::ofstream& os)
{
    const int npts = myPEenv.n_mpi_tasks();
    if (onpe0 && os.is_open())
    {
        os << "SCALARS " << name << " FLOAT" << endl;
        os << "LOOKUP_TABLE default" << endl;
    }
    // gather data onto processor 0
    int data        = cluster_indexes_.size();
    int* buff       = (int*)malloc(npts * sizeof(int));
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.gather(&data, 1, &buff[0], npts, 0);

    // write data
    if (onpe0 && os.is_open())
    {
        for (int iz = 0; iz < myPEenv.n_mpi_task(2); iz++)
        {
            for (int iy = 0; iy < myPEenv.n_mpi_task(1); iy++)
            {
                for (int ix = 0; ix < myPEenv.n_mpi_task(0); ix++)
                {
                    os << buff[myPEenv.xyz2task(ix, iy, iz)] << endl;
                }
            }
        }
        os << " " << endl;
    }

    free(buff);
}
