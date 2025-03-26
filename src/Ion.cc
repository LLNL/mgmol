// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Ion.h"
#include "KBprojectorSparse.h"
#include "MGmol_blas1.h"

#include <iomanip>
#include <iostream>

#include <mpi.h>

unsigned short isqrt(unsigned value)
{
    return static_cast<unsigned short>(sqrt(static_cast<double>(value)));
}

static unsigned int _nlproj_gid = 0;
static unsigned int _index      = 0;

void Ion::resetIndexCount() { _index = 0; }

Ion::Ion(const Species& species, const std::string& name, const double crds[3],
    const double velocity[3], const bool lock)
    : name_(name), species_(species), index_(_index), nlproj_gid_(_nlproj_gid)
{
    assert(name.size() > 0);

    kbproj_.reset(new KBprojectorSparse(species));

    _index++;
    _nlproj_gid += nProjectors();

    init(crds, velocity, lock);
}

Ion::Ion(const Species& species, const std::string& name, const double crds[3],
    const double velocity[3], const unsigned int index,
    const unsigned int nlproj_gid, const bool lock)
    : name_(name), species_(species), index_(index), nlproj_gid_(nlproj_gid)
{
    assert(name.size() > 0);

    kbproj_.reset(new KBprojectorSparse(species));

    init(crds, velocity, lock);
}

Ion::Ion(const Species& species, IonData data)
    : name_(data.ion_name),
      species_(species),
      index_(data.index),
      nlproj_gid_(data.nlproj_id)
{
    kbproj_.reset(new KBprojectorSparse(species));

    const bool lock = data.atmove ? false : true;
    init(data.current_position, data.velocity, lock);
}

Ion::Ion(const Ion& ion)
    : name_(ion.name_),
      species_(ion.species_),
      index_(ion.index_),
      nlproj_gid_(ion.nlproj_gid_)
{
    kbproj_.reset(new KBprojectorSparse(species_));

    for (short i = 0; i < 3; i++)
    {
        position_[i]     = ion.position_[i];
        old_position_[i] = ion.old_position_[i];
        lpot_start_[i]   = ion.lpot_start_[i];
        lstart_[i]       = ion.lstart_[i];
        force_[i]        = ion.force_[i];
        velocity_[i]     = ion.velocity_[i];
    }
    locked_ = ion.locked_;
    map_nl_ = ion.map_nl_;
    map_l_  = ion.map_l_;
    here_   = ion.here_;
}

// initialize data members
void Ion::init(const double crds[3], const double velocity[3], const bool lock)
{
    for (short i = 0; i < 3; i++)
    {
        position_[i]     = crds[i];
        old_position_[i] = crds[i];
        lpot_start_[i]   = 0;
        lstart_[i]       = 0.;
        force_[i]        = 0.;
        velocity_[i]     = velocity[i];
    }
    locked_ = lock;
#if DEBUG
    (*MPIdata::sout) << " position:" << position_[0] << "," << position_[1]
                     << "," << position_[2] << std::endl;
    if (locked_) (*MPIdata::sout) << name_ << " locked" << std::endl;
#endif
    here_   = false;
    map_nl_ = false;
    map_l_  = false;

    // initialize random state seed
    int idx             = index_ + 1; // shift by 1 to avoid divide by zero
    unsigned short root = isqrt(idx);
    unsigned short rem  = (short)(idx % root);
    rand_state_[0]      = root;
    rand_state_[1]      = root;
    rand_state_[2]      = rem;
}

void Ion::setup()
{
    kbproj_->setup(position_);

    map_nl_ = kbproj_->overlapPE();
}

void Ion::bcast(MPI_Comm comm)
{
    MPI_Bcast(old_position_, 3, MPI_DOUBLE, 0, comm);
    MPI_Bcast(position_, 3, MPI_DOUBLE, 0, comm);
    MPI_Bcast(force_, 3, MPI_DOUBLE, 0, comm);
    MPI_Bcast(velocity_, 3, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&locked_, 1, MPI_CHAR, 0, comm);
}

void Ion::set_lstart(const int i, const double cell_origin, const double hgrid)
{
    // t1 = nb of nodes between the boundary and crds
    double t1 = (position_[i] - cell_origin) / hgrid;

    // get the integral part of t1 in ic
    double t2;
    t1     = modf(t1, &t2);
    int ic = (int)t2;
    if (t1 > 0.5) ic++;
    if (t1 < -0.5) ic--;

    lpot_start_[i] = ic - (species_.dim_l() >> 1);
    lstart_[i]     = cell_origin + hgrid * lpot_start_[i];
    //(*MPIdata::sout)<<"lpot_start_[i]="<<lpot_start_[i]<<std::endl;
    //(*MPIdata::sout)<<"lstart_[i]    ="<<lstart_[i]<<std::endl;
    //(*MPIdata::sout)<<"species_.dim_l()="<<species_.dim_l()<<std::endl;
}

double Ion::minimage(const Ion& ion2, const double cell[3],
    const short periodic[3], double xc[3]) const
{
    for (short i = 0; i < 3; i++)
    {
        assert(cell[i] > 0.);
        xc[i] = position_[i] - ion2.position(i);

        if (periodic[i])
        {
            const double half_lattice = 0.5 * cell[i];
            xc[i]                     = fmod(xc[i], cell[i]);
            if (xc[i] > half_lattice)
                xc[i] -= cell[i];
            else if (xc[i] < -half_lattice)
                xc[i] += cell[i];
            assert(xc[i] <= half_lattice);
        }
    }

    const double rmin = sqrt(xc[0] * xc[0] + xc[1] * xc[1] + xc[2] * xc[2]);
    return rmin;
}

double Ion::distance(const Ion& ion2, double* xc, double* yc, double* zc) const
{
    *xc = position_[0] - ion2.position(0);
    *yc = position_[1] - ion2.position(1);
    *zc = position_[2] - ion2.position(2);

    const double r = (*xc) * (*xc) + (*yc) * (*yc) + (*zc) * (*zc);

    return sqrt(r);
}

double Ion::minimage(
    const Ion& ion2, const double cell[3], const short periodic[3]) const
{
    double xc[3] = { 0., 0., 0. };

    return minimage(ion2, cell, periodic, xc);
}

double Ion::distance(
    const Ion& ion2, const double ax, const double ay, const double az) const
{
    assert(ax > 0.);
    assert(ay > 0.);
    assert(az > 0.);

    const double dx = position_[0] - ion2.position(0);
    const double dy = position_[1] - ion2.position(1);
    const double dz = position_[2] - ion2.position(2);
    const double x  = fmod(dx, ax);
    const double y  = fmod(dy, ay);
    const double z  = fmod(dz, az);
    const double r  = x * x + y * y + z * z;

    return sqrt(r);
}

double Ion::minimage(
    const double point[3], const double cell[3], const short periodic[3]) const
{
    double xc[3];
    for (short i = 0; i < 3; i++)
    {
        assert(cell[i] > 0.);
        xc[i] = position_[i] - point[i];

        if (periodic[i] == 1)
        {
            const double half_lattice = 0.5 * cell[i];
            xc[i]                     = fmod(xc[i], cell[i]);
            if (xc[i] > half_lattice)
                xc[i] -= cell[i];
            else if (xc[i] < -half_lattice)
                xc[i] += cell[i];
            assert(xc[i] <= half_lattice);
        }
    }

    const double rmin = sqrt(xc[0] * xc[0] + xc[1] * xc[1] + xc[2] * xc[2]);
    return rmin;
}

bool operator<(const Ion& A, const Ion& B)
{
    return (A.position(0) < B.position(0));
}

bool operator==(Ion A, Ion B) { return (&A == &B); }

void Ion::printPosition(std::ostream& os) const
{
    os << "$$ ";
    if (locked_)
        os << "*";
    else
        os << " ";
    os << std::setw(4) << name_ << std::setw(10) << std::setprecision(4)
       << std::fixed << position_[0] << std::setw(10) << position_[1]
       << std::setw(10) << position_[2] << std::endl;
}

void Ion::printPositionAndForce(std::ostream& os) const
{
    os << "## ";
    if (locked_)
        os << "*";
    else
        os << " ";
    os << std::setw(4) << name_ << setiosflags(std::ios::right) << std::setw(10)
       << std::setprecision(4) << std::fixed << position_[0] << std::setw(10)
       << position_[1] << std::setw(10) << position_[2] << std::setprecision(7)
       << std::scientific << std::setw(16) << force_[0] << std::setw(16)
       << force_[1] << std::setw(16) << force_[2] << std::endl;
}

void Ion::getIonData(IonData& idata) const
{
    idata.ion_name   = name_;
    idata.atomic_num = atomic_number();
    idata.index      = index_;
    idata.nlproj_id  = nlproj_gid_;
    idata.pmass      = getMass();

    if (locked_)
        idata.atmove = 0;
    else
        idata.atmove = 1;

    for (int pos = 0; pos < 3; pos++)
    {
        idata.rand_state[pos]       = rand_state_[pos];
        idata.initial_position[pos] = 0.;
        idata.old_position[pos]     = old_position_[pos];
        idata.current_position[pos] = position_[pos];
        idata.velocity[pos]         = velocity_[pos];
        idata.force[pos]            = force_[pos];
    }
}

void Ion::resetPositionsToPrevious()
{
    for (int pos = 0; pos < 3; pos++)
        position_[pos] = old_position_[pos];
}

void Ion::setFromIonData(const IonData& data)
{
    // random state
    setRandomState(data.rand_state[0], data.rand_state[1], data.rand_state[2]);
    // previous position
    setPreviousPosition(
        data.old_position[0], data.old_position[1], data.old_position[2]);
    // velocity
    setVelocity(data.velocity[0], data.velocity[1], data.velocity[2]);
    // force
    setForce(data.force[0], data.force[1], data.force[2]);
}

// get gids (row index for <KB,Psi> matrix) for all the projectors
// of Ion
void Ion::getGidsNLprojs(std::vector<int>& gids) const
{
    // get global indexes by adding nlproj_gid_
    gids.clear();

    short nproj = nProjectors();
    for (short i = 0; i < nproj; ++i)
    {
        const int gid = nlproj_gid_ + i;
        gids.push_back(gid);
    }
}

// get sign factor for <KB,Psi> matrix for all the projectors
// of Ion
void Ion::getKBsigns(std::vector<short>& kbsigns) const
{
    kbproj_->getKBsigns(kbsigns);
}

void Ion::getKBcoeffs(std::vector<double>& coeffs)
{
    kbproj_->getKBcoeffs(coeffs);
}

void Ion::get_Ai(std::vector<int>& Ai, const int gdim, const short dir) const
{
    int idx = lpot_start(dir);

    int diml = species_.dim_l();
    diml     = std::min(diml, gdim - 1);
    Ai.resize(diml);

    const int sizeAi = (int)Ai.size();
    for (int ix = 0; ix < sizeAi; ix++)
    {
        Ai[ix] = idx % gdim;

        while (Ai[ix] < 0)
        {
            Ai[ix] += gdim; // periodic BC
        }
        assert(Ai[ix] >= 0);
        assert(Ai[ix] < gdim);

        idx++;
    }
}

double Ion::energyDiff(
    Ion& ion, const double lattice[3], const short bc[3]) const
{
    const double r = minimage(ion, lattice, bc);
    assert(r > 0.);

    return species_.ediff(ion.getSpecies(), r);
}
