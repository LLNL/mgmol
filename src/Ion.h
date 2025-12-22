// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_ION_H
#define MGMOL_ION_H

#include "GridFunc.h"
#include "IonData.h"
#include "KBprojector.h"
#include "Mesh.h"
#include "Species.h"

#include <cstring>
#include <memory>
#include <vector>

// Ion structure
class Ion
{
    // Name of Ion
    const std::string name_;

    const Species& species_;

    // unique index
    const unsigned int index_;

    // unique start index for nl projectors
    const unsigned int nlproj_gid_;

    // actual coordinates
    double position_[3];
    double old_position_[3];

    // atom locked (1) or not (0)
    char locked_;

    // tells if the nl projectors of this ion are non-zero for this task
    bool map_nl_;

    // tells if the local potential of this ion is non-zero for this task
    bool map_l_;

    // here=true if global contribution has to be computed on this process
    bool here_;

    std::shared_ptr<KBprojector> kbproj_;

    // indices of first grid point where the local pot. is non-zero
    short lpot_start_[3];

    // Coordinates of the corner of the grid that the local
    // difference potential is nonzero on.
    double lstart_[3];

    // Forces on the ion
    double force_[3];

    // Velocity of the ion
    double velocity_[3];

    // random state for thermostat
    unsigned short rand_state_[3];

    void shift_position(const short i, const double shift)
    {
        old_position_[i] = position_[i];

        position_[i] += shift;
    }

public:
    Ion(const Species& species, const std::string& name, const double crds[3],
        const double velocity[3], const bool lock = false);

    Ion(const Species& species, const std::string& name, const double crds[3],
        const double velocity[3], const unsigned int index,
        const unsigned int nlproj_id, const bool lock = false);
    Ion(const Species& species, IonData data);
    Ion(const Ion& ion);

    ~Ion() {};

    void init(const double crds[3], const double velocity[3], const bool lock);
    void setup();

    void setPreviousPosition(const double x, const double y, const double z)
    {
        old_position_[0] = x;
        old_position_[1] = y;
        old_position_[2] = z;
    }

    std::shared_ptr<KBprojector> kbproj() { return kbproj_; }
    const std::shared_ptr<KBprojector> kbproj() const { return kbproj_; }

    double radiusNLproj() const { return kbproj_->maxRadius(); }

    double computeRadiusVl() const
    {
        Mesh* mymesh           = Mesh::instance();
        const pb::Grid& mygrid = mymesh->grid();

        double radius = 0.;

        for (short dir = 0; dir < 3; dir++)
        {
            const double hgrid = mygrid.hgrid(dir);

            const double left = position_[dir] - lstart_[dir];
            const double right
                = lstart_[dir] + species_.dim_l() * hgrid - position_[dir];

            radius = left > radius ? left : radius;
            radius = right > radius ? right : radius;
        }

        return radius;
    }

    const Species& getSpecies() const { return species_; }

    double getRadiusLocalPot() const { return species_.lradius(); }

    double getRC() const { return species_.rc(); }
    bool isMass28() const { return species_.isMass28(); }
    bool isMassLargerThan1() const { return species_.isMassLargerThan1(); }
    const RadialInter& getLocalPot() const { return species_.local_pot(); }

    void get_Ai(std::vector<int>& Ai, const int gdim, const short dir) const;

    bool operator<(const Ion& ion) const
    {
        // return (10000.*position_[0]+100.*position_[1]+position_[2]
        //      <
        //      10000.*ion.position_[0]+100.*ion.position_[1]+ion.position_[2]);
        return (index_ < ion.index_);
    }

    std::string name() const { return name_; }
    bool compareName(const std::string& name) const
    {
        return (name_.compare(name) == 0);
    }
    unsigned int index() const { return index_; }
    bool compareIndex(const unsigned int index) const
    {
        return (index == index_);
    }
    unsigned int nlprojid() const { return nlproj_gid_; }

    int atomic_number() const { return species_.getAtomicNumber(); }

    bool map_nl() const { return map_nl_; }
    bool map_l() const { return map_l_; }
    void set_map_l(bool val) { map_l_ = val; }
    short lpot_start(const unsigned short i) const { return lpot_start_[i]; }
    double lstart(const unsigned short i) const { return lstart_[i]; }

    char locked(void) const { return locked_; }
    void lock(void) { locked_ = 1; }
    bool here(void) const { return here_; }
    void set_here(bool val) { here_ = val; }

    double eself() const { return species_.eself(); }

    double getZion() const
    {
        assert(species_.rc() > 1.e-10);
        return (double)species_.zion();
    }

    double position(const short i) const { return position_[i]; }
    double getPreviousPosition(const short i) const { return old_position_[i]; }
    void setPosition(const double x, const double y, const double z)
    {
        old_position_[0] = position_[0];
        old_position_[1] = position_[1];
        old_position_[2] = position_[2];

        position_[0] = x;
        position_[1] = y;
        position_[2] = z;

        kbproj_->clear();
    }

    void resetPositionsToPrevious();

    void shiftPositionXLBOMDTest(Vector3D shift)
    {
        for (short dir = 0; dir < 3; dir++)
            shift_position(dir, shift[dir]);

        kbproj_->clear();
    }
    void shiftPosition(const double shift[3])
    {
        for (short dir = 0; dir < 3; dir++)
            shift_position(dir, shift[dir]);
    }
    double norm2F() const
    {
        double ff = force_[0] * force_[0] + force_[1] * force_[1]
                    + force_[2] * force_[2];
        return ff;
    }

    short nProjectors() const { return kbproj_->nProjectors(); }

    short nProjectorsSubdomain() const
    {
        return kbproj_->nProjectorsSubdomain();
    }

    double velocity(const short i) const { return velocity_[i]; }
    void setVelocity(const double v0, const double v1, const double v2)
    {
        velocity_[0] = v0;
        velocity_[1] = v1;
        velocity_[2] = v2;
    }
    void rescaleVelocity(const double factor)
    {
        velocity_[0] *= factor;
        velocity_[1] *= factor;
        velocity_[2] *= factor;
    }
    double force(const short i) const
    {
        assert(i < 3);
        return force_[i];
    }
    void set_force(const short i, const double val) { force_[i] = val; }
    void add_force(const double f0, const double f1, const double f2)
    {
        force_[0] += f0;
        force_[1] += f1;
        force_[2] += f2;
    }
    void getPosition(double* const position) const
    {
        position[0] = position_[0];
        position[1] = position_[1];
        position[2] = position_[2];
    }
    void getForce(double* const force) const
    {
        force[0] = force_[0];
        force[1] = force_[1];
        force[2] = force_[2];
    }
    void setForce(const double v0, const double v1, const double v2)
    {
        force_[0] = v0;
        force_[1] = v1;
        force_[2] = v2;
    }

    double distance(const Ion& other_ion) const
    {
        double d2 = 0.;
        for (short i = 0; i < 3; i++)
        {
            double d = (position_[i] - other_ion.position_[i]);
            d2 += d * d;
        }
        return sqrt(d2);
    }
    void printPosition(std::ostream& os) const;
    void printPositionAndForce(std::ostream& os) const;

    void resetForce()
    {
        for (short i = 0; i < 3; i++)
        {
            force_[i] = 0.;
        }
    }
    void bcast(MPI_Comm);

    double minimage(const Ion&, const double cell[3], const short periodic[3],
        double xc[3]) const;
    double distance(const Ion&, double*, double*, double*) const;
    double minimage(
        const Ion&, const double cell[3], const short periodic[3]) const;
    double distance(const Ion&, const double, const double, const double) const;
    double minimage(
        const double point[3], const double cell[3], const short bc[3]) const;

    void set_lstart(const int, const double, const double);

    double getMass() const { return species_.getMass(); }
    unsigned short randomState(const short i) const { return rand_state_[i]; }
    void setRandomState(const unsigned short s0, const unsigned short s1,
        const unsigned short s2)
    {
        rand_state_[0] = s0;
        rand_state_[1] = s1;
        rand_state_[2] = s2;
    }
    void setFromIonData(const IonData& data);

    void getIonData(IonData& data) const;

    void getGidsNLprojs(std::vector<int>& gids) const;
    void getKBsigns(std::vector<short>& kbsigns) const;
    void getKBcoeffs(std::vector<double>& coeffs);
    double energyDiff(
        Ion& ion, const double lattice[3], const short bc[3]) const;

    static void resetIndexCount();
};

#endif
