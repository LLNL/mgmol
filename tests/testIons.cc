#include "Control.h"
#include "Ions.h"
#include "MGmol_MPI.h"
#include "Mesh.h"
#include "Species.h"

#include <random>

int main(int argc, char** argv)
{
    int mpirc = MPI_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD;

    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    MGmol_MPI::setup(comm, std::cout);
    Control::setup(comm, false, 0.);

    // create a domain [0.10.]^3
    const double origin[3]  = { 0., 0., 0. };
    const double ll         = 10.;
    const double lattice[3] = { ll, ll, ll };
    const unsigned ngpts[3] = { 32, 24, 20 };
    short lap_type          = 0;

    Mesh::setup(comm, ngpts, origin, lattice, lap_type);

    const double h[3] = { ll / (double(ngpts[0])), ll / (double(ngpts[1])),
        ll / (double(ngpts[2])) };

    // random number generator
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dis(0.0, 1.0);

    // create one species
    Species sp(MPI_COMM_WORLD);

    // read species info from pseudopotential file
    std::string file_path = argv[1];
    std::string filename(file_path + "/pseudo.C_ONCV_PBE_SG15");
    std::cout << "Potential = " << filename << std::endl;

    sp.read_1species(filename);
    sp.set_dim_nl(h[0]);
    sp.set_dim_l(h[0]);
    sp.initPotentials('f', h[0], true);

    // put species into a vector
    std::vector<Species> vsp;
    vsp.push_back(sp);

    Ions ions(lattice, vsp);
    ions.setupListIonsBoundaries(10000.);

    double velocity[3] = { 0., 0., 0. };

    // set "na" atoms coordinates and add them to "ions"
    const int na = 10;
    for (int i = 0; i < na; i++)
    {
        double x[3] = { origin[0] + lattice[0] * dis(gen),
            origin[1] + lattice[1] * dis(gen),
            origin[2] + lattice[2] * dis(gen) };
        if (myrank == 0)
            std::cout << "x,y,z = " << x[0] << ", " << x[1] << ", " << x[2]
                      << std::endl;

        // set all x to the values of PE0
        MPI_Bcast(&x[0], 3, MPI_DOUBLE, 0, comm);

        // make a name for atom based on species and order of reading in
        std::string stri = std::to_string(i);
        std::string aname("C" + stri);

        ions.addIonToList(sp, aname, &x[0], velocity, false);
    }

    ions.setup();

    std::vector<Ion*>& new_local_ions(ions.local_ions());

    int nlocal = new_local_ions.size();
    std::cout << "PE " << myrank << ", nlocal = " << nlocal << std::endl;

    int ntotal = 0;
    MPI_Allreduce(&nlocal, &ntotal, 1, MPI_INT, MPI_SUM, comm);
    mpirc = MPI_Finalize();
    if (mpirc != MPI_SUCCESS)
    {
        std::cerr << "MPI Finalize failed!!!" << std::endl;
        return 1;
    }

    if (ntotal != na)
    {
        std::cout << "ntotal = " << ntotal << std::endl;
        return 1;
    }

    return 0;
}
