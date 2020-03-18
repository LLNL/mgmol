// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Timer.h"

#include <mpi.h>
#include <string>
#include <vector>
class Vector3D;

class MDfiles
{
public:
    MDfiles();

    void createSnapshotDir(const int mdstep, std::string& md_print_dir);
    void printDataInFiles(std::vector<std::string>& ions_names,
        std::vector<double>& tau, std::vector<double>& forces,
        std::vector<double>& taum, std::vector<Vector3D>& centers,
        std::vector<float>& spreads, std::vector<int>& gids, const int mdstep,
        const double dt);

    static void printTimers(std::ostream& os) { print_data_tm_.print(os); }

private:
    static Timer print_data_tm_;

    // MPI group color
    unsigned color_;

    // MPI task within group
    unsigned key_;

    // MPI communicator for group
    MPI_Comm sub_comm_;

    // number of tasks in group
    int size_comm_;

    void appendTaskNumberToName(std::string& name);
    void gatherStrings(
        std::vector<std::string>& ions_names, std::vector<char>& recvbufs);
    void gatherVector3D(
        std::vector<Vector3D>& tau, std::vector<double>& recvbuf);
};
