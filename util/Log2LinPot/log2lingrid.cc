// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

void spline(double*, double*, const int, double, double, double*);
void splint(double*, double*, double*, const int, double, double*);

int main(int argc, char** argv)
{
    string input_file = argv[1];
    ifstream* tfile   = new ifstream(input_file.data(), ios::in);

    string str;
    int count = 0;

    int npot = 0;
    int npts = -1;
    while (count < 14)
    {
        getline(*tfile, str);
        istringstream ss(str);
        if (str[0] != '#')
        {
            count++;
            if (count == 9) ss >> npot;
            if (count == 13) ss >> npts;
            if (count == 14) str = '0';
        }
        std::cout << str << std::endl;
    }

    std::cout << "#npts : " << npts << std::endl;
    std::cout << "#npot : " << npot << std::endl;

    vector<double> x(npts);
    vector<double> psi(npts);
    vector<double> pot(npts);

    for (int l = 0; l < npot; l++)
    {
        for (int i = 0; i < npts; i++)
        {
            do
            {
                getline(*tfile, str);
                // cout<<str<<endl;
            } while (str[0] == '#');
            istringstream ss(str);
            ss >> x[i] >> pot[i] >> psi[i];
            // psi[i]/=x[i];
        }
        x[0]       = 0.;
        double ypr = -1. / (x[npts - 1] * x[npts - 1]);
        vector<double> ypp_pot(npts);
        vector<double> ypp_psi(npts);
        double zero  = 0.;
        double large = 1.e32;
        spline(&x[0], &pot[0], npts, zero, ypr, &ypp_pot[0]);
        spline(&x[0], &psi[0], npts, large, zero, &ypp_psi[0]);

        double yp_val, ypp_val;
        std::cout << "#l=" << l << std::endl;
        for (int i = 0; i < npts; i++)
        {
            double r = i * 0.01;
            double vpot;
            splint(&x[0], &pot[0], &ypp_pot[0], npts, r, &vpot);
            double vpsi;
            splint(&x[0], &psi[0], &ypp_psi[0], npts, r, &vpsi);
            std::cout << setw(5) << setprecision(8) << r << "\t" << setw(10)
                      << vpot << "\t" << vpsi << std::endl;
        }
    }
}

