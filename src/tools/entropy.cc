// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "entropy.h"
#include "fermi.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include "MPIdata.h"

using namespace std;

/* Evaluate the entropy function from occupations, f: 
   s(f) = -f ln(f) - (1-f) ln(1-f).
   Returns the entropy = trace[s(f)] = sum {s(f)}_i
*/
double entropy_eval(const std::vector<double>& f, std::vector<double>& s, const double occ_factor)
{
    const int dim = f.size();
    assert( s.size() == dim );

    const double tol=1.e-15;

#ifndef NDEBUG
    const double tol_interval = 1.e-6;
#endif

    double entropy = 0.;

    for (int st = 0; st < dim; st++)
    {
        const double fi = (double)f[st];
//        if( fi>1.+tol_interval )
//            (*MPIdata::sout)<<setprecision(15)<<scientific<<"f["<<st<<"]="<<fi<<endl;
        assert( fi>=0.-tol_interval );
        assert( fi<=1.+tol_interval );
        if(fi<tol){
            s[st] = (1.-fi)*log(1.-fi);
        }else if( fi>1.-tol ){
            s[st] = fi*log(fi);
        }else{
            s[st] = fi*log(fi)+(1.-fi)*log(1.-fi);
        }
        entropy += s[st];
        
//        if(onpe0)std::cout<<"("<<f[st]<<", "<<s[st]<<")"<<std::endl;
    }
//    if(onpe0)std::cout<<"entropy_val = "<<entropy<<endl;
    return (double)(-occ_factor)*entropy; // in units of kbt
}

double entropy_eval(const std::vector<float>& f, std::vector<float>& s, const double occ_factor)
{
    const int dim = f.size();
    assert( s.size() == dim );

    const double tol=1.e-15;

#ifndef NDEBUG
    const double tol_interval = 1.e-6;
#endif
    
    double entropy = 0.;

    for (int st = 0; st < dim; st++)
    {
        const double fi = (double)f[st];
//        if( fi>1.+tol_interval )
//            (*MPIdata::sout)<<setprecision(15)<<scientific<<"f["<<st<<"]="<<fi<<endl;
        assert( fi>=0.-tol_interval );
        assert( fi<=1.+tol_interval );
        if(fi<tol){
            s[st] = (float)((1.-fi)*log(1.-fi));
        }else if( fi>1.-tol ){
            s[st] = (float)(fi*log(fi));
        }else{
            s[st] = (float)(fi*log(fi)+(1.-fi)*log(1.-fi));
        }
        entropy += s[st];
    }
    
    return (double)(-occ_factor)*entropy; // in units of kbt
}

// Evaluate entropy given 'energies'
double entropy_evalFromEnergies(const double mu,
                          const int max_occ,
                          const double kBT,
                          const std::vector<double>& energies, 
                          std::vector<double>& s, 
                          const double occ_factor)
{
   std::vector<double> f((int)energies.size(),0.);
   
   // calculate occupations based on fermi-dirac function
   fermi_distribution(mu,max_occ,kBT,energies,f );
   
   // call entropy function
   double ent = entropy_eval(f,s,occ_factor);   
   
   return ent;
}

// Evaluate entropy given 'energies'
double entropy_evalFromEnergies(const double mu,
                          const int max_occ,
                          const double kBT,
                          const std::vector<float>& energies, 
                          std::vector<float>& s, 
                          const double occ_factor)
{
   std::vector<float> f((int)energies.size(),0.);
   
   // calculate occupations based on fermi-dirac function
   fermi_distribution(mu,max_occ,kBT,energies,f );
   
   // call entropy function
   double ent = entropy_eval(f,s,occ_factor);   
   
   return ent;
}
