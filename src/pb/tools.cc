// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// Some general tools
// $Id: tools.cc,v 1.6 2008/12/10 01:05:37 jeanluc Exp $
#include <iostream>
#include <ctime>
#include <sys/time.h>
using namespace std;
#include <stdlib.h>

#include "tools.h"

namespace pb{

FILE* open_file(const string filename, const string mode)
{
    FILE*   file_ptr;
    
    if( (file_ptr = fopen(filename.data(), mode.data())) == NULL){
        cout<<"\n cannot open "<<filename<<endl;
        exit(2);
    }else{
        cout<<" open "<<filename<<endl;
    }
    
    return file_ptr;
}

double timer(void)
{
  struct timeval  tt;

  gettimeofday( &tt, NULL );
  double val1 = (double)tt.tv_usec;
  val1 /= 1000000.;
  double val = (double)tt.tv_sec;
  val += val1;
  return val;
}

} // namespace pb

