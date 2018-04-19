// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#include "Vector3D.h"
using namespace std;

double Vector3D::minimage ( const Vector3D& w,
                            const Vector3D& l,
                            const short bc[3] )const
{
    //const Vector3D d=(*this)-w;
    Vector3D d;
    d.sub_and_assign((*this),w);

    Vector3D x(d);
    for(short i=0;i<3;i++){
        if( bc[i]==1 ){
            x[i]=remainder(d.x_[i],l.x_[i]);
        }
    }

    return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
}

Vector3D Vector3D::vminimage( const Vector3D& w, const Vector3D& l,
                              const short bc[3] )const
{
    assert( x_[0]==x_[0] );
    assert( x_[1]==x_[1] );
    assert( x_[2]==x_[2] );
    
    //const Vector3D d=(*this)-w;
    Vector3D d;
    d.sub_and_assign((*this),w);

    Vector3D r;

    for(short i=0;i<3;i++)
    {
        if( bc[i]==1 ){
            r.x_[i]=fmod(d.x_[i],l.x_[i]);

            const double ax2=0.5*l.x_[i];

            if( r.x_[i] > ax2 )r.x_[i]-=l.x_[i];
            else if( r.x_[i] < -ax2 )r.x_[i]+=l.x_[i];
            assert( r.x_[i] <= ax2 );
        }else{
            r.x_[i]=d.x_[i];
        }
    }

    return r;
}

void bcastvv3d(vector<Vector3D>& vv, MPI_Comm comm)
{
#ifdef USE_MPI
    const int n=(int)vv.size();
    double* buffer=new double[3*n];
    for(int j=0;j<n;j++)
    for(short i=0;i<3;i++){
        buffer[3*j+i]=vv[j].x_[i];
    }
    int retbcast=MPI_Bcast(buffer, 3*n, MPI_DOUBLE, 0, comm);
    if( retbcast!=MPI_SUCCESS )
    {
        cerr<<"ERROR!!!! bcastvv3d(), Failure in MPI_Bcast of 'radii_'!!!"<<endl;
        MPI_Abort(comm, 0);
    }
    for(int j=0;j<n;j++)
    for(short i=0;i<3;i++){
        vv[j].x_[i]=buffer[3*j+i];
    }
    delete[] buffer;
#endif
}
