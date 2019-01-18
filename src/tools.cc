// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "tools.h"
#include "Control.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"
#include "Mesh.h"
#include "SparseDistMatrix.h"
#include "Vector3D.h"

#include <string>
#include <sys/stat.h>
#include <time.h>

const int im_    = 2147483647;
const int iq_    = 127773;
const int mask_  = 123459876;
const double am_ = (1.0 / im_);

using namespace std;

// function to call if operator new can't allocate enough memory
void noMoreMemory()
{
    cerr << "Unable to satisfy request for memory for MPI task " << mype
         << endl;
    Control& ct = *(Control::instance());
    ct.global_exit(3);
}

// an atom name should start with a capital letter and end with a number
bool checkValidName(const std::string& name)
{
    if (name.empty()) return false;

    unsigned found = name.find_last_of("0123456789");
    if (name[found] != (*name.rbegin())) return false;

    found = name.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
    if (name[found] != (*name.begin())) return false;

    return true;
}

void stripLeadingAndTrailingBlanks(string& stringToModify)
{
    if (stringToModify.empty()) return;

    int startIndex = stringToModify.find_first_not_of(" ");
    // int endIndex   = stringToModify.find_last_not_of(" ");
    int endIndex = stringToModify.find_last_of("0123456789");
    // cerr<<"startIndex="<<startIndex<<endl;
    // cerr<<"endIndex="<<endIndex<<endl;
    string tempString(stringToModify);
    stringToModify.erase();
    stringToModify = tempString.substr(startIndex, (endIndex - startIndex + 1));
}

void stripName(string& stringToModify)
{
    if (stringToModify.empty()) return;

    // make new string with 1st two characters
    string tempString(stringToModify.substr(0, 2));
    int endIndex = tempString.find_first_of("_0123456789");
    // cout<<"endIndex="<<endIndex<<endl;
    if (endIndex > 0)
        stringToModify = tempString.substr(0, endIndex);
    else
        stringToModify = tempString;
}

void appendNumber(string& name, const int number)
{
    char extension[4];
    sprintf(extension, "%d", number);
    if (number < 10) name.append("0");
    if (number < 100) name.append("0");
    if (number < 1000) name.append("0");
    int l = 1;
    if (number > 9) l++;
    if (number > 99) l++;
    if (number > 999) l++;
    name.append(extension, l);
}

void read_comments(ifstream& tfile)
{
    // while( tfile.get()!='\n');
    // string str;
    // getline(tfile,str);

    char cc = (char)tfile.peek();
    while (cc == ('#') || (cc == '\n') || cc == ' ')
    {
        while (tfile.get() != '\n')
            ;
        cc = (char)tfile.peek(); // look at next character
    }
}

// rotate symmetric matrix mat
void rotateSym(dist_matrix::DistMatrix<DISTMATDTYPE>& mat,
    const dist_matrix::DistMatrix<DISTMATDTYPE>& rotation_matrix,
    dist_matrix::DistMatrix<DISTMATDTYPE>& work)
{
    work.symm('l', 'l', 1., mat, rotation_matrix, 0.);
    mat.gemm('t', 'n', 1., rotation_matrix, work, 0.);
}

void sqrtDistMatrix(dist_matrix::DistMatrix<DISTMATDTYPE>& u)
{
    (*MPIdata::sout) << "sqrtDistMatrix()" << endl;
    const int nst = u.m();
#if 0
    dist_matrix::DistMatrix<DISTMATDTYPE> u0(u);
    dist_matrix::DistMatrix<DISTMATDTYPE> test("t",bc,nst,nst);
#endif
#if 0
    test.gemm('t','n',1.,u0,u,0.);
    for(int i=0;i<test.mloc();i++)
    for(int j=0;j<i;j++){
        if( fabs( test.val(i+j*test.mloc()) )>1.e-8 )
            (*MPIdata::sout)<<"test["<<i<<"]["<<j<<"]="
                <<test.val(i+j*test.mloc())<<endl;
    }
    for(int i=0;i<test.mloc();i++){
        if( fabs( test.val(i+i*test.mloc())-1. )>1.e-8 )
            (*MPIdata::sout)<<"test["<<i<<"]["<<i<<"]="
                <<test.val(i+i*test.mloc())<<endl;
    }

#endif

    dist_matrix::DistMatrix<DISTMATDTYPE> w("w", nst, nst);
    dist_matrix::DistMatrix<DISTMATDTYPE> z("z", nst, nst);

    vector<DISTMATDTYPE> eigenvalues(nst);

    // a = (u+Id).
    dist_matrix::DistMatrix<DISTMATDTYPE> a(u);
    w.identity();
    a.axpy(1., w);

    // w = ( (2*Id + u**T + u) )**(-1/2)
    w.transpose(1., u, 2.);
    w.axpy(1., u);
    w.syev('v', 'l', eigenvalues, z);
    for (int i = 0; i < nst; i++)
    {
        //(*MPIdata::sout)<<"eigenvalues="<<eigenvalues[i]<<endl;
        eigenvalues[i] = 1. / sqrt(eigenvalues[i]);
    }
    // for(int i=0;i<nst;i++)
    //    (*MPIdata::sout)<<"eigenvalues="<<eigenvalues[i]<<endl;
    dist_matrix::DistMatrix<DISTMATDTYPE> g("g", &eigenvalues[0], nst, nst);

    // u = z * g * z**T
    w.symm('r', 'l', 1., g, z, 0.);
    g.gemm('n', 't', 1., w, z, 0.);

    // u = a * g
    u.symm('r', 'l', 1., g, a, 0.);

#if 0
    // verification
    test = u;
    w.gemm('n','n',1.,u,test,0.);
    w.axpy(-1.,u0);
    for(int i=0;i<w.mloc();i++)
    for(int j=0;j<w.nloc();j++){
        if( fabs( w.val(i+j*w.mloc()) )>1.e-8 )
            (*MPIdata::sout)<<"w["<<i<<"]["<<j<<"]="
                <<w.val(i+j*w.mloc())<<endl;
    }
#endif
}

void setSparseDistMatriConsolidationNumber(const int npes)
{
    int consolidation_number = 9;
    int remainder            = (npes % consolidation_number);
    // find consolidation number with max. or 0 remainder;
    if (remainder > 0)
        for (int i = consolidation_number - 2; i < consolidation_number + 8;
             i++)
        {
            if ((npes % i) == 0)
            {
                consolidation_number = i;
                break;
            }
            int new_remainder = npes % i;
            if ((new_remainder) > remainder)
            {
                consolidation_number = i;
                remainder            = new_remainder;
            }
        }
    dist_matrix::SparseDistMatrix<DISTMATDTYPE>::setConsolidationNumber(
        consolidation_number);

    if (onpe0)
        dist_matrix::SparseDistMatrix<DISTMATDTYPE>::printConsolidationNumber(
            *MPIdata::sout);
}

void reduceBytes(vector<char>& val, const MPI_Comm comm)
{
    assert(sizeof(char) == 1);

    const int nn = (int)val.size();

    // pack data into bits
    int nc = nn / 4;
    if (nc * 4 < nn) nc++;

    vector<char> bb(nc, 0);
    for (int i = 0; i < nn; i++)
    {
        if (val[i] == 1)
        {
            int ic = i >> 2;
            int ib = i - 4 * ic;
            switch (ib)
            {
                case 0:
                    bb[ic] = bb[ic] | 1; // set bit 0
                    break;
                case 1:
                    bb[ic] = bb[ic] | 2; // set bit 1
                    break;
                case 2:
                    bb[ic] = bb[ic] | 4; // set bit 2
                    break;
                case 3:
                    bb[ic] = bb[ic] | 8; // set bit 3
                    break;
                default:
                    break;
            }
        }
    }

    // transfer data
    vector<char> bbsum(nc, 0);
    MPI_Allreduce(&bb[0], &bbsum[0], nc, MPI_CHAR, MPI_BOR, comm);
    // cout<<"Allreduce terminated..."<<endl;

    // unpack data
    int i   = 0;
    int nn4 = nn / 4;
    for (int ic = 0; ic < nn4; ic++)
    {
        val[i++] = (char)((bbsum[ic] & 1) > 0);
        val[i++] = (char)((bbsum[ic] & 2) > 0);
        val[i++] = (char)((bbsum[ic] & 4) > 0);
        val[i++] = (char)((bbsum[ic] & 8) > 0);
    }
    for (int ic = nn4; ic < nc; ic++)
    {
        if (i < nn) val[i++] = (char)((bbsum[ic] & 1) > 0);
        if (i < nn) val[i++] = (char)((bbsum[ic] & 2) > 0);
        if (i < nn) val[i++] = (char)((bbsum[ic] & 4) > 0);
        if (i < nn) val[i++] = (char)((bbsum[ic] & 8) > 0);
    }
}

/////////////////////////////////////////////////////////////////////////////////////
// function add or subtract arrays
void arrayops(const double* const a, const double* const b, const double s,
    const double e, const int dim, double* result)
{
    for (int i = 0; i < dim; i++)
    {
        result[i] = e * (a[i] + s * b[i]);
    }
}

void printWithTimeStamp(const string string2print, ostream& os)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.barrier();
    if (onpe0)
    {
        struct tm* local_tm;
        time_t t = time(NULL);
        local_tm = localtime(&t);

        char* dateString = asctime(local_tm);
        char* c          = index(dateString, '\n');
        if (c != NULL) *c = (char)'\0';
        string date(dateString);

        os << date << ": " << string2print << endl;
        os.flush();
    }
#if 0
    static int s=0;
    int r=0;
    int mpierr=MPI_Reduce(&s, &r, 1, MPI_INT, MPI_SUM, 0, mmpi.commGlobal());
    if( mpierr!=MPI_SUCCESS )
    {
        cerr << " Error in MPI!!!" << endl;
        MPI_Abort(mmpi.commGlobal(),1);
    }
    if( r!=mmpi.size()*s && onpe0 )
    {
        cerr << " Error in barrier: "<<r<<"!="<<(mmpi.size()*s)<<endl;
        MPI_Abort(mmpi.commGlobal(),1);
    }
    s++;
    s=(s%128);
#endif
}

bool isOverlaping(const Vector3D& center, const float radius)
{
    Mesh* mymesh             = Mesh::instance();
    const pb::Grid& mygrid   = mymesh->grid();
    const pb::PEenv& myPEenv = mymesh->peenv();

    const double origin[3]
        = { mygrid.origin(0), mygrid.origin(1), mygrid.origin(2) };
    double div_lattice[3];
    for (short i = 0; i < 3; i++)
        div_lattice[i] = mygrid.ll(i) / (double)(myPEenv.n_mpi_task(i));

    double ll[3];
    for (short i = 0; i < 3; i++)
        ll[i] = (double)myPEenv.my_mpi(i) * div_lattice[i] + origin[i];
    // double  ur[3];
    // for(short i=0;i<3;i++)
    //    ur[i] = ll[i] + div_lattice[i];

    const double h[3] = { mygrid.hgrid(0), mygrid.hgrid(1), mygrid.hgrid(2) };

    // correct div_lattice to include only mesh points, not domain
    div_lattice[0] -= h[0];
    div_lattice[1] -= h[1];
    div_lattice[2] -= h[2];

    assert(radius > 1.);

    // subdomain center
    double cc[3] = { 0.5 * (2. * ll[0] + div_lattice[0]),
        0.5 * (2. * ll[1] + div_lattice[1]),
        0.5 * (2. * ll[2] + div_lattice[2]) };
    double t[3]  = { center[0], center[1], center[2] };
    for (short i = 0; i < 3; i++)
    {
        while (t[i] >= cc[i] + 0.5 * mygrid.ll(i))
            t[i] -= mygrid.ll(i);
        while (t[i] < cc[i] - 0.5 * mygrid.ll(i))
            t[i] += mygrid.ll(i);

        t[i] -= ll[i];
    }

    // find "p" closest point within subdomain (projection)
    double p[3] = { t[0], t[1], t[2] };
    if (p[0] < 0.) p[0] = 0.;
    if (p[1] < 0.) p[1] = 0.;
    if (p[2] < 0.) p[2] = 0.;
    if (p[0] > div_lattice[0]) p[0] = div_lattice[0];
    if (p[1] > div_lattice[1]) p[1] = div_lattice[1];
    if (p[2] > div_lattice[2]) p[2] = div_lattice[2];
    double d2 = (p[0] - t[0]) * (p[0] - t[0]) + (p[1] - t[1]) * (p[1] - t[1])
                + (p[2] - t[2]) * (p[2] - t[2]);
    if (d2 < radius * radius)
    {
        return true;
    }

    return false;
}

void exitWithErrorMessage(const string name)
{
    cerr << "Function " << name
         << " not implemented and should not be called!!!" << endl;
    MPI_Finalize();
    exit(0);
    // MPI_Abort(MPI_COMM_WORLD,0);
}

bool fileExists(const char* file)
{
    struct stat buf;
    return (stat(file, &buf) == 0);
}

//  Random number generator of Park and Miller
//  adapted from "Numerical Recipes in C", 2nd. Edition, p.279
//  Returns a uniform random deviate between 0. and 1.
double ran0()
{
    // Initialize the random number generator
    static int idum_ = 3356;

    idum_ ^= mask_;
    long k = idum_ / iq_;
    idum_  = 16807 * (idum_ - k * iq_) - 2836 * k;
    if (idum_ < 0) idum_ += im_;
    double ans = am_ * (double)(idum_);
    idum_ ^= mask_;
    return ans;
}

double minQuadPolynomial(const double e0, const double e1, const double de0,
    const bool print_flag, ostream& os)
{
    double x0 = -de0 / (2. * (e1 - de0 - e0));
    if (x0 < 0. && e1 < e0) x0 = 1.;
    x0 = min(x0, 1.);

    if (onpe0 && print_flag)
    {
        os << setprecision(12) << fixed << "x0=" << x0 << endl;
        double val = x0 * x0 * (e1 - e0 - de0) + x0 * de0 + e0;
        os << "Predicted E: " << val << endl;
    }
    return x0;
}

double minQuadPolynomialFrom3values(const double e0, const double e1,
    const double ehalf, const bool print_flag, ostream& os)
{
    double b = -e1 + 4 * ehalf - 3. * e0;
    double a = -4. * ehalf + 2. * e1 + 2. * e0;
    assert((0.25 * a + 0.5 * b + e0 - ehalf) < 1.e-8);

    if (a < 0.)
    {
        if (e1 < e0)
            return 1.;
        else
            return -1.; // failed to find minimum
    }

    double x0 = -b / (2. * a);
    x0        = min(x0, 1.);

    if (onpe0 && print_flag)
    {
        os << setprecision(12) << fixed << "x0=" << x0 << endl;
        double val = a * x0 * x0 + b * x0 + e0;
        // os<<"a="<<a<<", b="<<b<<endl;
        os << "Predicted E: " << val << endl;
    }
    return x0;
}
