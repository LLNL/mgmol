/******************************************************************************
$Id:$

 Yana:
   set hdf_dir = /usr/local/tools/hdf5-intel-serial-1.8.5
   icpc -I$hdf_dir/include -g -c -o read_hdf5.o read_hdf5.cc
   icpc -o read_hdf5 read_hdf5.o $hdf_dir/lib/libhdf5.a -openmp -lz
 IBM:
   set hdf_dir =  /usr/local/tools/hdf5/hdf5-1.6.5/serial
   mpxlC -DMPICH_IGNORE_CXX_SEEK -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -I$hdf_dir/include -g -c -o read_hdf5.o read_hdf5.cc
   mpxlC read_hdf5.o -L/bgl/local/lib  $hdf_dir/lib/libhdf5.a -openmp -lz -o read_hdf5
   
*******************************************************************************/

/*
  Usage example:
          read_hdf5 [-bov] file.hdf5 Function0003
*/
#include <iostream>
#include <fstream>
#include <malloc.h>
#include <vector>
#include <cstring>
#include <map>
#include <set>
#include <hdf5.h>

using namespace std;

map<int,double> cov_radii;
map<int,double> ball_radii;
map<int,string> colors;

const double bohr2ang = 0.529177;  // Bohr to Angstroem

const unsigned int MaxStrLength = 7;

struct FixedLengthString {
    char mystring[MaxStrLength];
};


void set_maps()
{
    // Hydrogen
    cov_radii.insert( pair<int,double>(1,0.6) );
    ball_radii.insert( pair<int,double>(1,0.4) );
    colors.insert( pair<int,string>(1,"White") );

    // Lithium
    cov_radii.insert( pair<int,double>(3,2.7) );
    ball_radii.insert( pair<int,double>(3,1.6) );
    colors.insert( pair<int,string>(3,"Ochre") );

     // Berylium
    cov_radii.insert( pair<int,double>(4,2.0) );
    ball_radii.insert( pair<int,double>(4,1.1) );
    colors.insert( pair<int,string>(4,"Ochre") );

   // Carbon
    cov_radii.insert( pair<int,double>(6,1.8) );
    ball_radii.insert( pair<int,double>(6,0.7) );
    colors.insert( pair<int,string>(6,"Cyan") );

    // Nitrogen
    cov_radii.insert( pair<int,double>(7,1.8) );
    ball_radii.insert( pair<int,double>(7,0.6) );
    colors.insert( pair<int,string>(7,"Blue") );

    // Oxygen
    cov_radii.insert( pair<int,double>(8,1.8) );
    ball_radii.insert( pair<int,double>(8,0.5) );
    colors.insert( pair<int,string>(8,"Red") );

    // F
    cov_radii.insert( pair<int,double>(9,1.8) );
    ball_radii.insert( pair<int,double>(9,0.4) );
    colors.insert( pair<int,string>(9,"Green") );

     // Na
    cov_radii.insert( pair<int,double>(11,3.2) );
    ball_radii.insert( pair<int,double>(11,1.9) );
    colors.insert( pair<int,string>(11,"Ochre") );

     // Mg
    cov_radii.insert( pair<int,double>(12,2.9) );
    ball_radii.insert( pair<int,double>(12,1.5) );
    colors.insert( pair<int,string>(12,"Ochre") );

    // Al
    cov_radii.insert( pair<int,double>(13,2.5) );
    ball_radii.insert( pair<int,double>(13,1.2) );
    colors.insert( pair<int,string>(13,"Ochre") );

    // Si
    cov_radii.insert( pair<int,double>(14,2.3) );
    ball_radii.insert( pair<int,double>(14,1.1) );
    colors.insert( pair<int,string>(14,"Ochre") );

    // P
    cov_radii.insert( pair<int,double>(15,2.3) );
    ball_radii.insert( pair<int,double>(15,1.1) );
    colors.insert( pair<int,string>(15,"Tan") );

    // S
    cov_radii.insert( pair<int,double>(16,2.2) );
    ball_radii.insert( pair<int,double>(16,1.1) );
    colors.insert( pair<int,string>(16,"Yellow") );

    // Cl
    cov_radii.insert( pair<int,double>(17,2.2) );
    ball_radii.insert( pair<int,double>(17,1.1) );
    colors.insert( pair<int,string>(17,"Ochre") );
}

/*
  Read dataset in hdf5 file. Returns a pointer to an array containing the data.
  Input:
          filename:    name of hdf5 file
          datasetname: name of dataset
  
  Output:
          dims:        dimensions of data (3D)
          lattice:     lattice parameters
          origin:      cell origin
          attributes:  address of array for attributes (orbital centers+radii)
          dima:        number of attributes
*/
double* get_function(char* filename, char* datasetname, 
                     hsize_t* dims,
                     double* origin, double* lattice, 
                     double** attributes, hsize_t* dima)
{   
    double* data;
    
    herr_t  status;
    htri_t  ishdf;
    hid_t   file_id, dset_id, filespace, attribute_id, attdataspace;

    hsize_t maxdims[3];
    hsize_t dimsa[2];
    int     rank;
    char    attname[50];
 
    /* check if file is in HDF5 format */
    ishdf=H5Fis_hdf5(filename);

    if( ishdf< 0){
        cerr<<"H5Fis_hdf5 unsuccessful"<<endl;
        return NULL;
    }else if( !ishdf ){
        cerr<<"Input file "<<filename<<" not in HDF5 format. Stop."<<endl;
        return NULL;
    }else{

        file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
        if(file_id<0){
            cerr<<"H5Fopen failed for "<<filename<<endl;
            return NULL;
        }
   
        /* Open dataset. */
        dset_id = H5Dopen2(file_id, datasetname,H5P_DEFAULT);
        if( dset_id<0 ){
            cerr<<"H5Dopen failed for datasetname "<<datasetname<<endl;
            return NULL;
        }

        filespace = H5Dget_space(dset_id);
        if( filespace<0 ){
            cerr<<"H5Dget_space failed."<<endl;
            return NULL;
        }

        /* Get dataspace size */
        rank=H5Sget_simple_extent_dims(filespace, dims, maxdims);
        H5Sclose(filespace);
        if( rank!=3 ){
            cerr<<"Problem with dataspace dimension, rank="<<rank<<endl;
            return NULL;
        }else{
            clog<<"Dataspace: dimension "<<(int)dims[0]<<" x "
                                         <<(int)dims[1]<<" x "
                                         <<(int)dims[2]<<endl;
            clog<<"Size: "<<(int)(dims[0]*dims[1]*dims[2])<<endl;
        }
        if( (int)(dims[0]*dims[1]*dims[2])< 1 )return NULL;

        /* Read data -> data */
        data=new double[dims[0]*dims[1]*dims[2]];
        status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, data);
        if( status<0 ){
            cerr<<"H5Dread failed."<<endl;
            return NULL;
        }

        /* Read lattice parameters */
        sprintf(attname,"%s","Lattice parameters");
        /*  Open a dataset attribute. */
        attribute_id = H5Aopen_name(dset_id, attname);
        if( attribute_id<0 ){
            cerr<<"H5Aopen failed for "<<attname<<endl;
        }else{
            status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, lattice);
            if( status<0 ){
                cerr<<"H5Aread failed."<<endl;
                return NULL;
            }
            status = H5Aclose(attribute_id);
            if( status<0 ){
                cerr<<"H5Aclose failed."<<endl;
                return NULL;
            }
            clog<<"Lattice parameters: "
                <<lattice[0]<<" x "<<lattice[1]<<" x "<<lattice[2]<<endl;
        }

        /* Read origin parameters */
        sprintf(attname,"%s","Cell origin");
        /*  Open a dataset attribute. */
        attribute_id = H5Aopen_name(dset_id, attname);
        if( attribute_id<0 ){
            cerr<<"WARNING: H5Aopen failed for "<<attname<<endl;
            clog<<0.<<"\t"<<0.<<"\t"<<0.<<" // cell origin"<<endl;
        }else{
            status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &origin[0]);
            if( status<0 ){
                cerr<<"H5Aread failed."<<endl;
                return NULL;
            }
            status = H5Aclose(attribute_id);
            if( status<0 ){
                cerr<<"H5Aclose failed."<<endl;
                return NULL;
            }
            clog<<"Cell Origin: "
                <<origin[0]<<","<<origin[1]<<","<<origin[2]<<endl;
        }

        clog<<"Mesh: "
            <<(int)dims[0]<<"\t"<<(int)dims[1]<<"\t"<<(int)dims[2]<<endl;

        
        if( datasetname[0]=='F' ){
            /* Read orbitals centers and radii */
            sprintf(attname,"%s","List of centers and radii");
            attribute_id = H5Aopen_name(dset_id, attname);
            if( attribute_id<0 ){
                cerr<<"WARNING: H5Aopen failed for "<<attname<<endl;
            }else{
                attdataspace=H5Aget_space(attribute_id);
                rank=H5Sget_simple_extent_dims(attdataspace, dimsa, maxdims);
                if( rank!=2 ){
                    cerr<<"Problem with attdataspace dimension, rank="<<rank<<endl;
                    return NULL;
                }
                if( dimsa[1]!=4 ){
                    cerr<<"Problem with attdataspace dimension, dima="<<(int)dimsa[1]<<endl;
                    return NULL;
                }
                (*dima)=dimsa[0];
                (*attributes)=new double[4*dimsa[0]];
                status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, (*attributes));
                if( status<0 ){
                    cerr<<"H5Aread failed."<<endl;
                    return NULL;
                }
                status=H5Sclose(attdataspace);
                if( status<0 ){
                    cerr<<"H5Sclose failed."<<endl;
                    return NULL;
                }
                status = H5Aclose(attribute_id);
                if( status<0 ){
                    cerr<<"H5Aclose failed."<<endl;
                    return NULL;
                }
                clog<<(int)(*dima)<<" Obital centers"<<endl;
                for (int i=0; i < (int)(*dima); i++){
                    clog<<"Obital center: ("
                        <<(*attributes)[4*i]<<","<<(*attributes)[4*i+1]<<","<<(*attributes)[4*i+2]
                        <<"), radius="<<(*attributes)[4*i+3]<<endl;
                }
            }
        }


        /* Close/release resources. */
        status = H5Dclose(dset_id);
        if( status<0 ){
            cerr<<"H5Dclose failed."<<endl;
            return NULL;
        }

        status = H5Fclose(file_id);
        if( status<0 ){
            cerr<<"H5Fclose failed."<<endl;
            return NULL;
        }
    }
    
    return data;
}

double* read_ionic_positions_hdf5(hid_t file_id)
{
    clog<<"Read ionic positions from hdf5 file"<<endl;

    // Open the dataset
    hid_t dataset_id = H5Dopen2(file_id, "/Ionic_positions",H5P_DEFAULT);
    if( dataset_id<0 ){
        cout<<"Ions::read_positions_hdf5() --- H5Dopen failed!!!"<<endl;
        return NULL;
    }

    int dim=(int)H5Dget_storage_size(dataset_id);
    
    double*  data=new double[dim];
    
    // ion i: data[3*i],data[3*i+1],data[3*i+2]
    herr_t  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                             H5P_DEFAULT, &data[0]);
    if( status<0 ){
        clog<<"read_ionic_positions_hdf5() --- H5Dread failed!!!"<<endl;
        return NULL;
    }
    H5Dclose(dataset_id);

    return data;
}

int read_atomic_numbers_hdf5(hid_t file_id, vector<int>& data)
{
    clog<<"Read atomic numbers from hdf5 file"<<endl;

    // Open the dataset
    hid_t dataset_id = H5Dopen2(file_id, "/Atomic_numbers",H5P_DEFAULT);
    if( dataset_id<0 ){
        cerr<<"Ions::read_atomic_numbers_hdf5() --- H5Dopen failed!!!"<<endl;
        return 0;
    }

    const int dim=(int)H5Dget_storage_size(dataset_id)/(int)sizeof(H5T_NATIVE_INT);
    
    data.resize(dim);

    clog<<"Read "<<dim<<" atomic numbers from hdf5 file"<<endl;
    herr_t  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                             H5P_DEFAULT, &data[0]);
    if( status<0 ){
        cerr<<"read_atomic_numbers_hdf5() --- H5Dread failed!!!"<<endl;
        return 0;
    }
    H5Dclose(dataset_id);    

    return dim;
}

int read_atomic_names_hdf5(hid_t file_id, vector<string>& data)
{
    clog<<"Read atomic names from hdf5 file"<<endl;

    // Open the dataset
    hid_t dataset_id = H5Dopen2(file_id, "/Atomic_names",H5P_DEFAULT);
    if( dataset_id<0 ){
        cerr<<"Ions::read_atomic_names_hdf5() --- H5Dopen failed!!!"<<endl;
        return 0;
    }

    // create type for strings of length MaxStrLength
    hid_t strtype =  H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, MaxStrLength);

    const int dim=(int)H5Dget_storage_size(dataset_id)/H5Tget_size(strtype);
    
    clog<<"Read "<<dim<<" atomic names from hdf5 file"<<endl;
    vector<FixedLengthString> tc(dim);
    herr_t  status = H5Dread(dataset_id, strtype, H5S_ALL, H5S_ALL,
                             H5P_DEFAULT, &tc[0]);
    if( status<0 ){
        cerr<<"read_atomic_names_hdf5() --- H5Dread failed!!!"<<endl;
        return 0;
    }

     for (vector<FixedLengthString>::const_iterator i = tc.begin (),
                                                  end = tc.end ();
         i != end;
         i++)
     {
         string t(i->mystring);
         data.push_back (t);
     }

    H5Dclose(dataset_id);    

    return dim;
}

void writeBOVheader(const string bov_filename,
                    const string data_filename,
                    const double origin[3],
                    const double ll[3],
                    const int mesh[3])
{
    const double hmesh[3]={ll[0]/mesh[0],ll[1]/mesh[1],ll[2]/mesh[2]};
    
    ofstream tfile(bov_filename.data(), ios::out);
#ifndef NDEBUG
    clog<<" Write down BOV header..."<<endl;
#endif
    tfile<<"TIME: 0.0"<<endl;
    tfile<<"DATA_FILE: "<<data_filename<<endl;
    tfile<<"DATA_SIZE: "
         <<mesh[0]<<" "<<mesh[1]<<" "<<mesh[2]<<endl;
    tfile<<"DATA_FORMAT: FLOAT"<<endl;
    tfile<<"VARIABLE: data"<<endl;
    tfile<<"DATA_ENDIAN: LITTLE"<<endl;
    tfile<<"CENTERING: zonal"<<endl;
    // shift by -0.5*h because VisIt assumes cell centered data
    tfile<<"BRICK_ORIGIN: "
         <<(origin[0]-0.5*hmesh[0])*bohr2ang<<" "
         <<(origin[1]-0.5*hmesh[1])*bohr2ang<<" "
         <<(origin[2]-0.5*hmesh[2])*bohr2ang<<endl;
    tfile<<"BRICK_SIZE: "
         <<ll[0]*bohr2ang<<" "
         <<ll[1]*bohr2ang<<" "
         <<ll[2]*bohr2ang<<endl;
}

void map3d_header(ostream& tfile,
                  const double origin[3],
                  const double ll[3])
{
#ifndef NDEBUG
    clog<<" Write down map3d header..."<<endl;
#endif
    tfile<<" 600 600              // pixels"<<endl;
    tfile<<" 0  -100    40        // viewPoint"<<endl;
    tfile<<origin[0]+0.5*ll[0]<<"\t"
         <<origin[1]+0.5*ll[1]<<"\t"
         <<origin[2]+0.5*ll[2]<<" // screenCenter"<<endl;
    tfile<<" 1     0     0        // horizontal direction"<<endl;
    tfile<<" 100 -100 100 0.2 0.5 // lightSource,ambient,diff"<<endl;
    tfile<<" 8.000000 8.000000    // hScreenSize, vScreenSize"<<endl;
    tfile<<" LightBlue 0.6 1.0    // background color"<<endl;
    tfile<<" 0.3  White           // bondradius, bond color"<<endl;
}

void write_map3d_header(ostream& tfile,
                        char* filename,
                        const double origin[3],
                        const double lattice[3])
{   
    herr_t  status;
    hid_t   file_id;

    /* check if file is in HDF5 format */
    htri_t ishdf=H5Fis_hdf5(filename);

    if( ishdf< 0){
        cerr<<"H5Fis_hdf5 unsuccessful"<<endl;
        return;
    }else if( !ishdf ){
        cerr<<"Input file "<<filename<<" not in HDF5 format. Stop."<<endl;
        return;
    }else{
        file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
        if(file_id<0){
            cerr<<"H5Fopen failed for "<<filename<<endl;
            return;
        }
    }
    
    vector<int> at_numbers;
    int n = read_atomic_numbers_hdf5(file_id, at_numbers);
    
    double* coord=read_ionic_positions_hdf5(file_id);
    
    map3d_header(tfile, origin, lattice);

    tfile<<n<<"  // natoms"<<endl;
    for(int i=0;i<n;i++)
    {
        const int at=at_numbers[i];
        tfile<<coord[3*i]<<"\t"
             <<coord[3*i+1]<<"\t"
             <<coord[3*i+2]<<"\t"
             <<ball_radii.find(at)->second<<"\t"
             <<cov_radii.find(at)->second<<"\t"
             <<colors.find(at)->second
             <<endl;
    }
    tfile<<"1  // nfunctions"<<endl;
    tfile<<origin[0]<<"\t"
         <<origin[1]<<"\t"
         <<origin[2]<<"  // origin"<<endl;
    tfile<<lattice[0]<<"\t"
         <<lattice[1]<<"\t"
         <<lattice[2]<<"  // lattice"<<endl;
    tfile<<"SteelBlue         // color"<<endl;
    tfile<<"0.0004  0.        // flevel,thresh"<<endl;
}


// write list of atoms in XYZ format
void writeAtomsXYZ(const string xyz_filename,
                   char* filename)
{   
    ofstream tfile(xyz_filename.data(), ios::out);

    herr_t  status;
    hid_t   file_id;

    /* check if file is in HDF5 format */
    htri_t ishdf=H5Fis_hdf5(filename);

    if( ishdf< 0){
        cerr<<"H5Fis_hdf5 unsuccessful"<<endl;
        return;
    }else if( !ishdf ){
        cerr<<"Input file "<<filename<<" not in HDF5 format. Stop."<<endl;
        return;
    }else{
        file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
        if(file_id<0){
            cerr<<"H5Fopen failed for "<<filename<<endl;
            return;
        }
    }
    
    set<string> atomicsp;
    atomicsp.insert("H");
    atomicsp.insert("Li");
    atomicsp.insert("Be");
    atomicsp.insert("B");
    atomicsp.insert("C");
    atomicsp.insert("N");
    atomicsp.insert("O");
    atomicsp.insert("F");
    atomicsp.insert("Na");
    atomicsp.insert("Mg");
    atomicsp.insert("Al");
    atomicsp.insert("Si");
    atomicsp.insert("P");
    atomicsp.insert("S");
    atomicsp.insert("Cl");
    atomicsp.insert("K");
    atomicsp.insert("Ca");
    atomicsp.insert("Ni");
    atomicsp.insert("Cu");
    atomicsp.insert("Ga");
    atomicsp.insert("Ge");
    atomicsp.insert("Au");
    
    vector<string> at_names;
    int n = read_atomic_names_hdf5(file_id, at_names);
    
    double* coord=read_ionic_positions_hdf5(file_id);
    
    tfile<<n<<endl; // number of atoms
    tfile<<endl;
    for(int i=0;i<n;i++)
    {
        string sp(at_names[i].substr(0,2));
        set<string>::const_iterator p;
        p=atomicsp.find(sp);
        if( p==atomicsp.end() ){
            sp=at_names[i].substr(0,1);
        }
        tfile<<sp<<"\t"
             <<coord[3*i]*bohr2ang<<"\t"
             <<coord[3*i+1]*bohr2ang<<"\t"
             <<coord[3*i+2]*bohr2ang<<"\t"
             <<endl;
    }
}
   
int main (int argc, char **argv)
{
    char    datasetname[50];
    char    h5filename[50];
    char    buf[255];

    hsize_t dims[3];
    hsize_t dima=1;
    double  lattice[3];
    double  origin[3];
    
    double* attributes=NULL;
    
    bool tmap3d = false;
    bool tbov   = false;
    if ( argc>1 ){
        tmap3d = !strcmp(argv[1],"-map3d");
        tbov   = !strcmp(argv[1],"-bov");
    }
    
    int index=1;
    if( tmap3d || tbov )index=2;
    strcpy(h5filename, argv[index]);
    strcpy(datasetname, argv[index+1]);

    clog<<"Dataset: "<<datasetname<<endl;
    double* data
        =get_function(h5filename, datasetname, dims, origin, lattice, 
                      &attributes, &dima);
    if( data==NULL ){
        cerr<<"Read failed!";
        return -1;
    }
        
    string base_filename(h5filename);
    // remove extension
    size_t found=base_filename.find_last_of('.');
    size_t end=base_filename.size();
    base_filename.erase(found,end-found);
    
    string output_data_filename(base_filename);
    output_data_filename.append("_");
    output_data_filename.append(datasetname);
    output_data_filename.append(".dat");
    clog<<"output_data_filename="<<output_data_filename<<endl;
      
    const int dim[3]={(int)dims[0],
                      (int)dims[1],
                      (int)dims[2]};

    const int incx=dim[1]*dim[2];
    const int incy=dim[2];

    // header,including atoms list
    if( tmap3d ){
        set_maps();
        string header_filename(base_filename);
        header_filename.append("_header.map3d");
        ofstream tfile(header_filename.data(), ios::out);
        write_map3d_header(tfile,h5filename, origin, lattice);
    }else{
        string xyz_filename(base_filename);
        xyz_filename.append("_atoms.xyz");
    
        writeAtomsXYZ(xyz_filename,h5filename);
        
        if( tbov ){
            string bov_filename(base_filename);
            bov_filename.append(".bov");
            clog<<"bov_filename="<<bov_filename<<endl;
            writeBOVheader(bov_filename,
                           output_data_filename,
                           origin, lattice, dim);
        }
    }            
    
    // function
    if( tbov ){
        const int inczf=dim[1]*dim[0];
        const int incyf=dim[0];
        clog<<" Write BOV data..."<<endl;
        float* fdata=new float[dim[0]*dim[1]*dim[2]];
        //transpose data
        for (int i=0; i < dim[0]; i++){
        for (int j=0; j < dim[1]; j++){
        for (int k=0; k < dim[2]; k++){
            int index1=k*inczf+j*incyf+i;
            int index2=i*incx+j*incy+k;
            fdata[index1]=(float)data[index2]; 
        }}}
        ofstream binary(output_data_filename.data(), ios::binary);
        binary.write((char*)fdata, sizeof(float)*dim[0]*dim[1]*dim[2]);    
    }else{
        clog<<" Write data..."<<endl;
        ofstream tfile(output_data_filename.data(), ios::out);
        if( !tmap3d ){
            tfile<<origin[0]<<"\t"
                 <<origin[1]<<"\t"
                 <<origin[2]<<"\t"
                 <<origin[0]+lattice[0]<<"\t"
                 <<origin[1]+lattice[1]<<"\t"
                 <<origin[2]+lattice[2]<<"  // cell corners"<<endl;
        }
        tfile<<dim[0]<<"\t"<<dim[1]<<"\t"<<dim[2]
             <<" // mesh"<<endl;
        
        for (int i=0; i < dim[0]; i++){
        for (int j=0; j < dim[1]; j++){
        for (int k=0; k < dim[2]; k++){
            int row=i*incx+j*incy+k;
            tfile<<data[row]<<endl; 
        }}}
    }

    delete[] data;
    delete[] attributes;

    return 0;
}     
