# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
import sys
import argparse
import h5py
import numpy as np


# Read dataset in hdf5 file.
# Returns an array containing the data.
def get_function(filename, datasetname, dims):

    # Check If File is in HDF5 Format
    try:
        ishdf = h5py.is_hdf5(filename)
    except Exception:
        print('\nh5py.is_hdf5 unsucessful')
        return None
    
    # If File is not in HDF5 Format, Stop
    if( not(ishdf) ):
        print('\nInput File ' + filename + ' not in HDF5 Format. Stop.')
        return None

    # If Everything Goes Fine, Proceed
    else:
        
        # Open HDF5 File
        try:
            file_id = h5py.h5f.open(bytes(filename, encoding='utf-8'),
                          h5py.h5f.ACC_RDONLY, h5py.h5p.DEFAULT)
        except Exception:
            print('\nHDF5 File: ' + filename + ' Failed to Open')
            return None

        # Open Dataset  
        try:
            dset_id = h5py.h5d.open(file_id, bytes(datasetname, encoding='utf-8'))
        except Exception:
            print('\nHDF5 Dataset: ' + datasetname + ' Failed to Open')
            return None

        # Copy of Dataspace for Dataset 
        try:
            filespace = dset_id.get_space()

        except Exception:
            print('\ndset_id.get_space() Failed.')
            return None

        # Get Dataspace Dimension
        ndims = filespace.get_simple_extent_ndims()
        # If Dataspace Dimension is not 3, Stop.
        if( not(ndims == 3) ):
            print('\nProblem with Dataspace Dimension, ndims = ' + str(ndims))
            return None
        
        # Shape of Dataspace (dims)
        dims = dims.tolist()
        dims = filespace.get_simple_extent_dims()

        print('Dataspace: Dimensions ' + str( int(dims[0]) ) + ' x '
                                       + str( int(dims[1]) ) + ' x '
                                       + str( int(dims[2])) )

        print('Size: ' + str( int(dims[0] * dims[1] * dims[2]) ))

        # If Size < 1, Stop.
        if( int( dims[0] * dims[1] * dims[2] ) < 1 ):
            return None

        # Read data -> data
        data = np.array(0.0, h5py.h5t.NATIVE_FLOAT)
        data.resize( int(dims[0] * dims[1] * dims[2]) , refcheck = False)

        # Dump Data into Numpy Array (data)
        try:
            status = dset_id.read(h5py.h5s.ALL, h5py.h5s.ALL, data)
        except Exception:
            print('\ndataset_id.read Failed.')
            return 0

    return ( data, dims )


''' MAIN '''
# USAGE:
# python hdf5toMM.py file.hdf5

def main():

    h5filename = sys.argv[1]
    field = 'Function'

    # Remove File Extension ( .hdf5 )
    base_filename = h5filename.split('.')[0].strip()
    base_filename = base_filename.split('/')[-1]

    # Use base_filename to make .dat filename
    output_data_filename = base_filename + '.dat'
    print('\noutput_data_filename = ' + output_data_filename)

    columns = []

    i = 0
    while i<1000:
      # Get data - Call get_function
      number = str(i)
      while len(number)<4:
        number = '0'+number

      datasetname = field + number
      print('\nDataset: ' + datasetname)

      dims = np.arange(0, dtype = h5py.h5t.NATIVE_INT32)  # Turns Into a TUPLE

      try:
        column, dims = get_function(h5filename, datasetname, dims)
      except Exception:
        print('\nRead Failed. \nEither the HDF5 File ' +
              'or the Dataset are not Present. Stop.\n')
        break

      # If data Empty, Stop.
      if( column is None or dims is None ):
        print('\nRead Failed.')
        return -1

      #add data just read as a column in list of columns
      columns.append(column)

      dim = [ int(dims[0]), int(dims[1]), int(dims[2]) ]

      i = i+1

    #build numpy 2d array from all the columns
    matrix = columns[0]
    for i in range(len(columns)-1):
      matrix = np.column_stack((matrix,columns[i+1]))

    print('\nWrite data...\n')

    nrows = dim[0]*dim[1]*dim[2]
    ncols = len(columns)

    np.savetxt('matrix.dat', matrix, delimiter='\t', fmt='%le', header=str(nrows) + '\t' + str(ncols))

    return 0

# Executes Main Function
if __name__ == '__main__':

    main()
