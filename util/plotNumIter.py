# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
"""
plotNumIter.py

Description: Plots number of inner iterations at each 
  MD step in an mgmol output file.

Author: Ian Dunn, Lawrence Livermore National Laboratory
"""

import sys
import numpy
import matplotlib.pyplot as plt

if __name__ == '__main__':


    colors = ['b', 'k', 'r', 'm', 'b:', 'k:', 'r:']

    for i in range(1, len(sys.argv)):

        # Name of mgmol output file.
        file_name = sys.argv[i]
    
        # Read mgmol output file.
        with open(file_name, 'r') as f:
            lines = f.readlines()
    
        # Counters for SCF and MD iterations.
        md_iteration = 0
        scf_counter = 0
    
        # List of number of inner iterations at each MD step.
        scf_steps = []
    
        # Indicates whether we are in the part of the output file
        # documenting the MD run.
        start = False
    
        # Iterate over lines in file.
        for line in lines:
            
            # Find start of MD run.
            if 'Verlet MD' in line:
                start = True
    
            # Determine number of inner iterations at each MD step.
            if start:
    
                if '%%    0 SC ENERGY' in line:
                    scf_steps.append(int(scf_counter))
                    md_iteration = line.split()[1]
    
                if 'SC ENERGY =' in line:
                    scf_counter = line.split()[1]
    
        # Plot results.
        plt.plot(scf_steps, colors[i-1], label=file_name.replace('/mgmol.out', ''))
    #plt.title(file_name)
    plt.xlabel('MD Iteration')
    plt.ylabel('# of SCF Iterations')
    plt.legend()
    plt.savefig('numIter.png', dpi=100)
    plt.show()
