"""
checkBonds.py

Brief:  Check that all D atoms are bonded to one
other D atom in d144.

Author:  Ian S. Dunn, Lawrence Livermore National Laboratory
"""

# Import modules.
import sys
import numpy as np

# Load .xyz movie file.
with open(sys.argv[1], 'r') as f:
    lines = f.readlines()

ang2bohr=1.8897269
bohr2ang=1./ang2bohr

box = np.array([42.4813, 21.2406, 21.2406])  * bohr2ang

current_timestep = []
num_timesteps = 0

for i in range(int(sys.argv[2]), len(lines)):

    if 'H' in lines[i]:
        coords = lines[i].rstrip('\n').split()[1:4]
        coords = np.array([float(coord) for coord in coords])
        current_timestep.append(coords)

    elif ' ' not in lines[i] and lines[i] != '\n':
        if len(current_timestep) == 0:
            continue

        current_timestep = np.array(current_timestep)

        # Check for points with no neighbors.
        neighborless = 0

        for k in range(len(current_timestep)):

            smallest_distance = 100000000

            for j in range(len(current_timestep)):

                if j == k:
                    continue

                distance = current_timestep[k] - current_timestep[j]
                
                for l in range(3):
                    distance[l] = min(abs(distance[l]), abs(abs(distance[l]) - box[l]))

                distance = np.linalg.norm(distance)

                if distance < smallest_distance:
                    smallest_distance = distance

            if smallest_distance > 4.0:
                neighborless += 1
                print k

        print "Time step "+str(np.ceil(i / 218.0))+": "+str(neighborless)
        num_timesteps += 1
        current_timestep = []
