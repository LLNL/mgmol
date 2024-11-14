import numpy as np
import os

# coords.in
O1 = np.array([0.00, 0.00, 0.00])
ref_H1 = np.array([-0.45, 1.42, -1.07])
ref_H2 = np.array([-0.45, -1.48, -0.97])

# factors and increments for bond lengths and bond angle
bondlength1_factor = np.linspace(0.95, 1.05, 11)
bondlength2_factor = np.linspace(0.95, 1.05, 11)
bondangle_increment = np.linspace(-5, 5, 11)

# output directory
output_dir = "PinnedH2O_3dof_coords"

# utilities
def calculate_bondlength(atom1, atom2):
    return np.linalg.norm(atom1 - atom2)

def calculate_bondangle(atom1, atom2, atom3):
    vector1 = atom1 - atom2
    vector2 = atom3 - atom2
    dot_product = np.dot(vector1, vector2)
    magnitude_product = np.linalg.norm(vector1) * np.linalg.norm(vector2)
    angle_rad = np.arccos(dot_product / magnitude_product)
    angle_deg = np.degrees(angle_rad)
    return angle_deg

# Rodrigues' rotation formula
def rotation_matrix(axis, angle_degrees):
    angle = np.radians(angle_degrees)
    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)
    ux, uy, uz = axis
    return np.array([
        [cos_theta + ux**2 * (1 - cos_theta), ux * uy * (1 - cos_theta) - uz * sin_theta, ux * uz * (1 - cos_theta) + uy * sin_theta],
        [uy * ux * (1 - cos_theta) + uz * sin_theta, cos_theta + uy**2 * (1 - cos_theta), uy * uz * (1 - cos_theta) - ux * sin_theta],
        [uz * ux * (1 - cos_theta) - uy * sin_theta, uz * uy * (1 - cos_theta) + ux * sin_theta, cos_theta + uz**2 * (1 - cos_theta)]
    ])

# generation
os.makedirs(output_dir, exist_ok=True)

ref_bondlength1 = calculate_bondlength(ref_H1, O1)
ref_bondlength2 = calculate_bondlength(ref_H2, O1)
ref_bondangle = calculate_bondangle(ref_H1, O1, ref_H2)

normal_vector = np.cross(ref_H1, ref_H2)
normal_unit_vector = normal_vector / np.linalg.norm(normal_vector)

for d_bondangle in bondangle_increment:
    Q = rotation_matrix(normal_unit_vector, d_bondangle)
    Q_ref_H2 = np.dot(Q, ref_H2)
    for f_bondlength1 in bondlength1_factor:
        for f_bondlength2 in bondlength2_factor:
            H1 = f_bondlength1 * ref_H1
            H2 = f_bondlength2 * Q_ref_H2
            filename = f"{output_dir}/coords_{f_bondlength1:.2f}_{f_bondlength2:.2f}_{d_bondangle}.in"
            with open(filename, "w") as file:
                file.write(f"O1  1   {O1[0]:.2f}    {O1[1]:.2f}    {O1[2]:.2f}  0\n")
                file.write(f"H1  2   {H1[0]:.2f}    {H1[1]:.2f}    {H1[2]:.2f}  1\n")
                file.write(f"H2  2   {H2[0]:.2f}    {H2[1]:.2f}    {H2[2]:.2f}  1\n")
