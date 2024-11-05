import numpy as np

O1 = np.array([0.00, 0.00, 0.00])
H1 = np.array([-0.45, 1.42, -1.07])
H2 = np.array([-0.45, -1.48, -0.97])

def calculate_bondlength(atom1, atom2):
    return np.linalg.norm(atom1 - atom2)

def calculate_bondangle(atom1, atom2, atom3, radian):
    vector1 = atom1 - atom2
    vector2 = atom3 - atom2
    dot_product = np.dot(vector1, vector2)
    magnitude_product = np.linalg.norm(vector1) * np.linalg.norm(vector2)
    angle = np.arccos(dot_product / magnitude_product)
    if not radian:
        angle = np.degrees(angle)
    return angle

def rotation_matrix(axis, angle):
    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)
    ux, uy, uz = axis
    return np.array([
        [cos_theta + ux**2 * (1 - cos_theta), ux * uy * (1 - cos_theta) - uz * sin_theta, ux * uz * (1 - cos_theta) + uy * sin_theta],
        [uy * ux * (1 - cos_theta) + uz * sin_theta, cos_theta + uy**2 * (1 - cos_theta), uy * uz * (1 - cos_theta) - ux * sin_theta],
        [uz * ux * (1 - cos_theta) - uy * sin_theta, uz * uy * (1 - cos_theta) + ux * sin_theta, cos_theta + uz**2 * (1 - cos_theta)]
    ])

plane_normal = np.cross(H2, H1)
plane_normal = plane_normal / np.linalg.norm(plane_normal)

target_plane_normal = np.array([0, 0, 1])
axis_to_align = np.cross(plane_normal, target_plane_normal)
axis_to_align /= np.linalg.norm(axis_to_align)
angle_to_align = np.arccos(np.clip(np.dot(plane_normal, target_plane_normal), -1.0, 1.0))

rot_matrix_align_plane = rotation_matrix(axis_to_align, angle_to_align)
H1_rotated = np.dot(rot_matrix_align_plane, H1)
H2_rotated = np.dot(rot_matrix_align_plane, H2)

bondlength1 = calculate_bondlength(H1, O1)
bondlength2 = calculate_bondlength(H2, O1)
bondangle = calculate_bondangle(H1, O1, H2, False)

print('Original system')
print(f'H1 = {H1}')
print(f'H2 = {H2}')
print(f'Bondlength of O1-H1 = {bondlength1}')
print(f'Bondlength of O1-H2 = {bondlength2}')
print(f'Angle between O1-H1 and O1-H2 = {bondangle}')

bondlength1 = calculate_bondlength(H1_rotated, O1)
bondlength2 = calculate_bondlength(H2_rotated, O1)
bondangle = calculate_bondangle(H1_rotated, O1, H2_rotated, False)

print('Aligned system')
print(f'H1 = {H1_rotated}')
print(f'H2 = {H2_rotated}')
print(f'Bondlength of O1-H1 = {bondlength1}')
print(f'Bondlength of O1-H2 = {bondlength2}')
print(f'Angle between O1-H1 and O1-H2 = {bondangle}')
