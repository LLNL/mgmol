import numpy as np
import os

def process_data(s_l1, s_l2, s_theta, N_l, N_theta, rdim):
    try:
        data_dir = f'data/{s_l1}_{s_l2}_{s_theta}'
        offline_file = os.path.join(data_dir, 'offline_PinnedH2O.out')
        log_file = f'{s_l1}_{s_l2}_{s_theta}_{N_l}_{N_theta}_{rdim}.log'

        f_O1_fom = None
        f_H1_fom = None
        f_H2_fom = None

        if not os.path.exists(offline_file):
            print(f"Error: File not found: {offline_file}")
            return None

        with open(offline_file, 'r') as f:
            eks_lines = [line for line in f if " SC ENERGY      [Ha]" in line]
            if eks_lines:
                Eks_fom = float(eks_lines[-1].split()[-1])
            else:
                print(f"Error: Could not find 'SC ENERGY    [Ha]' line in {offline_file}")
                return None

        with open(offline_file, 'r') as f:
            for line in f:
                if "O1" in line:
                    f_O1_fom = np.array(line.split()[6:9], dtype=float)
                elif "H1" in line:
                    f_H1_fom = np.array(line.split()[5:8], dtype=float)
                elif "H2" in line:
                    f_H2_fom = np.array(line.split()[5:8], dtype=float)

                if f_O1_fom is not None and f_H1_fom is not None and f_H2_fom is not None:
                    break # Exit the loop once we found all lines

        if f_O1_fom is None or f_H1_fom is None or f_H2_fom is None:
            print(f"Error: Could not find O1, H1, or H2 lines in {offline_file}")
            return None

        f_H2_rom = None
        f_O1_rom = None
        f_H1_rom = None

        if not os.path.exists(log_file):
            print(f"Error: File not found: {log_file}")
            return None

        with open(log_file, 'r') as f:
            lines = f.readlines()

        eks_index = -1
        for i in range(len(lines) - 1, -1, -1):  # Search from the end
            if "Eks:" in lines[i]:
                eks_index = i
                break

        if eks_index == -1:
            print(f"Error: Could not find 'Eks:' line in {log_file}")
            return None

        Eks_rom = float(lines[eks_index].split()[-1])  # Get the last element

        force_index = -1
        for i in range(len(lines) - 1, -1, -1):  # Search from the end
            if "Forces:" in lines[i]:
                force_index = i
                break

        if force_index == -1:
            print(f"Error: Could not find 'Forces:' line in {log_file}")
            return None

        if force_index + 3 >= len(lines):
            print(f"Error: Not enough lines after 'Forces:' in {log_file}")
            return None

        f_H2_rom = np.array(lines[force_index + 1].split(), dtype=float)
        f_O1_rom = np.array(lines[force_index + 2].split(), dtype=float)
        f_H1_rom = np.array(lines[force_index + 3].split(), dtype=float)

        return Eks_fom, f_O1_fom, f_H1_fom, f_H2_fom, Eks_rom, f_O1_rom, f_H1_rom, f_H2_rom

    except Exception as e:
        print(f"An error occurred: {e}")
        return None

s_l1 = "1.02"
s_l2 = "0.98"
s_theta = "2.0"
N_l = 5
N_theta = 2
rdim = 8

output_filename = f'energy_and_forces_{s_l1}_{s_l2}_{s_theta}_{N_l}_{N_theta}_{rdim}.out'

with open(output_filename, 'w') as outfile:
    # Print parameters to file
    print("s_l1:", s_l1, file=outfile)
    print("s_l2:", s_l2, file=outfile)
    print("s_theta:", s_theta, file=outfile)
    print("N_l:", N_l, file=outfile)
    print("N_theta:", N_theta, file=outfile)
    print("rdim:", rdim, file=outfile)

    results = process_data(s_l1, s_l2, s_theta, N_l, N_theta, rdim)

    if results:
        Eks_fom, f_O1_fom, f_H1_fom, f_H2_fom, Eks_rom, f_O1_rom, f_H1_rom, f_H2_rom = results

        print("Eks_fom:", Eks_fom, file=outfile)
        print("Eks_rom:", Eks_rom, file=outfile)

        print("f_O1_fom:", f_O1_fom, file=outfile)
        print("f_O1_rom:", f_O1_rom, file=outfile)

        print("f_H1_fom:", f_H1_fom, file=outfile)
        print("f_H1_rom:", f_H1_rom, file=outfile)

        print("f_H2_fom:", f_H2_fom, file=outfile)
        print("f_H2_rom:", f_H2_rom, file=outfile)

        def calculate_differences(fom, rom, name):
            diff = fom - rom
            abs_diff = np.linalg.norm(diff)
            rel_diff = abs_diff / np.linalg.norm(fom) if np.linalg.norm(fom)!= 0 else float('inf')

            print(f"Absolute difference in {name}:", abs_diff, file=outfile)
            print(f"Relative difference in {name}:", rel_diff, file=outfile)

        calculate_differences(Eks_fom, Eks_rom, "Eks")
        calculate_differences(f_O1_fom, f_O1_rom, "f_O1")
        calculate_differences(f_H1_fom, f_H1_rom, "f_H1")
        calculate_differences(f_H2_fom, f_H2_rom, "f_H2")

    else:
        print("Error occurred during data processing. Differences cannot be calculated.", file=outfile)

