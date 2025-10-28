import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import argparse

parser = argparse.ArgumentParser(description="Calculate and plot differences in energies and forces.")

parser.add_argument("--N_l", type=int, default=2, help="Sampling frequency of bond length points")
parser.add_argument("--N_theta", type=int, default=2, help="Sampling frequency of bond angle points")
parser.add_argument("--rdim", type=int, default=18, help="Dimension of the ROM basis")

args = parser.parse_args()

N_l = args.N_l
N_theta = args.N_theta
rdim = args.rdim

bondlength_min = 0.95
bondlength_max = 1.05
bondlength_num_increments = 10

bondangle_min = -5.0
bondangle_max = 5.0
bondangle_num_increments = 10

output_dir = f'results_{N_l}_{N_theta}_{rdim}'

def process_data(s_l1, s_l2, s_theta, N_l, N_theta, rdim):
    try:
        data_dir = f'/usr/workspace/nlrom/MGmol/PinnedH2O_3DOF/data_8/{s_l1}_{s_l2}_{s_theta}'
        offline_file = os.path.join(data_dir, 'offline_PinnedH2O.out')
        log_file = f'{output_dir}/{s_l1}_{s_l2}_{s_theta}.log'

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

Eks_diff_reproductive = []
f_O1_diff_reproductive = []
f_H1_diff_reproductive = []
f_H2_diff_reproductive = []

Eks_diff_predictive = []
f_O1_diff_predictive = []
f_H1_diff_predictive = []
f_H2_diff_predictive = []

for i in range(bondlength_num_increments + 1):
    bondlength_one = round(bondlength_min + i * (bondlength_max - bondlength_min) / bondlength_num_increments, 2)
    for j in range(i + 1):
        bondlength_two = round(bondlength_min + j * (bondlength_max - bondlength_min) / bondlength_num_increments, 2)
        for k in range(bondangle_num_increments + 1):
            bondangle = round(bondangle_min + k * (bondangle_max - bondangle_min) / bondangle_num_increments, 1)

            s_l1 = f"{bondlength_one:.2f}"
            s_l2 = f"{bondlength_two:.2f}"
            s_theta = f"{bondangle:.1f}"
            tag = f'{s_l1}_{s_l2}_{s_theta}'

            output_filename = f'{output_dir}/energy_and_forces_{tag}.out'

            with open(output_filename, 'w') as outfile:
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
                        abs_diff = np.linalg.norm(fom - rom)
                        rel_diff = abs_diff / np.linalg.norm(fom) if np.linalg.norm(fom)!= 0 else float('inf')
                        print(f"Absolute difference in {name}:", abs_diff, file=outfile)
                        print(f"Relative difference in {name}:", rel_diff, file=outfile)
                        return abs_diff

                    Eks_diff = calculate_differences(Eks_fom, Eks_rom, "Eks")
                    f_O1_diff = calculate_differences(f_O1_fom, f_O1_rom, "f_O1")
                    f_H1_diff = calculate_differences(f_H1_fom, f_H1_rom, "f_H1")
                    f_H2_diff = calculate_differences(f_H2_fom, f_H2_rom, "f_H2")

                    if i * N_l % bondlength_num_increments == 0 and j * N_l % bondlength_num_increments == 0 and k * N_theta % bondangle_num_increments == 0:
                        Eks_diff_reproductive.append(Eks_diff)
                        f_O1_diff_reproductive.append(f_O1_diff)
                        f_H1_diff_reproductive.append(f_H1_diff)
                        f_H2_diff_reproductive.append(f_H2_diff)
                    else:
                        Eks_diff_predictive.append(Eks_diff)
                        f_O1_diff_predictive.append(f_O1_diff)
                        f_H1_diff_predictive.append(f_H1_diff)
                        f_H2_diff_predictive.append(f_H2_diff)
                else:
                    print("Error occurred during data processing. Differences cannot be calculated.", file=outfile)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
# (Assuming the rest of your imports and script setup remain the same)

def plot_histogram(data, quantity, test_case):
    # Increased figure size for better spacing (Retained)
    plt.figure(figsize=(10, 7.5))
    
    # Calculate the histogram data first
    counts, bin_edges, _ = plt.hist(data, bins=20, color='skyblue', edgecolor='black')

    if quantity == "Eks":
        quantity_name = "absolute difference in total energy"
    elif quantity.startswith("f_H1"):
        quantity_name = r"$\| \mathbf{F}_1 - \widetilde{\mathbf{F}}_1 \|_2$"
    elif quantity.startswith("f_H2"):
        quantity_name = r"$\| \mathbf{F}_2 - \widetilde{\mathbf{F}}_2 \|_2$"
    elif quantity.startswith("f_"):
        quantity_name = f"magnitude of difference in force on {quantity[2:]}"
    else:
        raise ValueError("Invalid input quantity")

    plt.title(f'Histogram of {quantity_name}', fontsize=22)
    plt.xlabel('Difference', fontsize=22)
    plt.ylabel('Frequency', fontsize=22)

    # 1. Y-axis: Max 5 labels and starting at 0
    max_frequency = np.max(counts)
    
    # Calculate the ideal step size to yield at most 5 labels
    # The ceiling function ensures we cover the max value.
    # The number of steps will be up to 5 (including 0).
    max_y = np.ceil(max_frequency / 5.0) * 5.0  # Round up to the nearest multiple of 5
    if max_y == 0:
        max_y = 1 # Handle case where all counts are 0
    
    # Calculate the step size to get at most 5 ticks (excluding 0)
    # The target number of intervals is 4-5.
    num_intervals = 4 
    y_step = np.ceil(max_frequency / num_intervals)
    
    # We round the step size to a "nice" number (e.g., 1, 2, 5, 10, 20, 50, 100...)
    # This is a common plotting requirement. We'll use a simple heuristic for large numbers.
    if y_step <= 5:
        y_step = max(1, y_step)
    elif y_step <= 10:
        y_step = 10
    elif y_step <= 25:
        y_step = 25
    elif y_step <= 50:
        y_step = 50
    elif y_step <= 100:
        y_step = 100
    else:
        # For very high numbers, round to nearest 50 or 100
        y_step = np.ceil(y_step / 100.0) * 100.0
        
    y_step = int(y_step) # Ensure it's an integer step

    # Generate y-axis ticks
    # Use arange to get ticks starting at 0 up to max_frequency plus the step
    y_ticks = np.arange(0, np.ceil(max_frequency) + y_step, y_step)
    
    # Adjust to ensure the highest tick is not far above the max bar height
    y_ticks = y_ticks[y_ticks <= np.ceil(max_frequency) + y_step * 0.5]
    if y_ticks[-1] < np.ceil(max_frequency):
        y_ticks = np.append(y_ticks, y_ticks[-1] + y_step)
        
    # Ensure there is a 0 tick if it was somehow missed
    if y_ticks[0] != 0:
        y_ticks = np.insert(y_ticks, 0, 0)

    # Use unique ticks and convert to int for cleaner labels
    y_ticks = np.unique(y_ticks.astype(int))

    plt.yticks(y_ticks, fontsize=22)
    # Set y-limits from 0 up to the final highest tick
    plt.ylim(0, y_ticks[-1] + 0.5) 

    # 2. X-axis tick and scientific notation font control
    min_val, max_val = np.min(data), np.max(data)
    plt.xlim(min_val - 0.05 * (max_val - min_val), max_val + 0.05 * (max_val - min_val))
    
    # Use 6 ticks for less dense x-axis
    num_ticks = 6
    xticks = np.linspace(min_val, max_val, num_ticks)
    plt.xticks(xticks, fontsize=22)
    
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((0, 0))
    plt.gca().xaxis.set_major_formatter(formatter)

    # Adjusting the scientific notation font size (the 'x10^-4' part)
    plt.gca().ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    ax = plt.gca()
    
    # Explicitly set the font size of the scientific notation exponent
    try:
        # Matplotlib's way to find and set the exponent text's font size
        # This is a more robust way to increase the font size of the exponent
        ax.xaxis.get_offset_text().set_fontsize(22) 
    except Exception as e:
        print(f"Error setting scientific notation exponent font size: {e}")

    # Stats box positioning (Retained)
    total_count = len(data)
    mean_val = np.mean(data)
    stats_text = (f"Total {test_case} cases: {total_count}\n"
                  f"Mean: {mean_val:.3e}")
    
    # Placed in the top right corner
    plt.text(0.95, 0.95, stats_text, transform=plt.gca().transAxes, 
             fontsize=22, verticalalignment='top', horizontalalignment='right',
             bbox=dict(facecolor='white', alpha=0.7, edgecolor='black'))

    plt.tight_layout()
    plt.savefig(f"{output_dir}/{quantity}_difference_histogram_{test_case}.png") # Uncomment in final script

plot_histogram(Eks_diff_reproductive, "Eks", "reproductive")
plot_histogram(f_O1_diff_reproductive, "f_O1", "reproductive")
plot_histogram(f_H1_diff_reproductive, "f_H1", "reproductive")
plot_histogram(f_H2_diff_reproductive, "f_H2", "reproductive")

plot_histogram(Eks_diff_predictive, "Eks", "predictive")
plot_histogram(f_O1_diff_predictive, "f_O1", "predictive")
plot_histogram(f_H1_diff_predictive, "f_H1", "predictive")
plot_histogram(f_H2_diff_predictive, "f_H2", "predictive")
