import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
import h5py
import matplotlib.animation as animation
import sys

verbose = 1
vectorize = 0
temperature = 300.0 # Parameter
dt = 9.6
n_procs = 64
b = 32
Nt = 100

tag = f"H2O_64_{temperature}"
input_data_path = f"/p/lustre2/cheung26/H2O_64_9.6fs_new/data/{tag}"
output_data_path = f"/p/lustre2/cheung26/H2O_64_9.6fs_new/data/{tag}/postprocess"
n_blocks = round(n_procs ** (1/3))
mesh = b * n_blocks

if verbose:
    print(f"input_data_path = {input_data_path}")
    print(f"output_data_path = {output_data_path}")
    print(f"n_blocks = {n_blocks}")
    print(f"mesh = {mesh}")

os.makedirs(output_data_path, exist_ok=True)

for nt in range(0, Nt + 1):
    print(f"\nProcessing time step nt = {nt}")
    if nt == 0:
        filename = f"{tag}/orbital"
    else:
        filename = f"MD_mdstep{nt}"
            
    block_snapshots = []
    
    for rank in range(n_procs):
        snapshot_filename = filename + f"_snapshot.{rank:06}"
        snapshot_path = os.path.join(input_data_path, snapshot_filename)
        
        if not os.path.exists(snapshot_path):
            print(f"Warning: Snapshot file not found: {snapshot_path}. Skipping rank {rank}.")
            continue
            
        try:
            snapshot_data = h5py.File(snapshot_path, 'r')
        except Exception as e:
            print(f"Error opening HDF5 file {snapshot_path}: {e}. Skipping rank {rank}.")
            if 'snapshot_data' in locals() and snapshot_data is not None:
                snapshot_data.close()
            continue

        if verbose:
            print(f"  Processing rank = {rank}")

        if f'snapshot_matrix_num_cols' not in snapshot_data.keys():
            print(f"  Error: Key 'snapshot_matrix_num_cols' not found in {snapshot_path}. Skipping rank {rank}.")
            snapshot_data.close()
            continue

        n_cols_snapshot = snapshot_data[f'snapshot_matrix_num_cols'][0]
        
        if f'snapshot_matrix' not in snapshot_data.keys():
            print(f"  Error: Key 'snapshot_matrix' not found in {snapshot_path}. Skipping rank {rank}.")
            snapshot_data.close()
            continue
        
        block_snapshot = np.array(snapshot_data[f'snapshot_matrix'])
        snapshot_data.close()

        if verbose:
            print("  Block snapshot raw shape: ", block_snapshot.shape)

        expected_elements = b * b * b * n_cols_snapshot
        if block_snapshot.size != expected_elements:
            print(f"  Error: block_snapshot size ({block_snapshot.size}) does not match expected elements ({expected_elements}) for rank {rank}. Skipping.")
            continue

        block_snapshot = np.reshape(block_snapshot, (b, b, b, n_cols_snapshot))
        if verbose:
            print("  Block snapshot reshaped shape: ", block_snapshot.shape)

        block_snapshots.append(block_snapshot)
            
    if len(block_snapshots) != n_procs:
        print(f"Warning: Only {len(block_snapshots)} out of {n_procs} snapshots were successfully loaded for nt = {nt}. Cannot reconstruct full snapshot. Skipping this time step.")
        continue

    try:
        full_snapshots = np.concatenate(
            [ np.concatenate(
                [ np.concatenate(block_snapshots[(i*n_blocks+j)*n_blocks:(i*n_blocks+j+1)*n_blocks], axis=2)
                for j in range(n_blocks)], axis=1)
            for i in range(n_blocks)], axis=0
        )
    except Exception as e:
        print(f"Error concatenating blocks for nt = {nt}: {e}. This might indicate an issue with block_snapshots list or dimensions after loading.")
        continue

    if vectorize:
        full_snapshots = np.reshape(full_snapshots, (mesh**3, n_cols_snapshot))

    if verbose:
        print("  Full snapshots shape: ", full_snapshots.shape)
            
    density = np.sum(np.square(full_snapshots), axis=-1)

    if verbose:
        print("  Density shape: ", density.shape)

    try:
        pickle.dump(full_snapshots.astype('float32'), open(f"{output_data_path}/snapshots_{nt}.p",'wb'))
        pickle.dump(density.astype('float32'), open(f"{output_data_path}/density_{nt}.p", 'wb'))
    except Exception as e:
        print(f"Error saving pickle files for nt = {nt}: {e}")

print("\n--- All density data processed and saved as pickle files. ---")

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
mid_slice = mesh // 2

try:
    initial_density_path = f"{output_data_path}/density_{0}.p"
    initial_density = pickle.load(open(initial_density_path, 'rb'))
    im0 = axes[0].imshow(initial_density[:, :, mid_slice], cmap='viridis', origin='lower', animated=True)
    im1 = axes[1].imshow(initial_density[:, mid_slice, :], cmap='viridis', origin='lower', animated=True)
    im2 = axes[2].imshow(initial_density[mid_slice, :, :], cmap='viridis', origin='lower', animated=True)

    axes[0].set_title(f'Density (Z = {mid_slice})')
    axes[0].set_xlabel('X')
    axes[0].set_ylabel('Y')
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    axes[1].set_title(f'Density (Y = {mid_slice})')
    axes[1].set_xlabel('X')
    axes[1].set_ylabel('Z')
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    axes[2].set_title(f'Density (X = {mid_slice})')
    axes[2].set_xlabel('Y')
    axes[2].set_ylabel('Z')
    fig.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

    suptitle_obj = plt.suptitle(f'Density at Time 0 fs, Step 0', fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

except Exception as e:
    print(f"Error initializing animation with first frame data: {e}")
    print("Animation generation aborted.")
    plt.close(fig)
    sys.exit(1)

def animate(i):
    nt_current = i + 1
    density_path = f"{output_data_path}/density_{nt_current}.p"
    
    try:
        current_density = pickle.load(open(density_path, 'rb'))
    except FileNotFoundError:
        print(f"Density file not found for nt = {nt_current}. Skipping frame.")
        return im0, im1, im2, suptitle_obj
    except Exception as e:
        print(f"Error loading density for nt = {nt_current}: {e}. Skipping frame.")
        return im0, im1, im2, suptitle_obj

    im0.set_array(current_density[:, :, mid_slice])
    im1.set_array(current_density[:, mid_slice, :])
    im2.set_array(current_density[mid_slice, :, :])

    suptitle_obj.set_text(f'Density at Time {nt_current * dt} fs, Step {nt_current}')

    return im0, im1, im2, suptitle_obj

anim = animation.FuncAnimation(fig, animate, frames=Nt, interval=200, blit=True)

video_filename = os.path.join(output_data_path, "density.mp4")

print(f"Creating video: {video_filename} using matplotlib.animation...")
try:
    anim.save(video_filename, writer='ffmpeg', dpi=200)
    print(f"Video created successfully at: {video_filename}")
except Exception as e:
    print(f"Error saving animation: {e}")

plt.close(fig)
