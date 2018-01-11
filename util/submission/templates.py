md_template_d144 = """verbosity=0
xcFunctional=PBE
FDtype=4th
[Mesh]
nx=160
ny=80
nz=80
[Domain]
ox=0.
oy=0.
oz=0.
lx=42.4813
ly=21.2406
lz=21.2406
[Potentials]
pseudopotential=pseudo.D_tm_pbe
[Poisson]
solver=@
max_steps_initial=@50
max_steps=@50
reset=@
bcx=periodic
bcy=periodic
bcz=periodic
[Run]
type=MD
[MD]
type=@
num_steps=@
dt=@15.
[XLBOMD]
dissipation=@5
align=@
[Quench]
max_steps=@5
max_steps_tight=@
atol=1.e-@10
num_lin_iterations=3
ortho_freq=100
[SpreadPenalty]
type=@energy
damping=@
target=@1.75
alpha=@0.01
[Orbitals]
initial_type=Gaussian
initial_width=1.5
overallocate_factor=@2.
[ProjectedMatrices]
solver=@short_sighted
[LocalizationRegions]
radius=@8.
auxiliary_radius=@
move_tol=@0.1
[Restart]
input_filename=wave.out
input_level=3
interval=@
"""

md_template_H2O_64 = """verbosity=1
xcFunctional=PBE
FDtype=4th
[Mesh]
nx=128
ny=128
nz=128
[Domain]
ox=0.
oy=0.
oz=0.
lx=23.4884
ly=23.4884
lz=23.4884
[Potentials]
pseudopotential=pseudo.O_ONCV_PBE_SG15
pseudopotential=pseudo.D_ONCV_PBE_SG15
[Poisson]
solver=@
max_steps=@
[Run]
type=MD
[Quench]
max_steps=1000
atol=1.e-@
[MD]
type=@
num_steps=@
dt=10.
print_interval=5
[XLBOMD]
dissipation=@
align=@
[Restart]
input_filename=wave.out
input_level=4
output_level=4
interval=@
"""

quench_template_H2O_64 = """verbosity=1
xcFunctional=PBE
FDtype=4th
[Mesh]
nx=128
ny=128
nz=128
[Domain]
ox=0.
oy=0.
oz=0.
lx=23.4884
ly=23.4884
lz=23.4884
[Potentials]
pseudopotential=pseudo.O_ONCV_PBE_SG15
pseudopotential=pseudo.D_ONCV_PBE_SG15
[Run]
type=QUENCH
[Quench]
max_steps=1000
atol=1.e-8
[Orbitals]
initial_type=Fourier
[Restart]
output_level=4
"""

quench_template_d144 = """verbosity=1
xcFunctional=PBE
FDtype=4th
[Mesh]
nx=160
ny=80
nz=80
[Domain]
ox=0.
oy=0.
oz=0.
lx=42.4813
ly=21.2406
lz=21.2406
[Potentials]
pseudopotential=pseudo.D_tm_pbe
[Poisson]
solver=@
max_steps_initial=@50
max_steps=@50
bcx=periodic
bcy=periodic
bcz=periodic
[Run]
type=QUENCH
[Quench]
max_steps=200
atol=1.e-7
num_lin_iterations=3
ortho_freq=100
[SpreadPenalty]
type=@energy
damping=@
target=@1.75
alpha=@0.01
[Orbitals]
initial_type=Gaussian
initial_width=1.5
[ProjectedMatrices]
solver=@short_sighted
[LocalizationRegions]
radius=@8.
[Restart]
output_type=distributed
"""

H2O_64_params={
    'nodes':                    '32',
    'ntasks':                   '256',
    'omp_num_threads':          8 if omp_num_threads == 4 else omp_num_threads,
    'cores_per_task':           '2',
    'potentials': 'ln -s $maindir/potentials/pseudo.O_ONCV_PBE_SG15\nln -s $maindir/potentials/pseudo.D_ONCV_PBE_SG15',
    'lrs':        '',
    'jobname':                  'H2O_64',
}

d144_params={
    'nodes':                    '8',
    'walltime':                 '01:30:00',
    'ntasks':                   '125',
    'omp_num_threads':          omp_num_threads,
    'cores_per_task':           '1',
    'potentials': 'ln -s $maindir/potentials/pseudo.D_tm_pbe',
    'lrs':        '-l lrs.in',
    'jobname':                  'd144',
}

vulcan_params={
    'queue':                    'psmall',
    'scratch_path':              '/p/lscratchv/mgmolu/dunn27/mgmol/',
    'gres':                     'lscratchv',
    'exe':                      'mgmol-bgq',
}

cab_params={
    'queue':                    'pbatch',
    'scratch_path':              '/p/lscratchd/dunn27/mgmol/',
    'gres':                     'lscratchd',
    'omp_num_threads':          '1',
    'exe':                      'mgmol-pel',
    'walltime':                 '01:30:00',
}

runfile_quench_template="""#!/bin/tcsh
#MSUB -l nodes={nodes},walltime={walltime}
#MSUB -o mgmol.out
#MSUB -q {queue}
#MSUB -A comp
#MSUB -l gres={gres}
#MSUB -N {jobname}

rm -f queued
echo ' ' > running

use boost-nompi-1.55.0
export BOOST_ROOT=/usr/local/tools/boost-nompi-1.55.0
export Boost_NO_SYSTEM_PATHS=ON

setenv OMP_NUM_THREADS {omp_num_threads}

set ntasks = {ntasks}

set maindir = $home/mgmol

set exe     = $maindir/bin/{exe}

set datadir = `pwd`

set scratchdir = {scratch_path}`basename $datadir`
mkdir $scratchdir
cd $scratchdir

echo ' ' > running

set cfg_quench = mgmol_quench.cfg

cp $datadir/$cfg_quench .
cp $datadir/coords.in .
cp $datadir/lrs.in .

{potentials}

#1st run
srun -n $ntasks -c {cores_per_task} $exe -c $cfg_quench -i coords.in {lrs}

#restart
rm -f wave.out
set restart_file=`ls -ld * | awk '/snapshot0/ {{ print $9 }}' | tail -n1`
ln -s -f $restart_file wave.out

rm -f running
echo ' ' > queued
"""

runfile_md_template="""#!/bin/tcsh
#MSUB -l nodes={nodes},walltime={walltime}
#MSUB -o mgmol.out
#MSUB -q {queue}
#MSUB -A comp
#MSUB -l gres={gres}
#MSUB -N {jobname}

rm -f queued
echo ' ' > running

use boost-nompi-1.55.0
export BOOST_ROOT=/usr/local/tools/boost-nompi-1.55.0
export Boost_NO_SYSTEM_PATHS=ON

setenv OMP_NUM_THREADS {omp_num_threads}

set ntasks = {ntasks}

set maindir = $home/mgmol

set exe     = $maindir/bin/{exe}

set datadir = `pwd`

set scratchdir = {scratch_path}`basename $datadir`
mkdir $scratchdir
cd $scratchdir

echo ' ' > running

set cfg_md = mgmol_md.cfg

cp $datadir/$cfg_md .

#restart
rm -f wave.out
set restart_file=`ls -ld * | awk '/snapshot0/ {{ print $9 }}' | tail -n1`
ln -s -f $restart_file wave.out

#MD run
srun -n $ntasks -c {cores_per_task} $exe -c $cfg_md

#restart
rm -f wave.out
set restart_file=`ls -ld * | awk '/snapshot0/ {{ print $9 }}' | tail -n1`
ln -s -f $restart_file wave.out

rm -f running
echo ' ' > queued
"""
