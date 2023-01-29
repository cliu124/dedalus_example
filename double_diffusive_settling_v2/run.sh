conda activate dedalus2
mpiexec -n 4 python3 double_diffusive_settling.py
mpiexec -n 4 python3 -m dedalus merge_procs snapshots
mpiexec -n 4 python3 plot_slices.py snapshots/*.h5
cp -r frames /mnt/d/Data/dedalus_example/double_diffusive_settling_v2