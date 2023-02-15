conda activate dedalus2
mpiexec -n 4 python3 -m dedalus merge_procs snapshots
mpiexec -n 4 python3 plot_slices.py snapshots/*.h5