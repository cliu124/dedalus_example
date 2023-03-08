conda activate dedalus2
mpiexec -n 4 python3 double_diffusive_settling.py
mpiexec -n 4 python3 -m dedalus merge_procs snapshots
mpiexec -n 4 python3 plot_slices.py snapshots/*.h5
name="frames_$(date +"%Y%m%d%H%M%S")"
mv frames "$name"
mv snapshots "snapshots_$(date +"%Y%m%d%H%M%S")"
cp -r "$name" /mnt/d/Data/dedalus_example/double_diffusive_settling_v2