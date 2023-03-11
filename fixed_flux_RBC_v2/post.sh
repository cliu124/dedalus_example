conda activate dedalus2
mpiexec -n 4 python3 -m dedalus merge_procs snapshots
mpiexec -n 4 python3 plot_slices.py snapshots/*.h5
name="frames_$(date +"%Y%m%d%H%M%S")"
mv frames "$name"
cp -r "$name" /mnt/d/Data/dedalus_example/fixed_flux_RBC_v2