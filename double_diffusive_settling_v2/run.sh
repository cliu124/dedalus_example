conda activate dedalus2
mpiexec -n 4 python3 double_diffusive_settling.py
mpiexec -n 4 python3 -m dedalus merge_procs snapshots
mpiexec -n 4 python3 plot_slices.py snapshots/*.h5
today="$(date +"%Y%m%d%H%M%S")"
frames_name="frames_$today"
mv frames "$frames_name"
cp -r "$frames_name" /mnt/d/Data/dedalus_example/double_diffusive_settling_v2
snapshots_name="snapshots_$today"
cp double_diffusive_settling.py snapshots/
cp plot_slices.py snapshots/
mv snapshots "$snapshots_name"
cp -r "$snapshots_name" /mnt/d/Data/dedalus_example/double_diffusive_settling_v2

