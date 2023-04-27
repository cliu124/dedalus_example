conda activate dedalus2
mpiexec -n 4 python3 -m dedalus merge_procs analysis --cleanup
mpiexec -n 4 python3 plot_slices.py analysis/*.h5
today="$(date +"%Y%m%d%H%M%S")"
frames_name="frames_$today"
mv frames "$frames_name"
cp -r "$frames_name" /mnt/d/Data/dedalus_example/double_diffusive_settling_v2
cp double_diffusive_settling.py analysis/
cp plot_slices.py analysis/
analysis_name="analysis_$today"
mv analysis "$analysis_name"
cp -r "$analysis_name" /mnt/d/Data/dedalus_example/double_diffusive_settling_v2