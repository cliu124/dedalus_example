conda activate dedalus2
###ln -s /home/cliu124/dedalus_example/horizontal_convection_v2/analysis_20230326162736/analysis_s1.h5 restart.h5
mpiexec -n 4 python3 horizontal_convection.py
mpiexec -n 4 python3 -m dedalus merge_procs analysis --cleanup
mpiexec -n 4 python3 plot_slices.py analysis/*.h5
today="$(date +"%Y%m%d%H%M%S")"
frames_name="frames_$today"
mv frames "$frames_name"
cp -r "$frames_name" /mnt/d/Data/dedalus_example/horizontal_convection_v2
analysis_name="analysis_$today"
mv analysis "$analysis_name"
cp -r "$analysis_name" /mnt/d/Data/dedalus_example/horizontal_convection_v2
echo $today
echo "finished"