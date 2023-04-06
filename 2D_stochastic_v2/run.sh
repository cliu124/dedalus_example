conda activate dedalus2
###ln -s /home/cliu124/dedalus_example/2D_stochastic_v2/snapshots_20230326162736/snapshots_s1.h5 restart.h5
mpiexec -n 4 python3 stochastic_forcing_2D_v2.py
mpiexec -n 4 python3 -m dedalus merge_procs snapshots
mpiexec -n 4 python3 -m dedalus merge_procs scalar_data
mpiexec -n 4 python3 plot_slices.py snapshots/*.h5
today="$(date +"%Y%m%d%H%M%S")"
frames_name="frames_$today"
mv frames "$frames_name"
cp -r "$frames_name" /mnt/d/Data/dedalus_example/2D_stochastic_v2
snapshots_name="snapshots_$today"
mv snapshots "$snapshots_name"
cp -r "$snapshots_name" /mnt/d/Data/dedalus_example/2D_stochastic_v2
scalar_data_name="scalar_data_$today"
mv scalar_data "$scalar_data_name"
cp -r "$scalar_data_name" /mnt/d/Data/dedalus_example/2D_stochastic_v2
echo $today
echo "finished"
###rm -rf restart.h5