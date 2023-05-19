conda activate dedalus2
mkdir analysis
###ln -s /home/cliu124/dedalus_example/2D_stochastic_v2/analysis_20230326162736/analysis_s1.h5 restart.h5
mpiexec -n 1 python3 FPE.py
mpiexec -n 1 python3 -m dedalus merge_procs analysis --cleanup
mpiexec -n 1 python3 -m dedalus merge_procs scalar_data --cleanup
today="$(date +"%Y%m%d%H%M%S")"
frames_name="frames_$today"
mv frames "$frames_name"
cp -r "$frames_name" /mnt/d/Data/dedalus_example/FPE
analysis_name="analysis_$today"
mv analysis "$analysis_name"
cp -r "$analysis_name" /mnt/d/Data/dedalus_example/FPE
scalar_data_name="scalar_data_$today"
mv scalar_data "$scalar_data_name"
cp -r "$scalar_data_name" /mnt/d/Data/dedalus_example/FPE
echo $today
echo "finished"
###rm -rf restart.h5