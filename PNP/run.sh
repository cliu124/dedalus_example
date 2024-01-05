conda activate dedalus2
mkdir analysis
###ln -s /home/cliu124/dedalus_example/2D_stochastic_v2/analysis_20230326162736/analysis_s1.h5 restart.h5
mpiexec -n 1 python3 PNP.py
mpiexec -n 1 python3 -m dedalus merge_procs analysis --cleanup
today="$(date +"%Y%m%d%H%M%S")"
frames_name="frames_$today"
mv frames "$frames_name"
cp -r "$frames_name" /mnt/d/Data/dedalus_example/PNP
analysis_name="analysis_$today"
mv analysis "$analysis_name"
cp -r "$analysis_name" /mnt/d/Data/dedalus_example/PNP
scalar_data_name="scalar_data_$today"
echo $today
echo "finished"