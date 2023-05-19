conda activate dedalus2
mpiexec -n 1 python3 -m dedalus merge_procs snapshots --cleanup
today="$(date +"%m%d%H%M")"
frames_name="frames_$today"
mv frames "$frames_name"
cp -r "$frames_name" /mnt/d/Data/dedalus_example/FPE
snapshots_name="snapshots_$today"
mv snapshots "$snapshots_name"
cp -r "$snapshots_name" /mnt/d/Data/dedalus_example/FPE