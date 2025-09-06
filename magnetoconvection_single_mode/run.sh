conda activate dedalus2
now="dedalus_$(date +"%Y%m%d_%H%M%S")"
mkdir -p "$now"
cp -r *.py "$now"
cp -r *.h5 "$now"
cd "$now"
python3 magnetoconvection_v2.py 2>&1 | tee dedalus_output_$(date +%Y%m%d_%H%M%S).out
echo $now
echo "finished"