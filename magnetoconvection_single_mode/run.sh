conda activate dedalus2
now="dedalus_$(date +"%Y%m%d%H%M%S")"
mkdir -p "$now"
cp -r *.py "$now"
cd "$now"
python3 magnetoconvection_v2.py > "$now" 2>&1
echo $now
echo "finished"