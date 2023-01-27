mpiexec -n 4 python3 rayleigh_benard.py
mpiexec -n 4 python3 plot_snapshots.py snapshots/*.h5
cp -r frames /mnt/d/Data/dedalus_example/fixed_flux_RBC_v3