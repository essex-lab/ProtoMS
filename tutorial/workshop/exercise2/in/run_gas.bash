
cat << EOF > q_gas
cd \$PBS_O_WORKDIR
mpirun -np $1 protoms3 run_comb_gas.cmd
EOF
qsub -l walltime=0:30:00 -l nodes=1:ppn=$1 q_gas
