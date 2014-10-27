
cat << EOF > q_wat_dual
cd \$PBS_O_WORKDIR
mpirun -np 8 protoms3 run_free.cmd
EOF
qsub -l walltime=2:00:00 -l nodes=1:ppn=8 q_wat_dual
