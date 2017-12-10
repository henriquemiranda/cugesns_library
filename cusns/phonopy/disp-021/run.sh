#!/bin/bash
#PBS -N vasp21
#PBS -q main
#PBS -W group_list=2016038_HPCZENOBE_ULUX
#PBS -l walltime=5:00:00
#PBS -l select=48:mem=2625MB

#pre_run
exec > ${PBS_O_WORKDIR}/${PBS_JOBNAME}_${PBS_JOBID}.log
echo "------------------ Work dir --------------------"
cd ${PBS_O_WORKDIR} && echo ${PBS_O_WORKDIR}
echo "------------------ Job Info --------------------"
echo "jobid : $PBS_JOBID"
echo "jobname : $PBS_JOBNAME"
echo "job type : $PBS_ENVIRONMENT"
echo "submit dir : $PBS_O_WORKDIR"
echo "queue : $PBS_O_QUEUE"
echo "user : $PBS_O_LOGNAME"
echo "threads : $OMP_NUM_THREADS"

#commands
module load mkl/lp64/11.2.1.133 compiler/intel/2015.1.133 intelmpi/5.0.2.044/64
mpirun -np 48 vasp

#pos_run
echo 'done!'

