# Usage qsub sim_jobarray.sh

## This script runs the bis5-long_numSpine-5_diffDistribution.py script.
## Sarthak Sharma, Dec 13, 2016

echo $MODELS
#$ -S /bin/bash
#$ -cwd 
#$ -V
#$ -o /home/sarthak/out/
#$ -e /home/sarthak/err/
#$ -r y #Restarts the job in case of failiures
## Returning error code 99/100 to SGE from a program would schedule it for restart or mark it as E respectively

#DATA_OUT=/sarthak/diffInit/$JOB_NAME

echo "Starting job on: " `hostname`
echo "Starting date and time: " `date`
echo "Total cores demanded: " $NSLOTS
echo "Job name given: " $JOB_NAME
echo "Job ID: " $JOB_ID
echo "Task ID: " $SGE_TASK_ID

echo "Running Series simulation now."
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
export PYTHONPATH=$HOME/work/moose-core/_build/python

## Python main simulation file using MOOSE 
/usr/local/bin/python2.7 /home/sarthak/work/fhn/diffSpaces/codes/fhn_long_uniformSpacing.py # Iterations can also be an argument here. 
if [ "$?" = "0" ]; then
    echo "Simulation finished successfully"
else
    echo "Error in simulation !" 1>&2 # Dumping in stderr in this case.
fi
