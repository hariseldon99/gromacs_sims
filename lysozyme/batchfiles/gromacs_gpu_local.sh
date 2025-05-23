#!/bin/bash

export num_proc=64

#Use openmp. Set to false to switch to MPI
export USE_OPENMP=true

#Start time
start=`date +%s.%N`

#Load basic OHPC tools
module load ohpc
#Load cuda
module load cuda
#Load singularity
module load singularity
export SIFPATH=$HOME/images
export SIFIMG=gromacs_2022.3.sif

echo "Starting"
echo '---------------------------------------------'

export MDNAME=md_0_1
export GROFILE=npt.gro
export CPTFILE=npt.cpt
export MDPFILE=md.mdp
export TOPOL_FILE=topol.top

#Preprocessing 
LD_LIBRARY_PATH="" singularity run --nv -B ${PWD}:/host_pwd --pwd /host_pwd $SIFPATH/$SIFIMG gmx grompp -f $MDPFILE -c $GROFILE -t $CPTFILE -p $TOPOL_FILE -o ${MDNAME}.tpr

if [ "$USE_OPENMP" = true ]
then
    export OMP_NUM_THREADS=$num_proc
    export MPI_NUM_PROCS=1
else
    export OMP_NUM_THREADS=1
    export MPI_NUM_PROCS=$num_proc
fi

#Actual MD Dynamics: 
LD_LIBRARY_PATH="" singularity run --nv -B ${PWD}:/host_pwd --pwd /host_pwd $SIFPATH/$SIFIMG gmx mdrun -ntmpi $MPI_NUM_PROCS -nb gpu -pin on -v -deffnm -ntomp $OMP_NUM_THREADS $MDNAME


#End time
end=`date +%s.%N`

echo "OMP_NUM_THREADS= "$OMP_NUM_THREADS", MPI_NUM_PROCS= "$MPI_NUM_PROCS
export RUNTIME=$( echo "$end - $start" | bc -l )
echo '---------------------------------------------'
echo "Runtime: "$RUNTIME" sec"
echo '---------------------------------------------'
#----------------------------------------------------------
# Communicate job status to a telegram bot
#----------------------------------------------------------
# <-Create a telegram bot 
# <-Get TOKEN, CHATID from botfather
# <-See https://www.cytron.io/tutorial/how-to-create-a-telegram-bot-get-the-api-key-and-chat-id
# <-Put them into two environment variables TOKEN and CHATID 
# <-Store them in a config file and source it like below
#----------------------------------------------------------

source ${HOME}/.config/telegram/telegram.conf

LOGIN_NODE=kuhpchn
SSHBIN=/usr/bin/ssh
URL="https://api.telegram.org/bot${TOKEN}/sendMessage"
# Generate the telegram message  text
TEXT="${bell} Local Job launched on ${start} exiting @ ${HOSTNAME}:${PWD}"
CMD='curl -s --max-time 10 --retry 5 --retry-delay 2 --retry-max-time 10  -d '\""chat_id=${CHATID}&text=${TEXT}&disable_web_page_preview=true&parse_mode=markdown"\"" ${URL}"
$SSHBIN $LOGIN_NODE $CMD
