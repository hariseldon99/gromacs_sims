#!/bin/bash
#SBATCH --job-name=rsync_job
#SBATCH --ntasks=1
#SBATCH --qos=elevated
#SBATCH --time=90:00:00

# Set the source and destination directories
SRC_DIR="/home/admin/gitrepos/gromacs_sims/MMP2_Rolipram/Rolipram/"
DEST_DIR="/home/admin/mnt/gdrive/data/GROMACS/Rolipram/"

# Set the rsync options
RSYNC_OPTS="-avz --update"

# Set the sleep interval in seconds
SLEEP_INTERVAL=300

while true; do
    # Run the rsync command
    rsync $RSYNC_OPTS $SRC_DIR $DEST_DIR

    # Sleep for the specified interval
    sleep $SLEEP_INTERVAL
done
