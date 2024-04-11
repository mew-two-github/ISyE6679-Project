#!/bin/bash

# Request 2 nodes, 1 GPU, 32 GB memory.
#SBATCH -N2 --mem-per-cpu=16G --gres=gpu:1

# Request 10 minutes of wall-clock time.
#SBATCH -t10

# See https://gatech.service-now.com/technology?id=kb_article_view&sysparm_article=KB0042111 for more information on how to set resources

# Set the output log to be batch.log
#SBATCH -obatch.log

# Change to the submission directory
cd $SLURM_SUBMIT_DIR

# Load the cuda module so we can find the run-time libraries
module load cuda

# Run the code
./qr_solver > myoutput.log
