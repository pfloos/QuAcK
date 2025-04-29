#!/bin/bash

# Starting, ending, and step size
eta_start=0.00005
eta_end=0.007
eta_step=0.00005

# Job script name
job_script="$QUACK_ROOT/submission_scripts/submit_script_eta_scan.sh"

# Number format for eta
format="%.6f"

# Function to compare floating point values
float_ge() { awk -v n1="$1" -v n2="$2" 'BEGIN {exit !(n1 >= n2)}'; }

eta=$eta_start
while float_ge "$eta_end" "$eta"; do

    echo "Running with eta = $eta"

    # Write eta to eta_opt.dat (adjust if more complex format is needed)
    echo "$eta" > "$QUACK_ROOT/input/eta_opt.dat"

    # Submit job and capture job ID
    job_id=$(sbatch --parsable "$job_script")
    echo "Submitted job $job_id"

    # Wait for job to finish
    while squeue -j "$job_id" > /dev/null 2>&1 && squeue -j "$job_id" | grep -q "$job_id"; do
        sleep 10  # Check every 10 seconds
    done

    # Increase eta
    eta=$(awk -v a="$eta" -v s="$eta_step" 'BEGIN {printf "%.6f", a + s}')
done

