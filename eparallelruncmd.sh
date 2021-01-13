#!/bin/bash

set -e
set -x
module load nci-parallel/1.0.0

/apps/openmpi/4.0.2/bin/mpirun -map-by ppr:$tasks_per_node:node:PE=$cpupertask /apps/nci-parallel/1.0.0/bin/nci-parallel --poll 0 --shell /bin/bash --input-file $commandsfile
