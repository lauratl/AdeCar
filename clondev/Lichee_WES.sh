#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user lauratomaslopezslurm@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 10:00:00
#SBATCH --mem 50G
#SBATCH -p thin-shared,cola-corta

#usage Lichee.sh <Config file> <patient> [<depth threshold>]
# Reading config

source ReadConfig.sh $1
PATIENT=$2


# Loading modules

module load jdk/8u181


WORKDIR=${WORKDIR}/wes/clondev/

# Commands

depth=$3
if [ -z "$depth" ]; then
depth=0
fi

LICHEE=/mnt/netapp1/posadalab/APPS/lichee/LICHeE/release/

cd $LICHEE
./lichee -build \
	-i ${WORKDIR}/prepare_clonefinder_input/from_bams/with_python/${PATIENT}.LicheeInput \
	-minVAFPresent 0.05 \
	-maxVAFAbsent 0.025 \
	-maxVAFValid 0.6 \
	-maxClusterDist 0.15 \
	-e 0.1  \
	-minClusterSize 1 \
	-n 0 \
	-o ${WORKDIR}/lichee/${PATIENT}.Lichee \
	-v \
	-color \
	-dot \
	-showTree 4 


echo "FINISHED"
