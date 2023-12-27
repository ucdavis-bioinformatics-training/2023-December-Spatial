#!/bin/bash
#
#SBATCH --time=0-2:00:00  # days-hours:minutes:seconds
#SBATCH --job-name=spaceranger  # Job name
#SBATCH --array=1-2 # samples to run (reads lines in design.tsv)
#SBATCH --nodes=1
#SBATCH --ntasks=8 # Number of cores
#SBATCH --mem=16G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --account=workshop
#SBATCH --output=slurmout/spaceranger-count-%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/spaceranger-count-%A_%a.err # File to which STDERR will be written
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=user@domain.com # Email to which notifications will be sent

workshop=/share/workshop/Spatial_Workshop
export PATH=${workshop}/Software/spaceranger-2.1.1:$PATH
transcriptome=${workshop}/Reference/refdata-gex-mm10-2020-A
project=${workshop}/$USER
fastqs=${project}/00-RawData
outdir=${project}/01-Spaceranger
metadata=`sed "${SLURM_ARRAY_TASK_ID}q;d" ${project}/design.tsv`
sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" ${project}/design.tsv | cut -f1`
image=`sed "${SLURM_ARRAY_TASK_ID}q;d" ${project}/design.tsv | cut -f2`
serial=`sed "${SLURM_ARRAY_TASK_ID}q;d" ${project}/design.tsv | cut -f3`
area=`sed "${SLURM_ARRAY_TASK_ID}q;d" ${project}/design.tsv | cut -f4`

start=`date +%s`

mkdir -p ${outdir}

echo $HOSTNAME > ${outdir}/${sample}.log.out
ulimit -n 16000
ulimit -a >> ${outdir}/${sample}.log.out

echo "Sample ${sample}" >> ${outdir}/${sample}.log.out
echo "Allocated threads: " ${SLURM_NTASKS} >> ${outdir}/${sample}.log.out
echo "Allocated memory: " ${SLURM_MEM_PER_NODE} >> ${outdir}/${sample}.log.out

call="spaceranger count \
  --id=${sample} \
  --fastqs=${fastqs} \
  --sample=${sample} \
  --transcriptome=${transcriptome} \
  --image=${fastqs}/images/${image} \
  --slide=$serial \
  --area=$area \
  --reorient-images=true \
  --nosecondary \
  --localcores=${SLURM_NTASKS} \
  --localmem=${SLURM_MEM_PER_NODE} \
>> ${outdir}/${sample}.log.out \
2> ${outdir}/${sample}.log.err"

echo $call >> ${outdir}/${sample}.log.out

cd $outdir

eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime >> ${outdir}/${sample}.log.out
