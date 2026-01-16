#BSUB -J "RCTD_drug_resist[1-4]"
#BSUB -o %J_%I.out
#BSUB -e %J_%I.err
#BSUB -W 70:00
#BSUB -n 20  
#BSUB -R "span[hosts=1]"
#BSUB -M 500000
#BSUB -R "rusage[mem=500000]"

source /users/zhozp5/miniconda3/etc/profile.d/conda.sh 
conda activate hdWGCNA

SAMPLES=("4h_em" "18h_em" "32h_em" "44h_em")
SAMPLE_ID=${SAMPLES[$(($LSB_JOBINDEX-1))]}

# Run R script with sample name as argument
Rscript /data/xu_lab_projectsx/yuanzhou/drug_resist/scripts/02_rctd/rctd.R $SAMPLE_ID
