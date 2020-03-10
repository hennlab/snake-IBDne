DATASET=$1
BITS_RANGE=`seq 5 5 200`
ERRHOM_RANGE=`seq 1 2`

echo 'BITS ERRHOM ERRHET NRMSE' > results/other/${DATASET}_NormRMSE.txt
GENOME_LENGTH=`Rscript scripts/calc_genome_length.R results/other/${DATASET}_HWE-MAFfiltered.bim`

for BITS in ${BITS_RANGE[@]}
  do
    for ERRHOM in ${ERRHOM_RANGE[@]}
      do
        echo $BITS $ERRHOM NA `Rscript scripts/compare_IBD_kinship.R results/RoH_param_combs/${BITS}/IBD_${BITS}_${ERRHOM}_hap_segsJoined_QCed.ibd results/other/${DATASET}.genome $GENOME_LENGTH` >> results/other/${DATASET}_NormRMSE.txt
      done
  done

