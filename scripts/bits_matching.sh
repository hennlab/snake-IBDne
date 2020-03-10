DATASET=$1
BITS_RANGE=`seq 5 5 200` 
ERRHOM_RANGE=`seq 1 2`


echo Pop hap/notHap bits err_hom err_het diff/total match/total > results/other/${DATASET}_GL_IBD_3500kb.txt


for BITS in ${BITS_RANGE[@]}
  do
    for ERRHOM in ${ERRHOM_RANGE[@]}
      do
        grep -i $DATASET results/RoH_param_combs/${BITS}/IBD_${BITS}_${ERRHOM}_hap_segsJoined_QCed.ibd > results/RoH_param_combs/${BITS}/${DATASET}_IBD_${BITS}_${ERRHOM}_hap_segsJoined_QCed.ibd
        echo $DATASET hap $BITS $ERRHOM NA `Rscript scripts/compare_RoH.R results/calc_RoH/${DATASET}_RoH_segsJoined_QCed.hom results/RoH_param_combs/${BITS}/IBD_${BITS}_${ERRHOM}_hap_segsJoined_QCed.ibd 3500` >> results/other/${DATASET}_GL_IBD_3500kb.txt
      done
  done

