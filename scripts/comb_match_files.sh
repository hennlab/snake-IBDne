BITS_RANGE=`seq 5 5 200` # GERMLINE2 IBD segment seed
ERRHOM_RANGE=`seq 1 2`
CHRS=`seq 1 22`
DATASET=$1

for BITS in ${BITS_RANGE[@]}
  do 
    for ERRHOM in ${ERRHOM_RANGE[@]}
      do 
        for chr in ${CHRS[@]}
          do
            cat results/RoH_param_combs/${BITS}/${DATASET}.chr${chr}.phased_GERMLINE2_RoH_${BITS}_${ERRHOM}_hap.match >> results/RoH_param_combs/${BITS}/IBD_${BITS}_${ERRHOM}_hap.match
          done 
        done 
    done

