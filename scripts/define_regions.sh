FILE=$1
DATASET=$2
BITS=$(echo $FILE | cut -d '_' -f5)
ERRHOM=$(echo $FILE | cut -d '_' -f6)

N=`wc -l $FILE | cut -f1 -d" "`
if [ $N -gt 0 ]
  then
    Rscript scripts/get_outlier_regions.R results/RoH_param_combs/${BITS}/GERMLINE2_IBDdepth_${BITS}_${ERRHOM}_hap_IBDoutliers 3e6 2.5e5 > results/RoH_param_combs/${BITS}/GERMLINE2_IBDdepth_${BITS}_${ERRHOM}_hap_IBDoutliers_regions
    cat regions/exclude_regions_hg19.txt results/RoH_param_combs/${BITS}/GERMLINE2_IBDdepth_${BITS}_${ERRHOM}_hap_IBDoutliers_regions  results/calc_RoH/${DATASET}_lowDensityRegions.txt > results/RoH_param_combs/${BITS}/exclude_regions_hg19_forIBDNe_${BITS}_${ERRHOM}_hap.txt
  else
    cat regions/exclude_regions_hg19.txt results/calc_RoH/${DATASET}_lowDensityRegions.txt > results/RoH_param_combs/${BITS}/exclude_regions_hg19_forIBDNe_${BITS}_${ERRHOM}_hap.txt
  fi

