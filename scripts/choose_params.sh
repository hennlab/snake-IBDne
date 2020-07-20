# Choose best parameters for IBDne

# usage:
# bash choose_params.sh POPS_GL_IBD_3500kb.txt POPS_NormRMSE.txt POPS


IBD_3500=$1
NormRMSE=$2
P=$3

paste -d" " <(awk '{print $0}' $IBD_3500) <(awk '{print $4}' $NormRMSE) > temp
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, (1-$6 + $7 + 1-$8)/3}' temp > temp2
sed -r '1s/\S+//9' temp2 > temp
sed -i '1s/$/score/' temp
# Look only at haploid
mv temp ${P}_GL_IBD_3500kb_onlyHap.txt; rm temp2

head -n 1 ${P}_GL_IBD_3500kb_onlyHap.txt > header.txt
sed '1d' ${P}_GL_IBD_3500kb_onlyHap.txt > data.txt

sort -nk9 data.txt > sorted_data.txt ; rm data.txt
cat header.txt sorted_data.txt > ${P}_GL_IBD_3500kb_onlyHap.txt ; rm header.txt sorted_data.txt

cat ${P}_GL_IBD_3500kb_onlyHap.txt

PARAMS=`tail -n 1 ${P}_GL_IBD_3500kb_onlyHap.txt | cut -f 2-9 -d" "`
PARAMS=(${PARAMS//,/ })
echo $PARAMS
BITS=${PARAMS[1]}
ERRHOM=${PARAMS[2]}
PERC_DIFF=${PARAMS[4]}
PERC_MATCH=${PARAMS[5]}
KIN_RESIDS=${PARAMS[6]}
SCORE=${PARAMS[7]}

echo -e $P GERMLINE2 parameters"\n"Type: ${TYPE}"\n"Bits: ${BITS}"\n"Err_hom: ${ERRHOM}"\n"Percent difference: ${PERC_DIFF}"\n"Percent match: ${PERC_MATCH}"\n"Kin residuals: ${KIN_RESIDS}"\n"Score: ${SCORE} > results/other/${P}_GL_parameters.txt
# awk -v Pop="$P" -v IGNORECASE=1 '$1 ~ Pop && $3 ~ Pop {print}' RoH_param_combs/${BITS}/IBD_${BITS}_${ERRHOM}_hap_segsJoined_QCed.ibd | cut -f1-7,9 > ${P}_only_QCed.ibd


cat results/RoH_param_combs/${BITS}/IBD_${BITS}_${ERRHOM}_hap_segsJoined_QCed.ibd > results/other/${P}_only_QCed.ibd
