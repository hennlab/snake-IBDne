
file = $1 
char = sed 

if [ $chr -eq 1 ]
then
cat Himba_MEGA_10samples_chr${chr}.phased.vcf > Himba_MEGA_10samples_allchr.phased.vcf
else 
sed '/^#/ d' < Himba_MEGA_10samples_chr${chr}.phased.vcf >> Himba_MEGA_10samples_allchr.phased.vcf
fi
done
