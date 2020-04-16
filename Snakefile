"""
Title: snake-IBDne
Description: Custom pipeline to infer historical effective population sizes using IBDne
Authors: Shyamalika Gopalan, Mira Mastoras
"""

# -------
# SET-UP
# -------

import os

# Config file:
DATASET = config['dataset'] # name of dataset  ie Himba_MEGA_10samples
GENMAP_CHR=config['gmap_chr_dir'] # directory where genetic map files for individual chromosomes are located
PHASED=config['phased']
REF=config['ref']
# Variables
CHR = [i for i in range(1,23)]
BITS = [i for i in range(5,201,5)]
ERRHOM = [1,2]
IBDne_THREADs = 10 # specify number of threads to use for running the IBDne program

# ROH parameters
MAF=0.01
GENO=0.05
HWE=0.001
RoH_LEN=500
N_SNPS=30
MISSING=2
HET=1
cM=4
GMAX=100

# Define input files rules which combine
def chrom_combine_inputs(wildcards):
    files = expand("prepare/chroms_phased_vcf/{dataset}.chr{chrnum}.phased.vcf", chrnum=CHR, dataset=DATASET)
    return files

def concordance_input(wildcards):
    files = expand("results/RoH_param_combs/{bits}/IBD_{bits}_{err}_hap_segsJoined_QCed.ibd", bits=BITS, err=ERRHOM)
    return files

def comb_match_inputs(wildcards):
    files = expand("results/RoH_param_combs/{bits}/{dataset}.chr{chrnum}.phased_GERMLINE2_RoH_{bits}_{err}_hap.match", bits=BITS, chrnum=CHR, err=ERRHOM, dataset=DATASET)
    return files
def bits_matching_inputs(wildcards):
    files = expand("results/RoH_param_combs/{bits}/IBD_{bits}_{err}_hap_segsJoined_QCed.ibd", bits=BITS, err=ERRHOM)
    return files
# -------
# TARGETS
# -------

rule all:
  input:
    expand("results/IBDne_output/{dataset}_GERMLINE2_IBDNe_minibd-{cM}.log", dataset=DATASET, cM=cM),
    expand("results/IBDne_output/{dataset}_GERMLINE2_IBDNe_minibd-{cM}.pair.excl", dataset=DATASET, cM=cM),
    expand("results/IBDne_output/{dataset}_GERMLINE2_IBDNe_minibd-{cM}.region.excl", dataset=DATASET, cM=cM),
    expand("results/IBDne_output/{dataset}_GERMLINE2_IBDNe_minibd-{cM}.ne", dataset=DATASET, cM=cM),
    expand("results/IBDne_output/{dataset}_GERMLINE2_IBDNe_minibd-{cM}.boot", dataset=DATASET, cM=cM),
    expand("results/calc_RoH/{dataset}_RoHdepth_outliers", dataset=DATASET),
    expand("results/calc_RoH/{dataset}_RoHdepth_outlierRegions", dataset=DATASET),
    expand("results/other/{dataset}.genome", dataset=DATASET),
    expand("results/other/{dataset}_GL_parameters.txt", dataset=DATASET),
    expand("results/other/{dataset}.allchr.phased_lowDensityRegions.txt", dataset=DATASET),

# ------------------------
# PREPARE FILES FOR IBDne
# ------------------------

if PHASED=='FALSE':
  rule break_chrom:
    input:
      "data/{dataset}.bim"
    output:
      multiext("prepare/chroms/{dataset}.chr{chrnum}", ".bed", ".bim", ".fam")
    benchmark:
      "benchmarks/{dataset}/{chrnum}/break_chrom.txt"
    params:
      "prepare/chroms/{dataset}.chr{chrnum}"
    shell:
      """
      plink --bfile data/{wildcards.dataset} --chr {wildcards.chrnum} --make-bed --out {params}
      """

  rule phasing:
    input:
      multiext("prepare/chroms/{dataset}.chr{chrnum}", ".bed", ".bim", ".fam" )
    output:
      multiext("data/{dataset}.chr{chrnum}.phased", ".haps", ".sample"),
      log = "shapeit_log/{dataset}_{chrnum}.phased.log"
    params:
      inputmap = GENMAP_CHR+"chr{chrnum}.gmap.txt",
      in_pre="prepare/chroms/{dataset}.chr{chrnum}",
      out_pre="data/{dataset}.chr{chrnum}.phased"
    benchmark:
      "benchmarks/{dataset}/{chrnum}/phasing.txt"
    shell:
      """
      shapeit --input-bed {params.in_pre} --input-map {params.inputmap} --input-ref {REF}chr{chrnum}.hap.gz {REF}chr{chrnum}.legend {REF}.sample --output-max {params.out_pre} --output-log {output.log} --duohmm -W 5
      """

rule convert_germline:
  input:
    multiext("data/{dataset}.chr{chrnum}.phased", ".haps", ".sample")
  output:
    multiext("data/{dataset}.chr{chrnum}.phased", ".ped", ".map")
  params:
    prefix = "data/{dataset}.chr{chrnum}.phased",
    gmap = GENMAP_CHR+"chr{chrnum}.gmap.txt"
  benchmark:
    "benchmarks/{dataset}/{chrnum}/convert_germline.txt"
  shell:
    """
    python scripts/shapeit_to_germline.py {params.prefix} {params.gmap}
    """

rule make_vcf:
  input:
    multiext("data/{dataset}.chr{chrnum}.phased", ".haps", ".sample")
  output:
    vcf = "prepare/chroms_phased_vcf/{dataset}.chr{chrnum}.phased.vcf",
    log = "shapeit_log/{dataset}_{chrnum}.vcf.log"
  benchmark:
    "benchmarks/{dataset}/{chrnum}/make_vcf.txt"
  params:
    in_pre="data/{dataset}.chr{chrnum}.phased",
  shell:
    """
    shapeit -convert --input-haps {params.in_pre} --output-vcf {output.vcf} --output-log {output.log}
    """

rule combine_chrom:
  input:
    chrom_combine_inputs
  output:
    "prepare/all_chroms/{dataset}.allchr.phased.vcf"
  benchmark:
    "benchmarks/{dataset}/combine_chrom.txt"
  shell:
    """
    sed -n -e '/^#/ p' prepare/chroms_phased_vcf/{wildcards.dataset}.chr1.phased.vcf > {output} #initialize allchr.vcf with headers
    for i in {input}; do sed '/^#/d' $i >> {output} ; done
    """

rule make_allchr_map:
  input:
    "prepare/all_chroms/{dataset}.allchr.phased.vcf"
  output:
    "prepare/all_chroms/{dataset}.allchr.phased.map"
  params:
    "prepare/all_chroms/{dataset}.allchr.phased"
  shell:
    """
    plink --vcf {input} --recode --out {params}
    """

rule add_cm:
  input:
    "prepare/all_chroms/{dataset}.allchr.phased.map",
  output:
    "prepare/all_chroms/{dataset}.allchr.phased.cm.map"
  params:
    gmap = GENMAP_CHR+"chr@.gmap.txt",
    in_pre = "prepare/all_chroms/{dataset}.allchr.phased",
    out_pre = "prepare/all_chroms/{dataset}.allchr.phased.cm"
  shell:
    """
    plink --file {params.in_pre} --cm-map {params.gmap} --recode --out {params.out_pre}
    """

rule make_bim:
  input:
    "prepare/all_chroms/{dataset}.allchr.phased.vcf"
  output:
    "prepare/all_chroms/{dataset}.allchr.phased.bim"
  benchmark:
    "benchmarks/{dataset}/make_bim.txt"
  params:
    out_pre="prepare/all_chroms/{dataset}.allchr.phased"
  shell:
    """
    plink --vcf {input} --make-just-bim -out {params.out_pre}
    """

rule low_Density:
  input:
    "prepare/all_chroms/{dataset}.allchr.phased.bim"
  output:
    "results/other/{dataset}.allchr.phased_lowDensityRegions.txt"
  params:
    "results/other/{dataset}.allchr.phased"
  benchmark:
    "benchmarks/{dataset}/low_Density.txt"
  shell:
    """
    Rscript scripts/find_lowDensityRegs.R {input} {params} 1e6 0.05
    """

# ---------------------------------
# SYSTEMATICALLY SAMPLE PARAMETERS
# ---------------------------------

rule flip:
  input:
    "data/{dataset}.chr{chrnum}.phased.ped"
  output:
    "data/{dataset}.chr{chrnum}.phased.flip.ped"
  shell:
    """
    gawk ' {{ t = $1; $1 = $2; $2 = t; print; }} ' {input} > {output}
    """

rule gen_Match:
  input:
    map = "data/{dataset}.chr{chrnum}.phased.map",
    ped = "data/{dataset}.chr{chrnum}.phased.flip.ped"
  output:
    "results/RoH_param_combs/{bits}/{dataset}.chr{chrnum}.phased_GERMLINE2_RoH_{bits}_{err}_hap.match"
  params:
    prefix = "results/RoH_param_combs/{bits}/{dataset}.chr{chrnum}.phased_GERMLINE2_RoH_{bits}_{err}_hap",
  shell:
    """
    GERMLINE2 -mapfile {input.map} -pedfile {input.ped} -outfile {params.prefix} -err_hom {wildcards.err} -bits {wildcards.bits} -w_extend -haploid
    """

rule comb_match_files:
  input:
    comb_match_inputs
  output:
    expand("results/RoH_param_combs/{bits}/IBD_{bits}_{err}_hap.match", bits=BITS, err=ERRHOM)
  shell:
    """
    scripts/comb_match_files.sh {DATASET}
    """

rule prepare_ibd:
  input:
    "results/RoH_param_combs/{bits}/IBD_{bits}_{err}_hap.match"
  output:
    "results/RoH_param_combs/{bits}/IBD_{bits}_{err}_hap.ibd"
  shell:
    """
    awk '{{print $1 , "0" , $2 , "0" , $3 , $4 , $5, "20"}}' {input} > {output}
    """

rule repair_ibd:
  input:
    ibd = "results/RoH_param_combs/{bits}/IBD_{bits}_{err}_hap.ibd",
    vcf = expand("prepare/all_chroms/{dataset}.allchr.phased.vcf", dataset=DATASET),
    map = expand("prepare/all_chroms/{dataset}.allchr.phased.cm.map", dataset=DATASET)
  output:
    "results/RoH_param_combs/{bits}/IBD_{bits}_{err}_hap_segsJoined.ibd"
  shell:
    """
    cat {input.ibd} | java -jar progs/merge-ibd-segments.16May19.ad5.jar {input.vcf} {input.map} 0.6 1 | awk '{{print $1, $2, $3, $4, $5, $6, $7, $9}}' > {output}
    """

rule get_seg_depth:
  input:
    ibd = "results/RoH_param_combs/{bits}/IBD_{bits}_{err}_hap_segsJoined.ibd",
    bim = expand("prepare/all_chroms/{dataset}.allchr.phased.bim", dataset=DATASET)
  output:
    "results/RoH_param_combs/{bits}/GERMLINE2_IBDdepth_{bits}_{err}_hap"
  benchmark: "benchmarks/get_seg_depth_{bits}_{err}.log"
  shell:
    """
    Rscript scripts/get_seg_depthv3.R {input.ibd} {input.bim} > {output}
    """

rule get_outlier_SNPs:
  input:
    "results/RoH_param_combs/{bits}/GERMLINE2_IBDdepth_{bits}_{err}_hap"
  output:
    "results/RoH_param_combs/{bits}/GERMLINE2_IBDdepth_{bits}_{err}_hap_IBDoutliers"
  shell:
    """
    Rscript scripts/get_outlier_SNPs.R {input} 3 > {output}
    """

rule define_regions:
  input:
    "results/RoH_param_combs/{bits}/GERMLINE2_IBDdepth_{bits}_{err}_hap_IBDoutliers"
  output:
    "results/RoH_param_combs/{bits}/exclude_regions_hg19_forIBDNe_{bits}_{err}_hap.txt"
  shell:
    """
    bash scripts/define_regions.sh {input} {DATASET}
    """

rule filter_IBD:
  input:
    ibd = "results/RoH_param_combs/{bits}/IBD_{bits}_{err}_hap_segsJoined.ibd",
    txt = "results/RoH_param_combs/{bits}/exclude_regions_hg19_forIBDNe_{bits}_{err}_hap.txt"
  output:
    "results/RoH_param_combs/{bits}/IBD_{bits}_{err}_hap_segsJoined_QCed.ibd"
  shell:
    """
    Rscript scripts/filter_RoH-IBD_segs.R {input.ibd} {input.txt} 0.85 > {output}
    """

# ----------------------------
# CALCULATE ROH (in parallel)
# ----------------------------

rule calc_RoH:
  input:
    multiext("data/{dataset}", ".bed", ".bim", ".fam")
  output:
    "results/calc_RoH/{dataset}_RoH.hom"
  params:
    in_pre = "data/{dataset}",
    out_pre = "results/calc_RoH/{dataset}_RoH"
  benchmark:
    "benchmarks/{dataset}/calc_RoH.txt"
  shell:
    """
    plink --bfile {params.in_pre} --homozyg --homozyg-snp {N_SNPS} --homozyg-window-missing {MISSING} --homozyg-window-het {HET} --homozyg-kb {RoH_LEN} --out {params.out_pre}
    """

rule join_roh:
  input:
    bim = "data/{dataset}.bim",
    hom = "results/calc_RoH/{dataset}_RoH.hom"
  output:
    "results/calc_RoH/{dataset}_RoH_segsJoined.hom"
  benchmark:
    "benchmarks/{dataset}/join_roh.txt"
  shell:
    """
    Rscript scripts/join_RoHsegs.R {input.bim} {input.hom} 1e6 0.05 0.8
    """

# ----------------------
# Filter RoH segments
# ---------------------

rule roh_lowdens:
  input:
    "data/{dataset}.bim"
  output:
    ldr = "results/calc_RoH/{dataset}_lowDensityRegions.txt",
    reg = "results/calc_RoH/{dataset}_temp_regions.txt"
  benchmark:
    "benchmarks/{dataset}/roh_lowdens.txt"
  params:
    out_pre= "results/calc_RoH/{dataset}"
  shell:
    """
    Rscript scripts/find_lowDensityRegs.R {input} {params.out_pre} 1e6 0.05
    cat regions/exclude_regions_hg19.txt {output.ldr} > {output.reg}
    """

rule get_snp_depth:
  input:
    hom = "results/calc_RoH/{dataset}_RoH_segsJoined.hom",
    bim = "data/{dataset}.bim"
  output:
    "results/calc_RoH/{dataset}_RoHdepth"
  benchmark:
    "benchmarks/{dataset}/get_snp_depth.txt"
  shell:
    """
    Rscript scripts/get_seg_depthv2.R {input.hom} {input.bim} > {output}
    """

rule outliers:
  input:
    roh = "results/calc_RoH/{dataset}_RoHdepth",
    txt = "results/calc_RoH/{dataset}_temp_regions.txt"
  output:
    outliers = "results/calc_RoH/{dataset}_RoHdepth_outliers",
    regions = "results/calc_RoH/{dataset}_RoHdepth_outlierRegions",
    txt = "results/calc_RoH/{dataset}_exclude_regions_hg19.txt"
  benchmark:
    "benchmarks/{dataset}/outliers.txt"
  shell:
    """
    Rscript scripts/get_outlier_SNPs.R {input.roh} > {output.outliers}
    Rscript scripts/get_outlier_regions.R {output.outliers} > {output.regions}
    cat {input.txt} {output.regions} > {output.txt}
    """

rule qc_hom:
  input:
    hom = "results/calc_RoH/{dataset}_RoH_segsJoined.hom",
    txt = "results/calc_RoH/{dataset}_exclude_regions_hg19.txt"
  output:
    "results/calc_RoH/{dataset}_RoH_segsJoined_QCed.hom"
  benchmark:
    "benchmarks/{dataset}/qc_hom.txt"
  shell:
    """
    Rscript scripts/filter_RoH-IBD_segs.R {input.hom} {input.txt} 0.85 > {output}
    """

# -----------------------------------------
# # Calculate matching to RoH distribution
# -----------------------------------------

rule bits_matching:
  input:
    bits_matching_inputs,
    "results/calc_RoH/{dataset}_RoH_segsJoined_QCed.hom"
  output:
    "results/other/{dataset}_GL_IBD_3500kb.txt"
  shell:
    """
    bash scripts/bits_matching.sh {DATASET}
    """
# -----------------
# Kinship Analysis
# -----------------

rule kinship:
  input:
    bed = "data/{dataset}.bed"
  output:
    multiext("results/other/{dataset}_HWE-MAFfiltered", ".bed", ".bim", ".fam"),
    "results/other/{dataset}.genome"
  params:
    pre_one="data/{dataset}",
    pre_two = "results/other/{dataset}"
  benchmark:
    "benchmarks/{dataset}/kinship.txt"
  shell:
    """
    plink --bfile {params.pre_one} --maf {MAF} --make-bed --out {params.pre_two}_HWE-MAFfiltered
    plink --bfile {params.pre_two}_HWE-MAFfiltered --genome --out {params.pre_two}
    """

# ---------------------
# Calculate concordance
# ---------------------

### this rule might need an input function
rule calc_concordance:
  input:
    bim = "results/other/{dataset}_HWE-MAFfiltered.bim",
    ibd = concordance_input
  output:
    "results/other/{dataset}_NormRMSE.txt"
  benchmark:
    "benchmarks/{dataset}/calc_concordance.txt"
  shell:
    """
    bash scripts/calculate_concordance.sh {DATASET}
    """

# ----------
# RUN IBDne
# ----------

rule choose_params:
  input:
    ibd = "results/other/{dataset}_GL_IBD_3500kb.txt",
    norm = "results/other/{dataset}_NormRMSE.txt"
  output:
    "results/other/{dataset}_GL_parameters.txt",
    "results/other/{dataset}_only_QCed.ibd"
  benchmark:
    "benchmarks/{dataset}/choose_params.txt"
  shell:
    """
    bash scripts/choose_params.sh {input.ibd} {input.norm} {wildcards.dataset}
    """

rule ibdne:
  input:
    ibd = "results/other/{dataset}_only_QCed.ibd",
    map = "prepare/all_chroms/{dataset}.allchr.phased.cm.map"
  output:
    multiext("results/IBDne_output/{dataset}_GERMLINE2_IBDNe_minibd-{cM}", ".log", ".pair.excl", ".region.excl", ".ne", ".boot")
  benchmark:
    "benchmarks/{dataset}/{cM}_ibdne.txt"
  shell:
    """
    cat {input.ibd} | java -jar progs/ibdne.07May18.6a4.jar map={input.map} out=results/IBDne_output/{wildcards.dataset}_GERMLINE2_IBDNe_minibd-{cM} mincm={cM} nthreads={IBDne_THREADs} nboots=80 gmax={GMAX}
    """
