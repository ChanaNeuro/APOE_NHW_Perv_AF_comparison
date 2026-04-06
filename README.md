# APOE_NHW_Perv_AF_comparison

APOE_NHW_Peruvian_apoe3/4 vs4/4
#!/bin/bash
============================================================
============================================================
============================================================
#Task #Step 1. Data extraction
#• Extract a genomic region of ±1 MB around APOE. #• Create two separate datasets: (Peruvians and NHWs (Non-Hispanic Whites)
#Step 2. Allele frequency analysis #• Calculate allele frequencies for markers within this APOE region. #• Compare frequencies between the two datasets.
#Step 3. Stratified carrier analysis #• From the extracted region, identify ε3/ε4 (e34) and ε4/ε4 (e44) carriers. #• Re-run allele frequency analysis restricted to these carrier groups.
============================================================
#files used for the Cohorts: Peruvians & NHWs #/hihg/studies/AD/data/PeADI1_3/PeADI-wgs.chr19.snp.indel.recalibrated.vcf.gz #PeADI.WGS.Jan2025.txt #NHW #gcad.qc.compact_filtered.r4.wgs.36361.GATK.2023.06.06.biallelic.genotypes.chr19.ALL.vcf.bgz #vcf_samples_nhw.txt
============================================================
============================================================
APOE ±1Mb Region Analysis Pipeline
Cohorts: Peruvians & NHWs
Steps:
1. Extract ±1Mb APOE region for each cohort
2. Identify APOE-defining SNPs (rs429358, rs7412) genotypes
3. Classify carriers: e3/e4 and e4/e4
4. Subset VCFs by carrier group
5. Run enriched allele frequency analysis
6. Merge Peru vs NHW results for each carrier group
============================================================
You can check the available bcftools and load it
##frist run "module avail bcftools", then run teh second command and then the thirsd command module avail bcftools module load bcftools which bcftools
-----------------------------
INPUT FILES (must be bgzipped & indexed)
-----------------------------
# Peruvians
bcftools view
--regions chr19:43905781-45909394
PeADI-wgs.chr19.snp.indel.recalibrated.vcf.gz
-Oz -o apoe_pm1mb_peruvian.vcf.gz bcftools index apoe_pm1mb_peruvian.vcf.gz
NHWs
bcftools view
--regions chr19:43905781-45909394
gcad.qc.compact_filtered.r4.wgs.36361.GATK.2023.06.06.biallelic.genotypes.chr19.ALL.vcf.bgz
-Oz -o apoe_pm1mb_nhw.vcf.gz bcftools index apoe_pm1mb_nhw.vcf.gz
apoe_pm1mb_peruvian.vcf.gz -> Peruvians, ±1Mb around APOE
apoe_pm1mb_nhw.vcf.gz -> NHWs, ±1Mb around APOE
-----------------------------
SNP COORDINATES (GRCh38 build)
-----------------------------
rs429358: chr19:44908684 (REF=T, ALT=C)
rs7412: chr19:44908822 (REF=C, ALT=T)
============================================================
1. Extract genotypes for APOE SNPs in each cohort
============================================================
Peruvians
bcftools query -f '[%SAMPLE\t%GT\n]' -r chr19:44908684 apoe_pm1mb_peruvian.vcf.gz > rs429358_peru.txt bcftools query -f '[%SAMPLE\t%GT\n]' -r chr19:44908822 apoe_pm1mb_peruvian.vcf.gz > rs7412_peru.txt join -t $'\t' <(sort rs429358_peru.txt) <(sort rs7412_peru.txt) > apoe_genotypes_peru.txt
NHWs
bcftools query -f '[%SAMPLE\t%GT\n]' -r chr19:44908684 apoe_pm1mb_nhw.vcf.gz > rs429358_nhw.txt bcftools query -f '[%SAMPLE\t%GT\n]' -r chr19:44908822 apoe_pm1mb_nhw.vcf.gz > rs7412_nhw.txt join -t $'\t' <(sort rs429358_nhw.txt) <(sort rs7412_nhw.txt) > apoe_genotypes_nhw.txt
============================================================
2. Classify carriers (e3/e4 and e4/e4) for each cohort
============================================================
Function to classify carriers
classify_carriers () { local infile=$1 local e34_out=$2 local e44_out=$3
awk 'BEGIN{OFS="\t"} { sample=$1 g1=$2 # rs429358 genotype g2=$3 # rs7412 genotype
Plain Text
# Convert numeric GT to actual alleles
split(g1,a,"/")
split(g2,b,"/")
allele1_rs429358 = (a[1]==0 ? "T" : "C")
allele2_rs429358 = (a[2]==0 ? "T" : "C")
allele1_rs7412   = (b[1]==0 ? "C" : "T")
allele2_rs7412   = (b[2]==0 ? "C" : "T")
# Count ε4 alleles
e4_count = 0
if (allele1_rs429358=="C" && allele1_rs7412=="C") e4_count++
if (allele2_rs429358=="C" && allele2_rs7412=="C") e4_count++
if (e4_count==1) print sample > "'$e34_out'"
else if (e4_count==2) print sample > "'$e44_out'"
 
}' "$infile" }
Peruvians
classify_carriers apoe_genotypes_peru.txt peru_e34_ids.txt peru_e44_ids.txt
NHWs
classify_carriers apoe_genotypes_nhw.txt nhw_e34_ids.txt nhw_e44_ids.txt
============================================================
3. Subset VCFs to each carrier group
============================================================
subset_vcf () { local ids=$1 local invcf=$2 local outvcf=$3 bcftools view --samples-file "$ids" "$invcf" -Oz -o "$outvcf" bcftools index "$outvcf" }
Peruvians
subset_vcf peru_e34_ids.txt apoe_pm1mb_peruvian.vcf.gz apoe_pm1mb_peru_e34.vcf.gz subset_vcf peru_e44_ids.txt apoe_pm1mb_peruvian.vcf.gz apoe_pm1mb_peru_e44.vcf.gz
NHWs
subset_vcf nhw_e34_ids.txt apoe_pm1mb_nhw.vcf.gz apoe_pm1mb_nhw_e34.vcf.gz subset_vcf nhw_e44_ids.txt apoe_pm1mb_nhw.vcf.gz apoe_pm1mb_nhw_e44.vcf.gz
============================================================
4. Enriched allele frequency analysis for each subset
============================================================
freq_analysis () { local invcf=$1 local outfile=$2 bcftools +fill-tags "$invcf" -- -t AC,AN,AF,MAF,F_MISSING
| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/AC\t%INFO/AN\t%INFO/AF\t%INFO/MAF\t%INFO/F_MISSING\n'
| awk 'BEGIN{OFS="\t"; print "VAR_ID","RSID","AC","AN","AF","REF_AF","N_SAMPLES","MAF","F_MISSING"} { varid=$1":"$2":"$3":"$4 refaf=1-$8 nsamp=$7/2 print varid,$5,$6,$7,$8,refaf,nsamp,$9,$10 }' > "$outfile" }
Peruvians
freq_analysis apoe_pm1mb_peru_e34.vcf.gz freq_peru_e34_detailed.tsv freq_analysis apoe_pm1mb_peru_e44.vcf.gz freq_peru_e44_detailed.tsv
NHWs
freq_analysis apoe_pm1mb_nhw_e34.vcf.gz freq_nhw_e34_detailed.tsv freq_analysis apoe_pm1mb_nhw_e44.vcf.gz freq_nhw_e44_detailed.tsv
============================================================
5. Merge Peru vs NHW results for each carrier group
============================================================
merge_freqs () { local peru_file=$1 local nhw_file=$2 local outfile=$3
tail -n +2 "$peru_file" | sort -k1,1 > peru_sorted.tmp tail -n +2 "$nhw_file" | sort -k1,1 > nhw_sorted.tmp
join -t $'\t' peru_sorted.tmp nhw_sorted.tmp
| awk 'BEGIN{ OFS="\t"; print "VAR_ID","RSID_Peru","AC_Peru","AN_Peru","AF_Peru","REF_AF_Peru","N_SAMPLES_Peru","MAF_Peru","F_MISSING_Peru", "RSID_NHW","AC_NHW","AN_NHW","AF_NHW","REF_AF_NHW","N_SAMPLES_NHW","MAF_NHW","F_MISSING_NHW", "AF_Diff","REF_AF_Diff" } { af_diff=$5 - $13 refaf_diff=$6 - $14 print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,af_diff,refaf_diff }' > "$outfile"
rm peru_sorted.tmp nhw_sorted.tmp }
e3/e4 merge
merge_freqs freq_peru_e34_detailed.tsv freq_nhw_e34_detailed.tsv freq_e34_peru_nhw_merged.tsv
e4/e4 merge
merge_freqs freq_peru_e44_detailed.tsv freq_nhw_e44_detailed.tsv freq_e44_peru_nhw_merged.tsv
============================================================
OUTPUT FILES
============================================================
apoe_genotypes_peru.txt / apoe_genotypes_nhw.txt -> genotype table for rs429358 & rs7412
peru_e34_ids.txt / peru_e44_ids.txt -> sample IDs for each carrier group (Peru)
nhw_e34_ids.txt / nhw_e44_ids.txt -> sample IDs for each carrier group (NHW)
apoe_pm1mb_peru_e34.vcf.gz / apoe_pm1mb_peru_e44.vcf.gz -> subset VCFs (Peru)
apoe_pm1mb_nhw_e34.vcf.gz / apoe_pm1mb_nhw_e44.vcf.gz -> subset VCFs (NHW)
