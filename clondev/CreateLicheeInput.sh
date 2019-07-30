


source ReadConfig.sh $1
PATIENT=$2
WORKDIR=${WORKDIR}/wes/clondev/prepare_clonefinder_input/from_bams/with_python

#ORIDIR=$1
#WORKDIR=$2
#PATIENT=$3


#echo "Usage: $0 <raw_dir> <working_dir> <patient>"



tumors="adenoma carcinoma"


# Filter PASS variants


orisamples="${PATIENT}a_WES_normalFitlered.vcf ${PATIENT}b_WES_normalFitlered.vcf"

for f in $orisamples; do
awk '
{split($0,linea,"");
if(linea[1]=="#"){
	print $0}
else{if($7=="PASS"){
	print $0}}}
' ${ORIDIR}/variants_wes/$f > ${WORKDIR}/${f/.vcf/.PASS.vcf}

done

# Rename samples in the VCF header

module load gcc/6.4.0 bcftools/1.9

bcftools reheader \
	--samples <(printf "%s\n%s\n" ${PATIENT}_healthy ${PATIENT}_adenoma) \
	--output ${WORKDIR}/${PATIENT}_adenoma.malformed.vcf \
	${WORKDIR}/${PATIENT}b_WES_normalFitlered.PASS.vcf

bcftools reheader \
	--samples <(printf "%s\n%s\n" ${PATIENT}_healthy ${PATIENT}_carcinoma) \
	--output ${WORKDIR}/${PATIENT}_carcinoma.malformed.vcf \
	${WORKDIR}/${PATIENT}a_WES_normalFitlered.PASS.vcf

# Edit malformed line of the header


for TUMOR in $tumors; do

sed 's/END<POLARITY,Number=0,Type=null,Description="null"/END,POLARITY"/g' ${WORKDIR}/${PATIENT}_${TUMOR}.malformed.vcf > ${WORKDIR}/${PATIENT}_${TUMOR}.vcf

done


# Combine adenoma and carcinoma


module load gatk/3.7-0-gcfedb67

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CombineVariants \
    	-V:adenoma ${WORKDIR}/${PATIENT}_adenoma.vcf \
    	-V:carcinoma ${WORKDIR}/${PATIENT}_carcinoma.vcf \
    	-o ${WORKDIR}/${PATIENT}.missingreadcounts.vcf \
	-R ${RESDIR}/hs37d5.fa \
	--genotypemergeoption PRIORITIZE \
	--rod_priority_list adenoma,carcinoma


# Create bed with the non-diploid copy number regions

awk '
{if($6=="cn" || $6!=2){
	print $1"\t"$2"\t"$3}}
' ${ORIDIR}/cn_wes/${PATIENT}a_call.cns > ${WORKDIR}/${PATIENT}_carcinoma.cn.bed

awk '
{if($6=="cn" || $6!=2){
	print $1"\t"$2"\t"$3}}
' ${ORIDIR}/cn_wes/${PATIENT}b_call.cns > ${WORKDIR}/${PATIENT}_adenoma.cn.bed



 # Remove variants in non-diploid regions


module load gcc/6.4.0 vcftools/0.1.15


vcftools \
  --vcf ${WORKDIR}/${PATIENT}.missingreadcounts.vcf \
	--exclude-bed ${WORKDIR}/${PATIENT}_adenoma.cn.bed \
	--out ${WORKDIR}/${PATIENT}.tmp \
	--recode

vcftools \
  --vcf ${WORKDIR}/${PATIENT}.tmp.recode.vcf \
	--exclude-bed ${WORKDIR}/${PATIENT}_carcinoma.cn.bed \
	--out ${WORKDIR}/${PATIENT}.diploid \
	--recode

rm ${WORKDIR}/${PATIENT}.tmp.recode.vcf
mv ${WORKDIR}/${PATIENT}.diploid.recode.vcf ${WORKDIR}/${PATIENT}.diploid.vcf


# Remove indels

module load gatk/4.1.1.0

gatk SelectVariants \
        -V ${WORKDIR}/${PATIENT}.diploid.vcf \
        -R ${RESDIR}/hs37d5.fa \
         --select-type-to-include SNP \
        -O ${WORKDIR}/${PATIENT}.diploid.snvs.vcf



# Get list of positions to recover read counts

grep -v '^#' ${WORKDIR}/${PATIENT}.diploid.snvs.vcf | awk '{print $1"\t"$2-1"\t"$2}' > ${WORKDIR}/${PATIENT}.pos.bed


# Recover read counts



gatk CollectAllelicCounts \
          -I ${ORIDIR}/bams_wes/${PATIENT}Adenoma_recalibrated.bam \
          -R ${RESDIR}/hs37d5.fa \
          -L ${WORKDIR}/${PATIENT}.pos.bed \
          -O ${WORKDIR}/${PATIENT}_adenoma.allelicCounts.tsv

gatk CollectAllelicCounts \
          -I ${ORIDIR}/bams_wes/${PATIENT}Carcinoma_recalibrated.bam \
          -R ${RESDIR}/hs37d5.fa \
          -L ${WORKDIR}/${PATIENT}.pos.bed \
          -O ${WORKDIR}/${PATIENT}_carcinoma.allelicCounts.tsv


# Merge adenoma and carcinoma read counts


for TUMOR in $tumors; do

grep -v '^@' ${WORKDIR}/${PATIENT}_${TUMOR}.allelicCounts.tsv | sed "s/COUNT/COUNT_$TUMOR/g" | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$6 }' > ${WORKDIR}/${PATIENT}_${TUMOR}.Counts.tsv

done

join ${WORKDIR}/${PATIENT}_adenoma.Counts.tsv ${WORKDIR}/${PATIENT}_carcinoma.Counts.tsv > ${WORKDIR}/${PATIENT}.Counts


# Convert to LICHeE and CloneFinder input format

module load gcccore/6.4.0 python/2.7.15

python CreateCDInput_targeted.py \
	--input ${WORKDIR}/${PATIENT}.Counts \
	--lichee ${WORKDIR}/${PATIENT}.LicheeInput \
	--cloneFinder ${WORKDIR}/${PATIENT}.CloneFinderInput \
	--healthy "$HEALTHY" \
	--maxVafIfNotHealthy 0.9 \
	--minDepth 20 \
	--minVaf 0.05 \
	--germlineVaf 0.1
