


source ReadConfig.sh $1
PATIENT=$2
WORKDIR=${WORKDIR}/prepare_clonefinder_input/from_bams

#ORIDIR=$1
#WORKDIR=$2
#PATIENT=$3


#echo "Usage: $0 <raw_dir> <working_dir> <patient>"



tumors="adenoma carcinoma"


# Filter PASS variants

 
f=AC34a_WES_normalFitlered.vcf

awk '
{split($0,linea,"");
if(linea[1]=="#"){
	print $0}
else{if($7=="PASS"){
	print $0}}}
' ${ORIDIR}/variants_wes/$f > ${WORKDIR}/${f/.vcf/.PASS.vcf}


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
rm ${WORKDIR}/${PATIENT}.tmp
mv ${PATIENT}.diploid.recode.vcf ${PATIENT}.diploid.vcf


# Remove indels


awk '
{split($0,a,"");
if(a[1]=="#"){
	print $0}
else{
	if(length($4)==1 && length($5)==1){
		print $0}}}
' ${WORKDIR}/${PATIENT}.diploid.vcf > ${WORKDIR}/${PATIENT}.diploid.snvs.vcf


# Get list of positions to recover read counts

grep -v '^#' ${WORKDIR}/${PATIENT}.diploid.snvs.vcf | awk '{print $1"\t"$2-1"\t"$2}' > ${WORKDIR}/${PATIENT}.pos.bed


# Recover read counts

module load gatk/4.1.1.0

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

# Convert to LICHeE input format

awk -F " " '
BEGIN{print "#chr\tposition\tdescription\tNormal\tAdenoma\tCarcinoma"}
{if(NR==1){next}
if($5!=$9){
	if($5!="N" && $9!="N"){
		next}
	else{
		if($5!="N"){
			alt=$5}
		else{
			alt=$9}}}
else{alt=$5}
split($1,a,":");
print a[1]"\t"a[2]"\t"$4"/"$5"\t0.0\t"$3/($2+$3)"\t"$7/($6+$7)

}' ${WORKDIR}/${PATIENT}.counts > ${WORKDIR}/${PATIENT}.LicheeInput
