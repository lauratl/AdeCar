


source ReadConfig.sh $1

WORKDIR=${WORKDIR}/targeted/clondev/prepare_lichee_input/
SAMPLELIST=${ORIDIR}/${SAMPLELIST}
CONTROL=${ORIDIR}/${CONTROL}
LOG=${WORKDIR}/CreateLicheeInput_Targeted.log
TABLE=${WORKDIR}/CreateLicheeInput_Targeted.table



# Get the list of patients

PATIENTS=`cut -f2 ${SAMPLELIST} | sort | uniq | awk '{printf $0" "}'`



# Run the rest of the script for each patient 

echo "Starting Lichee input creation" > $LOG
echo "#patient;healthy;initial_variants;non-diploid-genes;diploidvariants;snvs;retrieved;enoughdepth;somatic;final" > $TABLE


for PATIENT in $PATIENTS; do

echo "Processing patient $PATIENT" >> $LOG


# Combine all samples from the patient


## Get also healthy sample if available (to remove the germline variants later)

SAMPLES=`awk -v patient=$PATIENT '{if($2==patient){printf $1" "}}' ${SAMPLELIST}`
HEALTHY=`awk -v patient=$PATIENT '{if($2==patient){print $1}}' ${CONTROL}`

echo "Samples are $SAMPLES" >> $LOG
echo "Healthy sample is $HEALTHY" >>$LOG

## Combine vcfs

module load gatk/3.7-0-gcfedb67 

SAMPLES_TO_COMBINE=`echo $SAMPLES | sed 's/ /\n/g' | awk -v workdir=$ORIDIR/variants_targeted/${PATIENT} '{split("abcdefghijklmnopqrstu",foosamplenames,"");printf "-V:"foosamplenames[NR]" "workdir"/"$1".vcf "}'`
ROD_PRIORITY_LIST=`echo $SAMPLES | sed 's/ /\n/g' | awk -v workdir=$ORIDIR/variants_targeted/${PATIENT} '{split("abcdefghijklmnopqrstu",foosamplenames,"");printf foosamplenames[NR]","}' | awk -F "" '{for(i=1;i<NF;i++){printf $i}}' `


java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CombineVariants \
	$SAMPLES_TO_COMBINE \
    	-o ${WORKDIR}/${PATIENT}.combined.vcf \
	-R ${RESDIR}/hg19.fasta \
	--genotypemergeoption PRIORITIZE \
	--rod_priority_list $ROD_PRIORITY_LIST


INITIAL=$(grep -v '^#' ${WORKDIR}/${PATIENT}.combined.vcf | wc -l)
echo "Initial variants: $INITIAL" >> $LOG


# Remove variants in non-diploid regions

## Extract non-diploid genes for that patient

GENES=`grep ";$PATIENT;" ${ORIDIR}/cn_targeted/targeted_cn.csv | cut -d ";" -f1 | awk '{printf $1" "}'` 

## Create bed file of the non-diploid regions from the genes

CNFILE=${WORKDIR}/${PATIENT}.cn.bed

if [ -f ${CNFILE} ]; then echo "Existing file removed"; rm $CNFILE; fi

awk 'BEGIN{print "chromosome\tstart\tend"}' > $CNFILE 

for GENE in $GENES; do
grep $GENE ${ORIDIR}/cn_targeted/design.bed | cut -f1,2,3  >> $CNFILE
done

echo "Genes with non-diploid copy number: $GENES" >>  $LOG


## Remove variants in non-diploid regions from the vcf

module load gcc/6.4.0 vcftools/0.1.15


vcftools \
    	--vcf ${WORKDIR}/${PATIENT}.combined.vcf \
	--exclude-bed ${WORKDIR}/${PATIENT}.cn.bed \
	--out ${WORKDIR}/${PATIENT}.diploid \
	--recode

mv ${WORKDIR}/${PATIENT}.diploid.recode.vcf ${WORKDIR}/${PATIENT}.diploid.vcf

DIPLOID=$(grep -v '^#' ${WORKDIR}/${PATIENT}.diploid.vcf | wc -l)
echo "Diploid variants: $DIPLOID" >> $LOG


# Remove indels

module load gatk/4.1.1.0

gatk SelectVariants \
        -V ${WORKDIR}/${PATIENT}.diploid.vcf \
        -R ${RESDIR}/hg19.fasta \
         --select-type-to-include SNP \
        -O ${WORKDIR}/${PATIENT}.diploid.snvs.vcf


SNVS=$(grep -v '^#' ${WORKDIR}/${PATIENT}.diploid.snvs.vcf | wc -l)
echo "Diploid SNVs: $SNVS" >> $LOG


# If no healthy sample, remove variants in dbSNP

if [ "$HEALTHY" == "" ]; then 


gatk SelectVariants \
	-V ${WORKDIR}/${PATIENT}.diploid.snvs.vcf \
	--discordance ${RESDIR}/dbsnp_138.hg19.vcf.gz \
	-R ${RESDIR}/hg19.fasta \
	-O ${WORKDIR}/${PATIENT}.diploid.snvs.nodbsnp.vcf

NEXTINPUT=${WORKDIR}/${PATIENT}.diploid.snvs.nodbsnp.vcf

else

NEXTINPUT=${WORKDIR}/${PATIENT}.diploid.snvs.vcf

fi


# Recover read counts

## Get list of positions to recover read counts

grep -v '^#' ${NEXTINPUT} | awk '{print $1"\t"$2-1"\t"$2}' > ${WORKDIR}/${PATIENT}.pos.bed


## Locate the corresponding bam file

for SAMPLE in $SAMPLES; do

BAMFILE=`find $ORIDIR/bams_targeted_new -name "*bam" | grep -e "$SAMPLE[_|.]"`


## Get read counts from the bam

gatk CollectAllelicCounts \
          -I ${BAMFILE} \
          -R ${RESDIR}/Hg19IonTorrentDefault.fa \
          -L ${WORKDIR}/${PATIENT}.pos.bed \
          -O ${WORKDIR}/${PATIENT}.${SAMPLE}.allelicCounts.tsv

done



# Merge read counts from the different samples

for SAMPLE in $SAMPLES; do

grep -v '^@' ${WORKDIR}/${PATIENT}.${SAMPLE}.allelicCounts.tsv | sed "s/COUNT/COUNT_$SAMPLE/g" | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$6 }' > ${WORKDIR}/${PATIENT}.${SAMPLE}.Counts.tsv

done

SAMPLESARRAY=($SAMPLES)



for SAMPLE in $SAMPLES; do

if [ $SAMPLE == ${SAMPLESARRAY[0]} ]; then
echo ""

elif [ $SAMPLE == ${SAMPLESARRAY[1]} ]; then
join ${WORKDIR}/${PATIENT}.${SAMPLESARRAY[0]}.Counts.tsv ${WORKDIR}/${PATIENT}.${SAMPLE}.Counts.tsv > ${WORKDIR}/${PATIENT}.tmp.${SAMPLE}.Counts 
PREVSAMPLE=${SAMPLE}

else
join ${WORKDIR}/${PATIENT}.tmp.${PREVSAMPLE}.Counts ${WORKDIR}/${PATIENT}.${SAMPLE}.Counts.tsv > ${WORKDIR}/${PATIENT}.tmp.${SAMPLE}.Counts
PREVSAMPLE=${SAMPLE}

fi
done

mv ${WORKDIR}/${PATIENT}.tmp.${PREVSAMPLE}.Counts ${WORKDIR}/${PATIENT}.Counts
rm ${WORKDIR}/${PATIENT}.tmp.*.Counts


RETRIEVED=$(grep -v '^CONTIG' ${WORKDIR}/${PATIENT}.Counts | wc -l)
echo "Retrieved read counts from  $RETRIEVED variants" >> $LOG


# Convert to LICHeE input format

module load gcccore/6.4.0 python/2.7.15

python CreateLicheeInput_targeted.py ${WORKDIR}/${PATIENT}.Counts ${WORKDIR}/${PATIENT}.LicheeInput $HEALTHY



FINAL=$(grep -v '^#' ${WORKDIR}/${PATIENT}.LicheeInput | wc -l)
echo "Finally kept $FINAL SNVs" >> $LOG

echo "$PATIENT;$HEALTHY;$INITIAL;$GENES;$DIPLOID;$SNVS;$RETRIEVED;$ENOUGH;$SOMATIC;$FINAL" >> $TABLE
done
