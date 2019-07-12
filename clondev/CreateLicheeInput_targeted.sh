


source ReadConfig.sh $1

WORKDIR=${WORKDIR}/targeted/clondev/prepare_lichee_input/
SAMPLELIST=${ORIDIR}/${SAMPLELIST}

#ORIDIR=$1
#WORKDIR=$2
#PATIENT=$3


# Get the list of patients

PATIENTS=`cut -f2 ${SAMPLELIST} | sort | uniq | awk '{printf $0" "}'`
PATIENTS="AC1"
PATIENTS="AC1 AC13 AC14 AC30 AC35 AC6"
#PATIENTS="AC30"
PATIENTS="AC30 AC35 AC6"
# Run the rest of the script for each patient 

for PATIENT in $PATIENTS; do

# Combine all samples from the patient

## Get only healthy samples by removing those that are AC<patient_number>c or AC<patient_number>d

SAMPLES=`awk -v patient=$PATIENT '{if($2==patient){print $0}}' ${SAMPLELIST} | grep -v "[0-9]c" | grep -v "[0-9]d" | awk '{printf $1" "}'`

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


# Create bed with the non-diploid copy number regions


## Extract non-diploid genes for that patient

GENES=`grep ";$PATIENT;" ${ORIDIR}/cn_targeted/targeted_cn.csv | cut -d ";" -f1 | awk '{printf $1" "}'` 


## Create bed file from the genes

CNFILE=${WORKDIR}/${PATIENT}.cn.bed

if [ -f ${CNFILE} ]; then echo "Existing file removed"; rm $CNFILE; fi

awk 'BEGIN{print "chromosome\tstart\tend"}' > $CNFILE 

for GENE in $GENES; do
grep $GENE ${ORIDIR}/cn_targeted/design.bed | cut -f1,2,3  >> $CNFILE
done



# Remove variants in non-diploid regions

module load gcc/6.4.0 vcftools/0.1.15


vcftools \
    	--vcf ${WORKDIR}/${PATIENT}.combined.vcf \
	--exclude-bed ${WORKDIR}/${PATIENT}.cn.bed \
	--out ${WORKDIR}/${PATIENT}.diploid \
	--recode

mv ${WORKDIR}/${PATIENT}.diploid.recode.vcf ${WORKDIR}/${PATIENT}.diploid.vcf


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

## Locate the corresponding bam file

for SAMPLE in $SAMPLES; do

BAMFILE=`find $ORIDIR/bams_targeted -name "*bam" | grep -e "$SAMPLE[_|.]"`


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


# Convert to LICHeE input format

HEADER=`echo "$SAMPLES" | sed 's/ /\t/g'  | awk -v samples="$SAMPLES" '{print "#chr\tposition\tdescription\tNormal\t"samples}'`

awk -v mindepth=20 -v samples="$SAMPLES" -v header="${HEADER}" -v nsamples=${#SAMPLESARRAY[@]} -F " " '
BEGIN{
split(samples,samplearray," ")
printf "#chr\tposition\tdescription\tNormal"
for (sample in samplearray){
printf "\t"samplearray[sample]
}
print ""
}
{if(NR==1){next}



split("",noN)
split(samples,samplearray," ")

k=5

for (sample in samplearray){

if($k!="N"){noN[k]=$k}
k=k+4
}



firstnuc=noN[5]

for (altnuc in noN){


if(firstnuc!=noN[altnuc]){next}
}




i=2
for (sample in samplearray){
j=i+1
if(($i+$j) < 20){next}
i=i+4
}

split($1,a,":");
printf a[1]"\t"a[2]"\t"$4"/"firstnuc"\t0.0"



i=2 
for (sample in samplearray){
j=i+1
printf "\t"$j/($i+$j)
i=i+4
}

print ""

}' ${WORKDIR}/${PATIENT}.Counts > ${WORKDIR}/${PATIENT}.LicheeInput



done
