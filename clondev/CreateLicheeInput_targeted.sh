


source ReadConfig.sh $1

WORKDIR=${WORKDIR}/targeted/clondev/prepare_lichee_input/
SAMPLELIST=${ORIDIR}/${SAMPLELIST}
CONTROL=${ORIDIR}/${CONTROL}
LOG=${WORKDIR}/CreateLicheeInput_Targeted.log
TABLE=${WORKDIR}/CreateLicheeInput_Targeted.table

#ORIDIR=$1
#WORKDIR=$2
#PATIENT=$3


# Get the list of patients

PATIENTS=`cut -f2 ${SAMPLELIST} | sort | uniq | awk '{printf $0" "}'`
#PATIENTS="AC1"
#PATIENTS="AC1 AC13 AC14 AC30 AC35 AC6"
PATIENTS="AC30"
#PATIENTS="AC30 AC35 AC6"
# Run the rest of the script for each patient 

echo "Starting Lichee input creation" > $LOG
echo "#patient;healthy;initial_variants;non-diploid-genes;diploidvariants;snvs;retrieved;enoughdepth;somatic;final" > $TABLE


for PATIENT in $PATIENTS; do

echo "Processing patient $PATIENT" >> $LOG

# Combine all samples from the patient

## Get only healthy samples by removing those that are AC<patient_number>c or AC<patient_number>d

#SAMPLES=`awk -v patient=$PATIENT '{if($2==patient){print $0}}' ${SAMPLELIST} | grep -v "[0-9]c" | grep -v "[0-9]d" | awk '{printf $1" "}'`

## Get also healthy sample if available (to remove the germline variants)

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

echo "Genes with non-diploid copy number: $GENES" >>  $LOG


# Remove variants in non-diploid regions

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


awk '
{split($0,a,"");
if(a[1]=="#"){
	print $0}
else{
	if(length($4)==1 && length($5)==1){
		print $0}
	else{
		split($5,b,",")
		for(alt in b){if(length(b[alt])!=1){next}}
		print $0
}}}
' ${WORKDIR}/${PATIENT}.diploid.vcf > ${WORKDIR}/${PATIENT}.diploid.snvs.vcf

SNVS=$(grep -v '^#' ${WORKDIR}/${PATIENT}.diploid.snvs.vcf | wc -l)
echo "Diploid SNVs: $SNVS" >> $LOG

# Get list of positions to recover read counts

grep -v '^#' ${WORKDIR}/${PATIENT}.diploid.snvs.vcf | awk '{print $1"\t"$2-1"\t"$2}' > ${WORKDIR}/${PATIENT}.pos.bed


# Recover read counts

module load gatk/4.1.1.0

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

HEADER=`echo "$SAMPLES" | sed 's/ /\t/g'  | awk -v samples="$SAMPLES" '{print "#chr\tposition\tdescription\tNormal\t"samples}'`

awk -v mindepth=20 -v minvaf=0.05 -v samples="$SAMPLES" -v header="${HEADER}" -v nsamples=${#SAMPLESARRAY[@]} -F " " '
{if(NR==1){
samplestring=""
for(i=2;i<=NF;i=i+4){
split($i,rc,"_")
samplestring=samplestring+"."+rc[3]"_"rc[4]"_"rc[5]
}
print samplestring"***"
split(samplestring,samplearray,".")

printf "#chr\tposition\tdescription\tNormal"
for (sample in samplearray){
printf "\t"samplearray[sample]
}
print ""
print length(samplearray)
}

split("",noN)

k=5

for (sample in samplearray){
#print sample
if($k!="N"){noN[k]=$k}
k=k+4
}




for (altnuc in noN){
firstnuc=noN[altnuc]
break
}

# Decide wether there are more than one alternative allele

flag=0
for (altnuc in noN){
if(firstnuc!=noN[altnuc]){
flag=1; uno =firstnuc; dos =noN[altnuc]
}
}

# If more than one alternative, check if one is 5-fold the other

if(flag==1){
unocount=0
doscount=0

for(j=5;j<=NF;j=j+4){
if($j==uno){unocount+=$(j-2)}
else{doscount+=$(j-2)}
}
if(unocount>doscount){major=uno;minor=dos;if((doscount/unocount)>0.2){next}}
else{major=dos;minor=uno;if((unocount/doscount)>0.2){next}
}

# If keep the major, set the minor read counts to 0


for(j=5;j<=NF;j=j+4){
if($j==minor){$(j-2)=0; $j="N"}
}

print "*********"$0
}


i=2
for (sample in samplearray){
j=i+1
if(($i+$j) < 20){next}
i=i+4
}

flag=0
i=2
for (sample in samplearray){
j=i+1

if($j/($i+$j) >=minvaf){flag=1}
i=i+4
}
if(flag==0){next}

split($1,a,":");
printf a[1]"\t"a[2]"\t"$4"/"firstnuc"\t0.0"



i=2 
for (sample in samplearray){
j=i+1
printf "\t"$j/($i+$j)
i=i+4
}

print ""

}' ${WORKDIR}/${PATIENT}.Counts > ${WORKDIR}/${PATIENT}.All

ENOUGH=$(grep -v '^#' ${WORKDIR}/${PATIENT}.All | wc -l)
echo "Kept $ENOUGH variants with enough depth" >> $LOG

# Filter out germline variants

awk -v healthy=$HEALTHY 'BEGIN{healthy_idx=0}{
if(NR==1){
for(i=5;i<=NF;i++){

if($i==healthy){healthy_idx=i}
}
if(healthy_idx==0){print $0}
else{printf $1"\t"$2"\t"$3"\t"$healthy_idx;
for(i=5;i<=NF;i++){
if(i!=healthy_idx){printf "\t"$i}
}
print healthy_idx "is healthy_idx"
print healthy " is healthy"
print ""
}
}
else{
if(healthy_idx==0){print $0
}else{
if($healthy_idx<0.2){
printf $1"\t"$2"\t"$3"\t"$healthy_idx;
for(i=5;i<=NF;i++){
if(i!=healthy_idx){printf "\t"$i}
}
print ""
}


}
}}' ${WORKDIR}/${PATIENT}.All > ${WORKDIR}/${PATIENT}.Somatic

SOMATIC=$(grep -v '^#' ${WORKDIR}/${PATIENT}.Somatic | wc -l)
echo "Kept $SOMATIC somatic SNVs" >> $LOG




awk '{
if($1=="chr19" && $2=="15285135"){next}
if($1=="chr21" && $2=="326385"){next}
if($1=="chr8" && $2=="13356818"){next}
print $0
}
' ${WORKDIR}/${PATIENT}.Somatic > ${WORKDIR}/${PATIENT}.LicheeInput

FINAL=$(grep -v '^#' ${WORKDIR}/${PATIENT}.LicheeInput | wc -l)
echo "Finally kept $FINAL SNVs" >> $LOG

echo "$PATIENT;$HEALTHY;$INITIAL;$GENES;$DIPLOID;$SNVS;$RETRIEVED;$ENOUGH;$SOMATIC;$FINAL" >> $TABLE
done
