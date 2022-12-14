#!/bin/bash
source /etc/environment
source /etc/profile

PIPEVERSION="3.1.0 vispa2"
STARTTIME=`date +'%y-%m-%d %H:%M:%S'`
DBHOSTIP="localhost"
STATSDBNAME='stats_summary_aav'

echo "
  +--------------------------------------------------------+
  |                                                        |
  |             Illumina Pipeline: VISPA 2 PE              |
  |                                                        |
  +--------------------------------------------------------+
  |  Author:   Andrea Calabria                             |
  |  Date:     April 2019                                  |
  |  Contact:  calabria.andrea@hsr.it                      |
  |  Version:  1.0 - based on VISPA2 3 - PE                |
  +--------------------------------------------------------+

  REQUIRED VARS and relative ORDER POSITION -> REMEMBER NO SPACES!!!!!!!!!
	1. disease id (MLD)
	2. patient id (MLD02)
	3. local server working path [/opt/NGS/results] -> it will be added DISEASE and PATIENT ID [/opt/NGS/results/MLD/MLD02/]
	4. R1
	5. R2
	6. pool id
	7. barcode list LTR
	8. barcode list LC
	9. reference genome (NB: indexed by BWA, hg19:: /opt/NGS/index/human/hg19/bwa/hg19.nochr0)
	10. tmp dir
	11. association file
	12. db host id
	13. db target schema
	14. db target table
	15. PhiX genome (bwa index, i.e.: PHIX=/opt/NGS/index/phix/phix174/bwa/phiX174.fasta)
	16. Vector genome (i.e.: LV 751: LVGENOME=/opt/NGS/index/virus/lv751wg/bwa/lv.fa)
	17. Contaminant FastA DB (i.e: /opt/NGS/vectors/UniVec_Tiget.fa)
	18. Number of MAXIMUM processes to apply to each step (max=number of CPU)
	19. GATK reference genome (i.e. /opt/NGS/index/human/hg19/bwa/hg19.nochr0.fa)
	20. Genome ID for CIGAR-MD program. NB: this option collects only whole genome so far, not single regions (even the program allows it!). Example: hg19, mm9, mfa5, ce, lv, lvarsa, lvwas, lvkana, lvamp, transposon, giada. [CIGARGENOMEID]
	21. Genome ID for CIGAR-MD program for VECTOR genome. Example: lv, lvarsa, lvwas, lvkana, lvamp, transposon, giada [VECTORCIGARGENOMEID]
	22. SUBOPTIMALTHRESHOLD to compare first vs second hit (AS vs XS tags), this is a percentage value (int).
	23. Remove TMP dir? values: {remove_tmp_yes, remove_tmp_no}
	24. LTR in forward
	25. LTR in reverse complement
	26. Linker Cassette in forward
	27. Linker Cassette in reverse
	28. FASTQ Quality Filter Method: {lam, slim, no_quality_filter}
	29. Repeats? values: {repeats_yes, repeats_no}
	30. Quality Filter Window (Quality in Phred Scale)

  NOTES:
  	- Input are now downloaded and processed locally as tmp files, then you can choose to delete them tmp files.
  	- LV is now removed from input data at before aligning to reference genome (e.g. human or mouse).
"
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Preprocessing input variables (delimiters:<>)"
## print input variables (check for log utils)
INPUTVARNUM=0
for INPUTVAR in "$@"; do
	let INPUTVARNUM++;
	printf -v INPUTNUM '%02d' $INPUTVARNUM;
    echo "  => Input Variable: Order number = <${INPUTNUM}> ; Var Content = <${INPUTVAR}>";
done

##### ================ start SETTINGS ====================== #####
GENOME="${9}";
SPECIE="`basename ${GENOME} | cut -f1 -d.`";
if [ ${SPECIE:0:2} = "hg" ]; then
	SPECIE="human";
else
	if [ ${SPECIE:0:2} = "mm" ]; then
	SPECIE="mouse";
	fi
fi

# LTR 32 bases
LTR="${24}";
LTR_rc="${25}";
LC_fwd="${26}";
LC_rev="${27}";
SEQCONTAM="/opt/applications/scripts/isatk/vector/known.seqs.contaminants.tsv"
# suboptimal threshold set up
SUBOPTIMALTHRESHOLD="${22}"

ASSOCIATIONFILE="${11}";
DBHOSTID="${12}";
DBSCHEMA="${13}"
DBTABLE="${14}"

DISEASE="${1}";
PATIENT="${2}";
POOL="${6}";

BARCODE_LTR="${7}";
BARCODE_LC="${8}";

R1_FASTQ="${4}";
R2_FASTQ="${5}";

NGSWORKINGPATH="${3}"; # WITHOUT LAST SLASH -> if not present, I will create it a new one with this name using mkdir
#OUTDIR_POOLS="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/pools";
OUTDIR_MERGE_BAM="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/bam";
OUTDIR_POOL_BAM="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/bam/"${POOL};
#OUTDIR_MERGE_SAM="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/sam";
#OUTDIR_POOL_SAM="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/sam/"${POOL};
OUTDIR_MERGE_BED="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/bed";
OUTDIR_POOL_BED="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/bed/"${POOL};
OUTDIR_MERGE_QUAL="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/quality";
OUTDIR_POOL_QUAL="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/quality/"${POOL};
OUTDIR_MERGE_QUANTIF="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/quantification";
OUTDIR_POOL_QUANTIF="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/quantification/"${POOL};
OUTDIR_MERGE_ISS="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/iss";
OUTDIR_POOL_ISS="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/iss/"${POOL};
OUTDIR_POOL_ISS_REPEATS="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/iss/"${POOL}/repeats;
OUTDIR_MERGE_BCMUXALL="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/bcmuxall";
OUTDIR_POOL_BCMUXALL="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/bcmuxall/"${POOL};
#OUTDIR_PHIX="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/phix"
#OUTDIR_VECTOR="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/vector"
#OUTDIR_VECTOR_POOL="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/vector/"${POOL};
TMPDIR="${10}";

FASTQ_QF="${28}";
REPEATS="${29}";
Q_FILTER="${30}";
FUSION_PRIMERS="${31}";


PHIXGENOME="${15}" ;
VECTORGENOME="${16}" ;
CONTAMINANTDB="${17}" ; # for example /opt/NGS/vectors/UniVec_Tiget.fa
GATKREFGENOME="${19}" ;
CIGARGENOMEID="${20}" ;
VECTORCIGARGENOMEID="${21}" ;
REMOVE_TMP_DIR="${23}";

MAXTHREADS="${18}" ;

#ANNOTATIONGTF=/home/carlo/workspace/bed/AAV-CBA-V5-mMecp2-iso1-UTR3.bed
BWA_MIN_ALN_LEN=30
SUBOPTIMALTHRESHOLD_rp='40';
minSM_fromN_rp=180
minSM_rp=30
##### ================ end SETTINGS ======================== #####

mkdir ${NGSWORKINGPATH}
mkdir ${NGSWORKINGPATH}/${DISEASE}
mkdir ${NGSWORKINGPATH}/${DISEASE}/${PATIENT}

#mkdir ${OUTDIR_POOLS};
mkdir ${OUTDIR_MERGE_BAM};
#mkdir ${OUTDIR_MERGE_SAM};
mkdir ${OUTDIR_MERGE_BED};
mkdir ${OUTDIR_MERGE_QUAL};
mkdir ${OUTDIR_POOL_QUAL};
#mkdir ${OUTDIR_PHIX};
#mkdir ${OUTDIR_VECTOR};
#mkdir ${OUTDIR_VECTOR_POOL};
mkdir ${OUTDIR_MERGE_QUANTIF};
mkdir ${OUTDIR_POOL_QUANTIF};
mkdir ${OUTDIR_MERGE_ISS};
mkdir ${OUTDIR_POOL_ISS};
mkdir ${OUTDIR_POOL_ISS_REPEATS};
mkdir ${OUTDIR_MERGE_BCMUXALL};
mkdir ${OUTDIR_POOL_BCMUXALL};

mkdir ${OUTDIR_POOL_BAM};
#mkdir ${OUTDIR_POOL_SAM};
mkdir ${OUTDIR_POOL_BED};
# TMPDIR="/storage/d2/ngs/pipelinetmpdir/breastcancer_sa46" ;
mkdir ${TMPDIR} ;

mkdir ${TMPDIR}/bcmuxall/
mkdir ${TMPDIR}/bed/
mkdir ${TMPDIR}/bam/
mkdir ${TMPDIR}/sam/
mkdir ${TMPDIR}/pools/
mkdir ${TMPDIR}/iss/
#mkdir ${OUTDIR_POOL_QUAL}/FragmentSize/;

OUTDIR_VECTOR=${TMPDIR}/pools
OUTDIR_VECTOR_POOL=${TMPDIR}/pools


echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Create stats_summary table into DB (if not exists)"
### create DB table (if not exists) of summary. it will update stats for all tags. and of run

mysql -u andrea --password=andrea -e "create database if not exists ${DBSCHEMA} ; " ;

mysql -u andrea --password=andrea -e "
	CREATE TABLE IF NOT EXISTS ${STATSDBNAME} (
		RUN_ID varchar(255) default null,
		RUN_NAME varchar(1000) default null,
		DISEASE varchar(255) default null,
		PATIENT varchar(255) default null,
		POOL varchar(255) default null,
		TAG varchar(255) default null,
		LTR_ID varchar(255) default null,
		LC_ID varchar(255) default null,
		RAW_READS int(30) DEFAULT NULL,
		QUALITY_PASSED int(30) DEFAULT NULL,
		PHIX_MAPPING int(30) DEFAULT NULL,
		PLASMID_MAPPED_BYPOOL int(30) DEFAULT NULL,
		RAW_NO_PLASMID int(30) DEFAULT NULL,
		BARCODE_MUX int(30) DEFAULT NULL,
		LTR_IDENTIFIED int(30) DEFAULT NULL,
		TRIMMING_LTRR1 int(30) DEFAULT NULL,
		TRIMMING_LTRR1R2 int(30) DEFAULT NULL,
		TRIMMING_LTRR1R2_LCR1 int(30) DEFAULT NULL,
		TRIMMING_FINAL_RESCUED int(30) DEFAULT NULL,
		TRIMMING_FINAL_LTRLC int(30) DEFAULT NULL,
		LV_MAPPED int(30) DEFAULT NULL,
		BWA_INPUT int(30) DEFAULT NULL,
		BWA_MAPPED int(30) DEFAULT NULL,
		BWA_MAPPED_PP int(30) DEFAULT NULL,
		BWA_MAPPED_ST int(30) DEFAULT NULL,
		BWA_MAPPED_OVERALL int(30) DEFAULT NULL,
		BWA_ALIGNED_R1 int(30) DEFAULT NULL,
		RECALIB_MAPPED int(30) DEFAULT NULL,
		RECALIB_MAPPED_PP int(30) DEFAULT NULL,
		RECALIB_MAPPED_ST int(30) DEFAULT NULL,
		RECALIB_MAPPED_OVERALL int(30) DEFAULT NULL,
		RECALIB_ALIGNED_R1 int(30) DEFAULT NULL,
		RECALIB_SOFCLIPPED_READS int(30) DEFAULT NULL,
		FILTER_MATE_TO_REMOVE int(30) DEFAULT NULL,
		FILTER_CIGARMD_TO_REMOVE int(30) DEFAULT NULL,
		FILTER_JOINT_MC_TO_REMOVE int(30) DEFAULT NULL,
		FILTER_JOINT_MC_PP int(30) DEFAULT NULL,
		FILTER_JOINT_MC_ST int(30) DEFAULT NULL,
		FILTER_JOINT_MC_OVERALL int(30) DEFAULT NULL,
		FILTER_JOINT_ALIGNED_R1 int(30) DEFAULT NULL,
		FILTER_ALMQUAL_PP int(30) DEFAULT NULL,
		FILTER_ALMQUAL_ST int(30) DEFAULT NULL,
		FILTER_ALMQUAL_OVERALL int(30) DEFAULT NULL,
		FILTER_ALMQUAL_ALIGNED_R1 int(30) DEFAULT NULL,
		ISS_FINAL int(30) DEFAULT NULL,
		ISS_MAPPED int(30) DEFAULT NULL,
		ISS_MAPPED_PP int(30) DEFAULT NULL,
		ISS_MAPPED_ST int(30) DEFAULT NULL,
		ISS_MAPPED_OVERALL int(30) DEFAULT NULL,
		ISS_ALIGNED_R1 int(30) DEFAULT NULL
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;
	" ${DBSCHEMA} ;

# checking vars
RUN_ID=`date +"%Y%m%d%H%M%S"`
RUN_NAME="${DISEASE}|${PATIENT}|${POOL}"

BCLTRf=`cut "${BARCODE_LTR}" -f1`;
BCLTR=(${BCLTRf}) ;
BCLCf=`cut "${BARCODE_LC}" -f1`;
BCLC=(${BCLCf}) ;

ASSOBCLIST=(`cut -f1 ${ASSOCIATIONFILE} | sort | uniq `); ## barcode list (as first colun of the association fTIGETile!) user defined!

##### ================ COPY DATA INTO TMP DIR ================== #####
cp ${R1_FASTQ} ${R2_FASTQ} ${TMPDIR}/
# identify base names
R1_NAME="`basename ${R1_FASTQ}`";
R2_NAME="`basename ${R2_FASTQ}`";
# now reassign raw input vars
R1_FASTQ=${TMPDIR}/${R1_NAME};
R2_FASTQ=${TMPDIR}/${R2_NAME};
RAW_READS=$((`zcat ${R1_FASTQ} | wc -l `/4));

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] New TMP files check"
ls -lh ${R1_FASTQ}
ls -lh ${R2_FASTQ}

##### ================ RAW DATA QUALITY ======================== #####
# 0. quality check ...! NOT WORKING FINE IIF FASTQ.GZ COMES FROM ILLUMINA MACHINE/SW...!!! YOU MUST CONVERT THEM ALL: gunzip first and then gzip them all again (their compression is not standard)
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Quality check -> FastQC reports and FastX stats"
fastqc --nogroup -o ${OUTDIR_POOL_QUAL} --contaminants ${SEQCONTAM} -t ${MAXTHREADS} -f fastq ${R1_FASTQ} ${R2_FASTQ}

# FASTA of RANDOM TAGs
if ! [ ${FASTQ_QF} = "no_quality_filter" ]; then
	if [ -z "$Q_FILTER" ]; then
		## remove low quality R1-R2 reads with fastq_qf program
		fastq_qf -a ${R1_FASTQ} -b ${R2_FASTQ} -o ${TMPDIR} -t ${MAXTHREADS} -m ${FASTQ_QF}
		R1_FASTQ="${TMPDIR}/QF.${R1_NAME}";
		R2_FASTQ="${TMPDIR}/QF.${R2_NAME}";
		QUALITY_PASSED=$((`zcat ${R1_FASTQ} | wc -l `/4));
		rm ${TMPDIR}/${R1_NAME} ${TMPDIR}/${R2_NAME};
	else
		## remove low quality R1-R2 reads with fastq_qf program
		fastq_qf -a ${R1_FASTQ} -b ${R2_FASTQ} -o ${TMPDIR} -t ${MAXTHREADS} -m ${FASTQ_QF} -f ${Q_FILTER}
		R1_FASTQ="${TMPDIR}/QF.${R1_NAME}";
		R2_FASTQ="${TMPDIR}/QF.${R2_NAME}";
		QUALITY_PASSED=$((`zcat ${R1_FASTQ} | wc -l `/4));
		rm ${TMPDIR}/${R1_NAME} ${TMPDIR}/${R2_NAME};
	fi
fi

# ##### ================ PHIX quantification START ======================== #####
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] PHIX Alignment to reference genome each single pair"

### old version (with sam output)
# bwa mem -k 16 -r 1 -M -T 15 -R "@RG\tID:Test\tSM:PHIX\tCN:TIGET" -t ${MAXTHREADS} ${PHIXGENOME} <(zcat ${R1_FASTQ} ) <(zcat ${R2_FASTQ} )> ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.phix.PE.sam
# samtools view -F 2308 -q 25 -f 35 -uS ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.phix.PE.sam > ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.phix.PE.bam;
# rm ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.phix.PE.sam

### NEW version (without sam output)
bwa mem -k 16 -r 1 -M -v 1 -T 15 -t ${MAXTHREADS} ${PHIXGENOME} <(zcat ${R1_FASTQ} ) <(zcat ${R2_FASTQ} ) | samtools view -F 2308 -q 25 -f 35 -uS - > ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.phix.PE.bam;

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filtering raw data"
samtools view ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.phix.PE.bam | cut -f 1 > ${TMPDIR}/PHIX.header.list;
sort --parallel=5 ${TMPDIR}/PHIX.header.list > ${TMPDIR}/PHIX.header.sorted.list;
rm ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.phix.PE.bam ${TMPDIR}/PHIX.header.list;
## rename files
BNAME_R1=`basename ${R1_FASTQ} | sed 's/.gz//g' | cut -d'.' -f1`;
BNAME_R2=`basename ${R2_FASTQ} | sed 's/.gz//g' | cut -d'.' -f1`;
zcat ${R1_FASTQ} | fqreverseextract_pureheader ${TMPDIR}/PHIX.header.sorted.list | pigz -f -c > ${TMPDIR}/${BNAME_R1}_R1_nophix.fastq.gz &
zcat ${R2_FASTQ} | fqreverseextract_pureheader ${TMPDIR}/PHIX.header.sorted.list | pigz -f -c > ${TMPDIR}/${BNAME_R2}_R2_nophix.fastq.gz &
wait
# checking VARS
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Checking Variables ";
PHIX_MAPPING=`wc -l ${TMPDIR}/PHIX.header.sorted.list | cut -d' ' -f1 ` ;
##PHIX_MAPPING="99999999" ;
## give new gzipped files without PHIX, change variables
rm ${R1_FASTQ} ${R2_FASTQ} ${TMPDIR}/PHIX.header.sorted.list
R1_FASTQ="${TMPDIR}/${BNAME_R1}_R1_nophix.fastq.gz" # BNAME_R1=`basename ${R1_FASTQ} | sed 's/.gz//g' | cut -d'.' -f1`;
R2_FASTQ="${TMPDIR}/${BNAME_R2}_R2_nophix.fastq.gz"
# ##### ================ PHIX quantification END  ======================== #####


##### ================ PLASMIDS quantification POOL based - start ======================== #####
BNAME_R1=`basename ${R1_FASTQ} | sed 's/.gz//g' | cut -d'.' -f1`;
BNAME_R2=`basename ${R2_FASTQ} | sed 's/.gz//g' | cut -d'.' -f1`;
for PLASMID in kana amp; do
	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Alignment to reference ${PLASMID} genome"
	## --- paired-end reads ---
	# low seed value due to small genome

	### old version (with sam output)
	# bwa mem -k 14 -r 1 -T 15 -c 1 -R "@RG\tID:plasmid\tSM:${PLASMID}\tCN:TIGET" -t ${MAXTHREADS} /opt/genome/vector/plasmid/plasmids.${PLASMID}.fa <(zcat ${R1_FASTQ} ) <(zcat ${R2_FASTQ} ) > ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${PLASMID}.sam
	# quality filter of seed-1 value
	# samtools view -F 2308 -q 20 -uS ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${PLASMID}.sam | samtools sort - ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${PLASMID}.sorted
	# rm ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${PLASMID}.sam

	### NEW version (without sam output)
	bwa mem -k 14 -r 1 -v 1 -T 15 -c 1 -t ${MAXTHREADS} /opt/genome/vector/plasmid/plasmids.${PLASMID}.fa <(zcat ${R1_FASTQ} ) <(zcat ${R2_FASTQ} ) | samtools view -F 2308 -q 20 -uS - | samtools sort - ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${PLASMID}.sorted

	# fastqc -o ${OUTDIR_POOL_QUAL} --contaminants ${SEQCONTAM} -t ${MAXTHREADS} -f bam ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.${PLASMID}.merge.bam
	samtools view ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${PLASMID}.sorted.bam | cut -f1 > ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${PLASMID}.list
done

sort --parallel=5 -u ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.*.list > ${OUTDIR_VECTOR_POOL}/${DISEASE}_${PATIENT}_${POOL}.plasmids.list

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Extract no PLASMIDS reads from raw data"
### ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.fastq.gz
zcat ${R1_FASTQ} | fqreverseextract_pureheader ${OUTDIR_VECTOR_POOL}/${DISEASE}_${PATIENT}_${POOL}.plasmids.list | pigz -f -c > ${TMPDIR}/${BNAME_R1}_noPlasmids.fastq.gz &
zcat ${R2_FASTQ} | fqreverseextract_pureheader ${OUTDIR_VECTOR_POOL}/${DISEASE}_${PATIENT}_${POOL}.plasmids.list | pigz -f -c > ${TMPDIR}/${BNAME_R2}_noPlasmids.fastq.gz &
wait
# var renaming
rm ${R1_FASTQ} ${R2_FASTQ}
R1_FASTQ="${TMPDIR}/${BNAME_R1}_noPlasmids.fastq.gz"
R2_FASTQ="${TMPDIR}/${BNAME_R2}_noPlasmids.fastq.gz"

# file compression
pigz -f ${OUTDIR_VECTOR_POOL}/${DISEASE}_${PATIENT}_${POOL}.plasmids.list
pigz -f ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.*.list
# get the number of remaining sequences
RAW_NO_PLASMID=$((`zcat ${R1_FASTQ} | wc -l | cut -d' ' -f1 `/4)) ;
##### ================ PLASMIDS quantification POOL based - end ======================== #####

# FASTA of RANDOM TAGs
if [ ${FASTQ_QF} = "slim" ]; then
	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Export Random TAGs from R2"
	trimmomatic SE -phred33 -threads $MAXTHREADS ${R2_FASTQ} "/dev/stdout" CROP:12 | fastq_to_fasta -Q 33 -n | pigz -f -c > ${OUTDIR_POOL_QUAL}/r2.${POOL}.qf.noPlasmids.noPhiX.TAGs.fa.gz
fi

## remove first N bases
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Remove first 12 bases"
trimmomatic SE -threads ${MAXTHREADS} -phred33 ${R1_FASTQ} ${TMPDIR}/r1.no12.fastq.gz HEADCROP:12 &
trimmomatic SE -threads ${MAXTHREADS} -phred33 ${R2_FASTQ} ${TMPDIR}/r2.no12.fastq.gz HEADCROP:12 &
wait

rm ${R1_FASTQ} ${R2_FASTQ}

# 1. demux barcode R1
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Demux barcodes of R1 -> LTR only"
fastq-multx -m 1 -B ${BARCODE_LTR} ${TMPDIR}/r1.no12.fastq.gz ${TMPDIR}/r2.no12.fastq.gz -o ${TMPDIR}/r1.no12.%.fastq.gz -o ${TMPDIR}/r2.no12.%.fastq.gz

# 2. demux barcode R2 (revert pairs!!)
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Demux barcodex of R2 -> LC only"
for k in ${BCLTR[@]}; do
  fastq-multx -m 1 -B ${BARCODE_LC} ${TMPDIR}/r2.no12.${k}.fastq.gz ${TMPDIR}/r1.no12.${k}.fastq.gz -o ${TMPDIR}/bcmuxall/r2.no12.${k}.%.fastq.gz -o ${TMPDIR}/bcmuxall/r1.no12.${k}.%.fastq.gz
  rm ${TMPDIR}/r1.no12.${k}.fastq.gz ${TMPDIR}/r2.no12.${k}.fastq.gz ${TMPDIR}/bcmuxall/r1.no12.${k}.unmatched.fastq.gz
done
rm ${TMPDIR}/r1.no12.fastq.gz ${TMPDIR}/r2.no12.fastq.gz ${TMPDIR}/r1.no12.unmatched.fastq.gz ${TMPDIR}/r2.no12.unmatched.fastq.gz

# align, filter and recognize ISs
# !!! WARNING AND TODO !!!
# this approach processes ALL combinations of LTR and LC -> if a specific combination does not exist... it wil do all steps producing empty files... this is NOT CORRECT but (unless some programs return errors) it works.
for b in ${BCLTR[@]}; do
  for k in ${BCLC[@]}; do
    TAG="${b}.${k}"
    # echo $TAG
    for a in ${ASSOBCLIST[@]}; do
    	if [ ${TAG} == ${a} ]; then echo "--- This tag exists in the association file -> analyze it ---";
   #  		##### ================ TRIMMING ======================== #####
   #  		echo $TAG
			# echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Trimming LTR (32 bp)"
			# #1) trim LTR from R1. If trimming this read will remove the whole read, DO NOT rescue it! (policy: everything that is trimmed in R1 BUT will lead to a zero or too small read is removed)
			# flexbar2.5 --reads ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.fastq.gz --reads2 ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.fastq.gz --target ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTR -a ${LTR} -f i1.8 --threads ${MAXTHREADS} -ae LEFT -m 18 -q 1 --max-uncalled 800 -ai -4 -ao 22 -at 3.6 -z GZ
			# #2) trim LTR from R2. If trimming this read will remove the whole read, RESCUE IT (file PREFIX_1_single.fastq, do not keep the PREFIX_2_single.fastq reads because these have been removed by a problem in R1)!
			# flexbar2.5 --reads ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTR_1.fastq.gz --reads2 ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTR_2.fastq.gz --target ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTRR2LTR -a ${LTR_rc} -f i1.8 --threads ${MAXTHREADS} -ae RIGHT -m 2 -q 1 --max-uncalled 800 -ai -4 -ao 22 -at 3.6 -s -z GZ
			# rm ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTR_2.fastq.gz;
			#3) trim LC from R1. If trimming this read will remove the whole read, DO NOT rescue it! (policy: everything that is trimmed in R1 BUT will lead to a zero or too small read is removed)
			#3.1) from the paired-end reads
			flexbar2.5 --reads ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.fastq.gz --reads2 ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.fastq.gz --target ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LC -f i1.8 -a ${LC_rev} --threads ${MAXTHREADS} -ae RIGHT -at 4 -ao 8 -ai -4 -m 2 -q 1 --max-uncalled 800 -z GZ

			#3.2) from the single end rescued reads
			#4) trim LC from R2. If trimming this read will remove the whole read, RESCUE IT (file PREFIX_1_single.fastq, do not keep the PREFIX_2_single.fastq reads because these have been removed by a problem in R1)!
			#flexbar2.5 --reads ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LC_1.fastq.gz --reads2 ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LC_2.fastq.gz --target ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC -f i1.8 -a ${LC_fwd} --threads ${MAXTHREADS} -ae LEFT -at 4 -ao 8 -ai -4 -m 2 -q 1 --max-uncalled 800 -s -z GZ

			flexbar -r ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LC_2.fastq.gz -t ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC -qf i1.8 -n ${MAXTHREADS} -a ${LC_fwd} -at LEFT -ao 14 -z GZ -o --adapter-error-rate 0.1 --adapter-trimmed-out ONLY -u 8000000
			zcat ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC.fastq.gz | fastq_to_fasta -Q 33 -n | fasta2csv | cut -d' ' -f1 > ${TMPDIR}/bcmuxall/r2.${b}.${k}.cleanR1LCR2LC.list
			sort --parallel=5 ${TMPDIR}/bcmuxall/r2.${b}.${k}.cleanR1LCR2LC.list > ${TMPDIR}/bcmuxall/r2.${b}.${k}.cleanR1LCR2LC.sorted.list
			zcat ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LC_1.fastq.gz | fqextract_pureheader <(cat ${TMPDIR}/bcmuxall/r2.${b}.${k}.cleanR1LCR2LC.sorted.list) | pigz --best -f -c > ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC_1.fastq.gz
			mv ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC.fastq.gz ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC_2.fastq.gz
			#5) rescue reads (from step 4 and 2) -> the output reads to align will be the following 3 sets:
			## (1) the paired-end reads correctly trimmed, seq id 0: zcat ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTRR2LTRR1LCR2LC_1|2.fastq.gz
			## (2) the rescued reads from step 2 (removed by R2 too small genomic part) clean by LC, seq id 6: zcat ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTRR2LTRR1LC_rescued.fastq.gz
			## (3) the reacued reads from step 4 (removed by R2 too small genomic part), natively clean, seq id 5: zcat ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTRR2LTRR1LCR2LC_1_single.fastq.gz
			### these are the ID: zcat ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTRR2LTRR1LC_2_single.fastq.gz | fastq_to_fasta -Q 33 | fasta2csv | cut -d' ' -f1 | sort | uniq
			# join the rescued reads (theoretically, no overlaps between these datasets because they are incremental, thus if you discarded a read in a previous step you will miss it in the subsequent, avoiding overlap) -> a simple append is enough
			# cat ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTRR2LTRR1LCR2LC_1_single.fastq.gz ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTRR2LTRR1LC_rescued.fastq.gz > ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTRR2LTRR1LCR2LC_rescued.fastq.gz
			# rm ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTRR2LTRR1LC_2.fastq.gz;

			# # NOW THE INPUT DATASETS ARE:
			# # 1) THE PAIRED-END READS R1 AND R2 CLEANED BY LTR-LC: ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTRR2LTRR1LCR2LC_1.fastq.gz ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTRR2LTRR1LCR2LC_2.fastq.gz
			# # 2) THE RESCUED READS (SINGLE END) CLEANED BY LTR (PARTIALLY BY LC): ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTRR2LTRR1LCR2LC_rescued.fastq.gz
			# # rename these datasets in conventional names:
			# mv ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTRR2LTRR1LCR2LC_1.fastq.gz ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.fastq.gz
			# mv ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTRR2LTRR1LCR2LC_2.fastq.gz ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.noLTRLC.fastq.gz
			# mv ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTRR2LTRR1LCR2LC_rescued.fastq.gz ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.rescued.fastq.gz

			# # ##### ================ PLASMIDS quantification TAG based - start ======================== #####
			# # for PLASMID in kana amp; do
			# # 	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Alignment to reference ${PLASMID} genome"
			# # 	## --- paired-end reads ---
			# # 	# low seed value due to small genome
			# # 	bwa-7.10 mem -k 16 -r 1 -T 15 -c 1 -R "@RG\tID:plasmid\tSM:${PLASMID}\tCN:TIGET" -t ${MAXTHREADS} /opt/genome/vector/plasmid/plasmids.${PLASMID}.fa <(zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.fastq.gz ) <(zcat ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.noLTRLC.fastq.gz ) > ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.${PLASMID}.sam
			# # 	# quality filter of seed-1 value
			# # 	samtools view -F 2308 -q 13 -f 64 -uS ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.${PLASMID}.sam | samtools sort - ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.${PLASMID}.sorted
			# # 	## --- single end (rescued ones) ---
			# # 	bwa-7.10 mem -k 16 -r 1 -T 15 -c 1 -R "@RG\tID:plasmid\tSM:${PLASMID}\tCN:TIGET" -t ${MAXTHREADS} /opt/genome/vector/plasmid/plasmids.${PLASMID}.fa <(zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.rescued.fastq.gz ) > ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.${PLASMID}.single.sam
			# # 	# quality filter of seed-1 value
			# # 	samtools view -F 2308 -q 13 -uS ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.${PLASMID}.single.sam | samtools sort - ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.${PLASMID}.single.sorted
			# # 	## --- merge files ---
			# # 	samtools merge -h ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.${PLASMID}.sam ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.${PLASMID}.merge.bam ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.${PLASMID}.single.sorted.bam ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.${PLASMID}.sorted.bam

			# # 	rm ${TMPDIR}/${TMPDIR}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.${PLASMID}.sam
			# # 	# fastqc -o ${OUTDIR_POOL_QUAL} --contaminants ${SEQCONTAM} -t ${MAXTHREADS} -f bam ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.${PLASMID}.merge.bam
			# # 	samtools view ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.${PLASMID}.merge.bam | cut -f1 > ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.${PLASMID}.list
			# # done
			# # cat ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.*.list | sort | uniq > ${OUTDIR_VECTOR_POOL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.plasmids.list
			# # echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Extract no PLASMIDS reads from raw data"
			# # ### ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.fastq.gz
			# # zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.fastq.gz | fqreverseextract_pureheader ${OUTDIR_VECTOR_POOL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.plasmids.list | gzip > ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.noPlasmids.fastq.gz
			# # zcat ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.noLTRLC.fastq.gz | fqreverseextract_pureheader ${OUTDIR_VECTOR_POOL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.plasmids.list | gzip > ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.noLTRLC.noPlasmids.fastq.gz
			# # zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.rescued.fastq.gz | fqreverseextract_pureheader ${OUTDIR_VECTOR_POOL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.plasmids.list | gzip > ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.noPlasmids.single.fastq.gz



            ##### ========================== REMOVING SEQUENCES THAT DON'T START WHITH FUSION  ================================= #####
			echo "REMOVING SEQUENCES THAT DON'T START WHITH FUSION"
			flexbar -r ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC_1.fastq.gz -b ${FUSION_PRIMERS}  -bu -bt LTAIL -t ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC_1_fusion

			awk 'NR%4==1 {print substr($1,2)}' ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC_1_fusion_barcode_unassigned.fastq > ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC_1_fusion.exclude.list

			zcat ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC_1.fastq.gz  | \
			python2 fqextract.pureheader.v2.py ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC_1_fusion.exclude.list exclude | pigz -f -c > ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.fastq.gz
            zcat ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC_2.fastq.gz  | \
            python2 fqextract.pureheader.v2.py ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC_1_fusion.exclude.list exclude | pigz -f -c > ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.fastq.gz

			rm ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC_1_fusion*

			# # ##### =========================================================================================================== #####






			##### ========================== FASTQ JOIN - start ================================= #####
			fastq-join <(zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.fastq.gz) <(zcat ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.fastq.gz) -o ${OUTDIR_POOL_BCMUXALL}/no12.${b}.${k}.%.fastq
		    # a questo punto i file sono tutti nella forma con suffissi:
			# _join.fastq.gz -> le unite correttamente
			# _un1.fastq.gz -> quelle derivanti da R1 non unite (e similmente _un2.fastq.gz quelle da R2 non unite)
			cp ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.fastq.gz ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.fastq.gz ${OUTDIR_POOL_BCMUXALL}/
			pigz -f ${OUTDIR_POOL_BCMUXALL}/no12.${b}.${k}.join.fastq;
			pigz -f ${OUTDIR_POOL_BCMUXALL}/no12.${b}.${k}.un1.fastq;
			pigz -f ${OUTDIR_POOL_BCMUXALL}/no12.${b}.${k}.un2.fastq;
			cat ${OUTDIR_POOL_BCMUXALL}/no12.${b}.${k}.join.fastq.gz ${OUTDIR_POOL_BCMUXALL}/no12.${b}.${k}.un1.fastq.gz > ${OUTDIR_POOL_BCMUXALL}/no12.${b}.${k}.joinedPlusUn1.fastq.gz
			rm ${OUTDIR_POOL_BCMUXALL}/no12.${b}.${k}.join.fastq.gz ${OUTDIR_POOL_BCMUXALL}/no12.${b}.${k}.un*;
			##### =========================== FASTQ JOIN - end ================================= #####

			##### ================ Vector quantification TAG based - start ======================== #####
		    echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Alignment to reference VECTOR genome"
		    echo -e "R1: ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.fastq.gz \nR2: ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.fastq.gz"
		    ## --- paired-end reads ---
			# low seed value due to small genome
		    bwa mem -k 14 -r 1 -v 1 -T ${BWA_MIN_ALN_LEN} -c 1 -R "@RG\tID:${b}_${k}\tSM:Vector\tCN:TIGET" -t ${MAXTHREADS} ${VECTORGENOME} <(zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.fastq.gz ) | samtools view -F 260 -q 12 -uS - | samtools sort - ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.sorted.md
			samtools index ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.sorted.md.bam ;
			# echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Import data into DB from BAM (new table different from the standard one)"
			# ## --- single end (rescued ones) ---
			# bwa mem -k 14 -r 1 -v 1 -T 15 -c 1 -R "@RG\tID:${b}_${k}\tSM:Vector\tCN:TIGET" -t ${MAXTHREADS} ${VECTORGENOME} <(zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.rescued.fastq.gz ) > ${TMPDIR}/bcmuxall/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.single.sam
			# samtools view -F 2308 -q 14 -uS ${TMPDIR}/bcmuxall/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.single.sam | samtools sort - ${TMPDIR}/bcmuxall/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.single.sorted.md
			# rm ${TMPDIR}/bcmuxall/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.single.sam
			# samtools index ${TMPDIR}/bcmuxall/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.single.sorted.md.bam ;
			# ## --- merge files ---
			# samtools merge -h ${TMPDIR}/bcmuxall/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.sam ${TMPDIR}/bcmuxall/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.merge.bam ${TMPDIR}/bcmuxall/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.single.sorted.md.bam ${TMPDIR}/bcmuxall/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.sorted.md.bam
			## import into DB
			# TAG="${b}.${k}"
			# import_iss --bam ${TMPDIR}/bcmuxall/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.merge.bam --locusRead read1 -p ${MAXTHREADS} -g ${VECTORCIGARGENOMEID} --dbhost ${DBHOSTIP} --dbschema ${DBSCHEMA} --dbtable ${DBTABLE}_lv --dbuser andrea --dbpassword andrea -a ${ASSOCIATIONFILE} --tag ${TAG} --bypassOutFile;
			# extract reads ID of the LTR mapping reads
			samtools view ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.sorted.md.bam | cut -f1 | sort | uniq > ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.sorted.md.vector.list
			# fastqc -o ${OUTDIR_POOL_QUAL} --contaminants ${SEQCONTAM} -t ${MAXTHREADS} -f bam ${TMPDIR}/bcmuxall/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.merge.bam
			echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Extract no LV reads from raw data"
			### ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.fastq.gz
			zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.fastq.gz | fqextract_pureheader ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.sorted.md.vector.list | pigz -f -c > ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.Vector.fastq.gz
			zcat ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.fastq.gz | fqextract_pureheader ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.sorted.md.vector.list | pigz -f -c > ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.Vector.fastq.gz
			# remove
		    # rm ${TMPDIR}/bcmuxall/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.sorted.md.bam
		    pigz -f ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.sorted.md.vector.list
			##### ================ LV quantification TAG based - end ======================== #####

		# else
		# 	echo "Tag not in  the association file, go to the next TAG please.";
		fi
	done
  done
done

#### ======================================= ALIGNMENT PAIRED ENDS ================================================
for b in ${BCLTR[@]}; do
  for k in ${BCLC[@]}; do
  	TAG="${b}.${k}"
    # echo $TAG
    for a in ${ASSOBCLIST[@]}; do
    	if [ ${TAG} == ${a} ]; then echo "--- This tag exists in the association file -> analyze it ---";
    	echo $TAG
		  echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Alignment to reference genome in PE and SE using BWA MEM"
			# test with  -P -M -c 2 -T 15 -> very stringent... too few reads in output
			# explaination: -c 1 -> keep only reads with 1 MEM; -P -> try to rescue some pairs with SW:: ATTENTION!!!!! This option will destroy all paired-end reads!!!!! Do NOT USE IT!; -T 15 -> min alignment score; -M -> flag short alignments (for Picard);
			bwa mem -k 18 -r 1 -M -v 1 -T 15 -c 1 -R "@RG\tID:${TAG}\tSM:${b}.${k}\tCN:TIGET.${DISEASE}.${PATIENT}.${POOL}" -t ${MAXTHREADS} ${GENOME} <(zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.Vector.fastq.gz ) <(zcat ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.Vector.fastq.gz ) > ${TMPDIR}/sam/${b}.${k}.sam

			# create BAM and sort them
			echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Creating BAM and indexes (filter from here the dataset using only valid reads: mapped and primary)"
			samtools view -F 4 -uS ${TMPDIR}/sam/${b}.${k}.sam | samtools sort - ${TMPDIR}/bam/${b}.${k}.F4.sorted.md;
			samtools index ${TMPDIR}/bam/${b}.${k}.F4.sorted.md.bam
			samtools view -F 2308 -uS ${TMPDIR}/sam/${b}.${k}.sam | samtools sort - ${TMPDIR}/bam/${b}.${k}.sorted.md;
			# Removing Reads in chrM and chrUn
			samtools view -h ${TMPDIR}/bam/${b}.${k}.sorted.md.bam | awk '{if($3 != "chrM" && $3 != "chrUn"){print $0}}' | samtools view -Sb - > ${TMPDIR}/bam/${b}.${k}.sorted.md.cleaned.bam
			samtools index ${TMPDIR}/bam/${b}.${k}.sorted.md.cleaned.bam

			# cp ${TMPDIR}/bam/${b}.${k}.sorted.md.cleaned.bam ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.pg.iss.bam

			## Repeats Identification
			if [ ${REPEATS} = "repeats_yes" ]
				then
				# Identifying Low Mapping Quality Reads after hg19 alignment
				bamtools filter -in ${TMPDIR}/bam/${b}.${k}.sorted.md.bam -mapQuality "<=5" > ${TMPDIR}/bam/${b}.${k}.sorted.md.lowMQ.bam;
				# Filter out the R2 from BAM and Converting to FASTA
				bamtools filter -in ${TMPDIR}/bam/${b}.${k}.sorted.md.lowMQ.bam -isFirstMate true > ${TMPDIR}/bam/${b}.${k}.sorted.md.lowMQ.R1.bam;
				bamtools convert -format fasta -in ${TMPDIR}/bam/${b}.${k}.sorted.md.lowMQ.R1.bam -out ${TMPDIR}/${b}.${k}.sorted.md.lowMQ.R1.fa;
				mv ${TMPDIR}/bam/${b}.${k}.sorted.md.lowMQ.R1.bam ${OUTDIR_POOL_ISS_REPEATS};
				# RUN RepeatMasker
				RepeatMasker -no_is -species ${SPECIE} -pa ${MAXTHREADS} -dir ${OUTDIR_POOL_ISS_REPEATS} -q ${TMPDIR}/${b}.${k}.sorted.md.lowMQ.R1.fa;
				rm ${TMPDIR}/${b}.${k}.sorted.md.lowMQ.R1.fa ${OUTDIR_POOL_ISS_REPEATS}/${b}.${k}.sorted.md.lowMQ.R1.fa.cat* ${OUTDIR_POOL_ISS_REPEATS}/*.masked;
				# Import RM BED file into DB, Smith-Waterman >250
				tail -n+4 ${OUTDIR_POOL_ISS_REPEATS}/${b}.${k}.sorted.md.lowMQ.R1.fa.out | awk '{ if ($1 >= 250) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11 }' | sort -uk5,5 | awk 'BEGIN{OFS="\t"}{if($9=="C"){strand="-"}else{strand="+"};print "R_"$11,$6-1,$7,$11"_"$5,".",strand,$10}' > ${OUTDIR_POOL_ISS_REPEATS}/${b}.${k}.RM.R1.bed
				isa_importrediss_frombed -b ${OUTDIR_POOL_ISS_REPEATS}/${b}.${k}.RM.R1.bed -a ${ASSOCIATIONFILE} --patient ${PATIENT} --pool ${POOL} --tag ${TAG} -d ${DBHOSTID} --dbschema ${DBSCHEMA} --dbtable ${DBTABLE}_repeats
				rm ${OUTDIR_POOL_ISS_REPEATS}/${b}.${k}.RM.R1.bed;
			fi

			rm ${TMPDIR}/sam/${b}.${k}.sam ;

		    ##### ======= remove TMP no vector files ==== avoid redundancy =====
		    # rm ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.noVector.fastq.gz
		    # rm ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.noLTRLC.noVector.fastq.gz
		 #   else
			# echo "Tag not in  the association file, go to the next TAG please.";
		fi
	done
  done
done

# #### ======================================= RECALIBRATION ================================================
# echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Clip contaminants (GATK)"
# for b in ${BCLTR[@]}; do
# 	for k in ${BCLC[@]}; do
# 	TAG="${b}.${k}"
# 	# echo $TAG
# 		for a in ${ASSOBCLIST[@]}; do
# 			if [ ${TAG} == ${a} ]; then echo "--- This tag exists in the association file -> analyze it ---";
# 				echo $TAG
# 				## remove/clip eventual residues of internal controls
# 				# gatk-GenomeAnalysisTK -T ClipReads --outputStatistics ${TMPDIR}/bam/${b}.${k}.sorted.md.pg.clipstats -CR SOFTCLIP_BASES -l ERROR -I ${TMPDIR}/bam/${b}.${k}.sorted.md.bam -o ${TMPDIR}/bam/${b}.${k}.sorted.md.pg.bam -R ${GATKREFGENOME} --clipSequencesFile ${CONTAMINANTDB} #--num_threads ${MAXTHREADS}
# 				# samtools index ${TMPDIR}/bam/${b}.${k}.sorted.md.pg.bam
# 				# head -10 ${TMPDIR}/bam/${b}.${k}.sorted.md.pg.clipstats > ${TMPDIR}/bam/${b}.${k}.sorted.md.pg.contaminamnts.clipstats
# 				# rm ${TMPDIR}/bam/${b}.${k}.sorted.md.pg.clipstats

# 				### TMP!!!
# 				ln -s ${TMPDIR}/bam/${b}.${k}.sorted.md.cleaned.bam ${TMPDIR}/bam/${b}.${k}.sorted.md.pg.bam
# 				samtools index ${TMPDIR}/bam/${b}.${k}.sorted.md.pg.bam

# 				#
# 			fi
# 		done
# 	done
# done

#### ======================================= FILTERING ================================================
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filtering data Branching Paired/Single-end mapped reads"
for b in ${BCLTR[@]}; do
  for k in ${BCLC[@]}; do
    TAG="${b}.${k}"
    # echo $TAG
    for a in ${ASSOBCLIST[@]}; do
    	if [ ${TAG} == ${a} ]; then echo "--- This tag exists in the association file -> analyze it ---";
    		echo $TAG

#     		#### ============== BY MATE AND CIGAR-MD TAGS ============== ######
# 		 #    ### ------ CUSTOM FILTERS ---------
# 		 #    echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filtering data BY MATE"
# 		 #    filter_by_mate --bam ${TMPDIR}/bam/${b}.${k}.sorted.md.pg.bam -o ${TMPDIR}/bam/${b}.${k}.sorted.md.toremovebymate.list -p ${MAXTHREADS} -g ${CIGARGENOMEID} --mateoverlap 3 --maxinsertsizethreshold 2 --maxsamapq 15 --suboptimalThreshold ${SUBOPTIMALTHRESHOLD} --deltaMaxSubOptimal 5 --SAalnAction ignore
# 		 #    echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filter by CIGAR and MD tag"
# 			# filter_by_cigar_bam --bam ${TMPDIR}/bam/${b}.${k}.sorted.md.pg.bam -o ${TMPDIR}/bam/${b}.${k}.sorted.md.toremovebycigar.list --minStartingMatches 5 -p ${MAXTHREADS} -g ${CIGARGENOMEID} --minStartingBasesNoIndels 3 --compareSubOptimal --suboptimalThreshold ${SUBOPTIMALTHRESHOLD} --SAalnAction ignore --XSlikeTag XS --ASlikeTag AS --endClipThreshold 12 --pruneByPair read1 --minStartingMatches_fromNbp 2
# 		 # 	# join the two list of reads to remove
# 		 # 	# awk '{print $0}'  is elegant but not reliable because is the fie does not exist, it ends with no output -> better cat
# 		 # 	cat ${TMPDIR}/bam/${b}.${k}.sorted.md.toremovebymate.list ${TMPDIR}/bam/${b}.${k}.sorted.md.toremovebycigar.list | sort | uniq > ${TMPDIR}/bam/${b}.${k}.sorted.md.toremove.list
# 			# CHIMERANROWS=`wc -l ${TMPDIR}/bam/${b}.${k}.sorted.md.toremove.list | cut -d' ' -f1`
# 			# if [ $CHIMERANROWS -gt 0 ]
# 			# 	then
# 			# 	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] --- found reads to remove, now removing them all from BAM file"
# 			# 	picard-tools FilterSamReads INPUT=${TMPDIR}/bam/${b}.${k}.sorted.md.pg.bam FILTER=excludeReadList RLF=${TMPDIR}/bam/${b}.${k}.sorted.md.toremove.list SO=coordinate O=${TMPDIR}/bam/${b}.${k}.sorted.md.filter.bam
# 			# else
# 			# 	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] --- reads to remove not found, go ahead"
# 			# 	cp ${TMPDIR}/bam/${b}.${k}.sorted.md.pg.bam ${TMPDIR}/bam/${b}.${k}.sorted.md.filter.bam
# 			# fi
# 			# samtools index ${TMPDIR}/bam/${b}.${k}.sorted.md.filter.bam


			#### ============== BY ALIGNMENT QUALITY ============== ######
			## questo step deve venire DOPO i miei custom filters se no mi elimina i pairs (poiche' guarda anche la qualita)!
		    ### ------ PAIRED END - PROPER PAIRS ------- BRANCH ---------
			echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filtering data (Bamtools)"
			bamtools filter -in ${TMPDIR}/bam/${b}.${k}.sorted.md.cleaned.bam -isMapped true -isMateMapped true -isPaired true -isProperPair true -isPrimaryAlignment true -mapQuality ">=12" -out ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.bam

			# reindex, otherwise GATK does not work...!
		    samtools index ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.bam
		    cp ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.bam ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.pg.iss.bam
		    samtools index ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.pg.iss.bam

			## create a BED file wiht only R2 reads -> this output file will be used to extract the product end (that will be the R2 start)
			echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Creating BED file of ONLY R2 reads (this file will be used to extract R2 start as product end)"
			bamtools filter -in ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.bam -isSecondMate true -isFirstMate false -isMapped true -isMateMapped true -isPaired true -isProperPair true -isPrimaryAlignment true | bedtools bamtobed > ${TMPDIR}/bed/${b}.${k}.sorted.allr2reads.bed
			rm ${TMPDIR}/bam/${b}.${k}*.reads;
		fi
	done
  done
done

#### ======================================= FOLDER ORGANIZATION ================================================
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Create BED files"
for b in ${BCLTR[@]}; do
  for k in ${BCLC[@]}; do
    TAG="${b}.${k}"
    # echo $TAG
    for a in ${ASSOBCLIST[@]}; do
    	if [ ${TAG} == ${a} ]; then echo "--- This tag exists in the association file -> analyze it ---";
			echo $TAG
			bedtools bamtobed -cigar -i ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.pg.iss.bam | grep '/1' > ${TMPDIR}/bed/${b}.${k}.sorted.md.rel.pg.bed
		# else
		# 	echo "Tag not in  the association file, go to the next TAG please.";
		fi
	done
  done
done

#### ======================================= DB IMPORT AND FINAL ISS BAM GENERATION ================================================
for b in ${BCLTR[@]}; do
	for k in ${BCLC[@]}; do
	TAG="${b}.${k}"
	# echo $TAG
		for a in ${ASSOBCLIST[@]}; do
			if [ ${TAG} == ${a} ]; then echo "--- This tag exists in the association file -> analyze it ---";
				echo $TAG
				# sort BEd file
				bedtools sort -i ${TMPDIR}/bed/${b}.${k}.sorted.md.rel.pg.bed > ${TMPDIR}/bed/${b}.${k}.sorted.md.rel.pg.sorted.bed

				### TODO: script to update product end (not read end!!) => product end is the R2 start. The program that we need to develop updates the end field of the bed file with suffix *.cigarmd.bed. Then the program if DB import must be updated!

				# echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Import BED data into DB (my script - running only for MY USER! so far)"
				TAG="${b}.${k}"
				#isa_importrediss_frombed -b ${TMPDIR}/bed/${b}.${k}.sorted.md.rel.pg.bed -a ${ASSOCIATIONFILE} --patient ${PATIENT} --pool ${POOL} --tag ${TAG} -d ${DBHOSTID} --dbschema ${DBSCHEMA} --dbtable ${DBTABLE}

				echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Create the BAM file of this final output with ONLY the valid ISs "
				# extract unique headers from the BED file
				cut -f 4 ${TMPDIR}/bed/${b}.${k}.sorted.md.rel.pg.bed | sed 's/\/1//g' | sed 's/\/2//g' | sort | uniq > ${TMPDIR}/bed/${b}.${k}.sorted.md.rel.pg.sorted.list
				BEDNROWS=`wc -l ${TMPDIR}/bed/${b}.${k}.sorted.md.rel.pg.sorted.list | cut -d' ' -f1`
				# if ROWS > 0 => do the filtering, else the BAM is already at the final version (: it means that no simple repeats are in the file)
				if [ $BEDNROWS -gt 0 ]
				then
					#extract the ISs of the BED file into a BAM file
					picard-tools FilterSamReads INPUT=${TMPDIR}/bam/${b}.${k}.sorted.md.rel.bam FILTER=includeReadList RLF=${TMPDIR}/bed/${b}.${k}.sorted.md.rel.pg.sorted.list SO=coordinate O=${TMPDIR}/bam/${b}.${k}.sorted.md.rel.pg.iss.bam

				else
					cp ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.bam ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.pg.iss.bam

				fi
				rm ${TMPDIR}/bam/${b}.${k}*.reads;
                filter_by_cigar_bam_noprint_onlyrep_stats --bam ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.pg.iss.bam\
                  -o ${TMPDIR}/bam/${b}.${k}.toremovebycigar.list --minStartingMatches ${minSM_rp} -p ${MAXTHREADS} -g ${CIGARGENOMEID}\
                       --minStartingBasesNoIndels 50 --compareSubOptimal --suboptimalThreshold ${SUBOPTIMALTHRESHOLD_rp}\
                        --SAalnAction ignore --XSlikeTag XS --ASlikeTag AS --endClipThreshold 12 --pruneByPair read1\
                         --minStartingMatches_fromNbp ${minSM_fromN_rp}
                if [ ! -s ${TMPDIR}/bam/${b}.${k}.toremovebycigar.list ]
                then
	                echo "No reads to be filtered, going on"
                else
    	            picard-tools FilterSamReads I=${TMPDIR}/bam/${b}.${k}.sorted.md.rel.pg.iss.bam\
    	             O=${TMPDIR}/bam/${b}.${k}.sorted.md.rel.pg.iss.filtered.bam FILTER=excludeReadList\
    	              RLF=${TMPDIR}/bam/${b}.${k}.toremovebycigar.list
    	            mv ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.pg.iss.filtered.bam ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.pg.iss.bam
                fi
                samtools index ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.pg.iss.bam
                rm ${TMPDIR}/bam/${b}.${k}.toremovebycigar.list.;
		#bamPEFragmentSize -b ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.pg.iss.bam -hist ${OUTDIR_POOL_QUAL}/FragmentSize/${b}.${k}.sorted.md.rel.pg.iss.histLengths.png -p ${MAXTHREADS} -T ${TAG} --samplesLabel ${TAG}

				### import data using the new DB structure
				echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Import data into DB from BAM (new table different from the standard one)"
				TAG="${b}.${k}"
				import_iss --bam ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.pg.iss.bam --locusRead read1 -p ${MAXTHREADS} -g ${CIGARGENOMEID} --dbhost ${DBHOSTIP} --dbschema ${DBSCHEMA} --dbtable ${DBTABLE}_refactored --dbuser andrea --dbpassword andrea -a ${ASSOCIATIONFILE} --tag ${TAG} --outfile ${OUTDIR_POOL_ISS}/${DBSCHEMA}_${DBTABLE}_refactored.tsv --bypassDB;
				import_iss --bam ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.pg.iss.bam --locusRead read2 -p ${MAXTHREADS} -g ${CIGARGENOMEID} --dbhost ${DBHOSTIP} --dbschema ${DBSCHEMA} --dbtable ${DBTABLE}_refactored --dbuser andrea --dbpassword andrea -a ${ASSOCIATIONFILE} --tag ${TAG} --outfile ${OUTDIR_POOL_ISS}/${DBSCHEMA}_${DBTABLE}_refactored.read2.tsv --bypassDB;
				rm ${TMPDIR}/bam/${b}.${k}*.reads;
				rm ${TMPDIR}/bed/${b}.${k}.sorted.md.rel.pg.sorted.list
			# else
			# 	echo "Tag not in  the association file, go to the next TAG please.";
			fi
		done
	done
 done
#### ======================================= REDO USING REVERSE COMPLEMENT ================================================
#for b in ${BCLTR[@]}; do
#	for k in ${BCLC[@]}; do
#	TAG="${b}.${k}"
#	# echo $TAG
#		for a in ${ASSOBCLIST[@]}; do
#			if [ ${TAG} == ${a} ]; then echo "--- This tag exists in the association file -> analyze it ---";
#				echo $TAG
#				echo "Reverse complement"
#				# extracting reads name for which we have an integration
#				samtools view ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.pg.iss.bam | cut -f1 | sort | uniq >\
#				         ${TMPDIR}/bam/${b}.${k}.reads_to_exclude.list
#
#				# excluding those reads from succcessive analysis. you can pass a second parameter to fqextract.v2
#				# if done, no matter the value, it will exclude reads in the list
#				# excluding and make the reverse complement now r1 becomes r2 and viceversa
#			    zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.fastq.gz  | python2 /home/carlo/workspace/filtered/fqextract.pureheader.v2.py \
#			        ${TMPDIR}/bam/${b}.${k}.reads_to_exclude.list exclude | fastx_reverse_complement -z \
#			         -o ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.revComp.fastq.gz
#
#               zcat ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.fastq.gz  | python2 /home/carlo/workspace/filtered/fqextract.pureheader.v2.py \
#			        ${TMPDIR}/bam/${b}.${k}.reads_to_exclude.list exclude | fastx_reverse_complement -z \
#			         -o ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.revComp.fastq.gz
#
#			    bwa mem -k 14 -r 1 -v 1 -T ${BWA_MIN_ALN_LEN} -c 1 -R "@RG\tID:${b}_${k}\tSM:Vector\tCN:TIGET" -t ${MAXTHREADS} \
#			     ${VECTORGENOME} <(zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.revComp.fastq.gz ) | samtools view -F 260 -q 12 -uS - | \
#			        samtools sort - ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.revComp.noLTRLC.sorted.md
#			    samtools index ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.revComp.noLTRLC.sorted.md.bam ;
#
#			    samtools view ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.revComp.noLTRLC.sorted.md.bam | cut -f1 | sort | #uniq > ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.revComp.sorted.md.vector.list
#
#			    echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Extract no LV reads from raw data"
#			    ### ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.fastq.gz
#			    zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.revComp.fastq.gz | fqextract_pureheader ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.revComp.sorted.md.vector.list | pigz -f -c > ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.revComp.Vector.fastq.gz
#			    zcat ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.revComp.fastq.gz | fqextract_pureheader ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.revComp.sorted.md.vector.list | pigz -f -c > ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.revComp.Vector.fastq.gz
#			    # remove
#		        # rm ${TMPDIR}/bcmuxall/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.sorted.md.bam
#		        pigz -f ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.revComp.sorted.md.vector.list
#
#			    echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Alignment to reference genome in PE and SE using BWA MEM"
#			    # test with  -P -M -c 2 -T 15 -> very stringent... too few reads in output
#			    # explaination: -c 1 -> keep only reads with 1 MEM; -P -> try to rescue some pairs with SW:: ATTENTION!!!!! This option will d#estroy all paired-end reads!!!!! Do NOT USE IT!; -T 15 -> min alignment score; -M -> flag short alignments (for Picard);
	#		    bwa mem -k 18 -r 1 -M -v 1 -T 15 -c 1 -R "@RG\tID:$j\tSM:${b}.${k}\tCN:TIGET.${DISEASE}.${PATIENT}.${POOL}" -t ${MAXTHREADS} ${GENOME} <(zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.revComp.Vector.fastq.gz ) <(zcat ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.revComp.Vector.fastq.gz ) > ${TMPDIR}/sam/${b}.${k}.revComp.sam
#
#			    # create BAM and sort them
#			    echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Creating BAM and indexes (filter from here the dataset using only valid reads: mapped and primary)"
#			    samtools view -F 4 -uS ${TMPDIR}/sam/${b}.${k}.revComp.sam | samtools sort - ${TMPDIR}/bam/${b}.${k}.revComp.F4.sorted.md;
#			    samtools index ${TMPDIR}/bam/${b}.${k}.revComp.F4.sorted.md.bam
#			    samtools view -F 2308 -uS ${TMPDIR}/sam/${b}.${k}.revComp.sam | samtools sort - ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md;
#			    # Removing Reads in chrM and chrUn
#			    samtools view -h ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.bam | awk '{if($3 != "chrM" && $3 != "chrUn"){print $0}}' | samtools view -Sb - > ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.cleaned.bam
#			    samtools index ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.cleaned.bam
#                rm ${TMPDIR}/sam/${b}.${k}.revComp.sam ;
#
#
#
#               echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filtering data (Bamtools)"
#			    bamtools filter -in ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.cleaned.bam -isMapped true -isMateMapped true -isPaired true -isProperPair true -isPrimaryAlignment true -mapQuality ">=12" -out ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.rel.bam
#
#			    # reindex, otherwise GATK does not work...!
#		        samtools index ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.rel.bam
#		        cp ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.rel.bam ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.rel.pg.iss.bam
#		        samtools index ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.rel.pg.iss.bam
#
#			    ## create a BED file wiht only R2 reads -> this output file will be used to extract the product end (that will be the R2 start)
#			    echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Creating BED file of ONLY R2 reads (this file will be used to extract R2 start as product end)"
#			    bamtools filter -in ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.rel.bam -isSecondMate true -isFirstMate false -isMapped \
#			     true -isMateMapped true -isPaired true -isProperPair true -isPrimaryAlignment true | \
#			        bedtools bamtobed > ${TMPDIR}/bed/${b}.${k}.revComp.sorted.allr2reads.bed
#			    rm ${TMPDIR}/bam/${b}.${k}*.reads;
#
#
#			    bedtools bamtobed -cigar -i ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.rel.pg.iss.bam | \
#			        grep '/1' > ${TMPDIR}/bed/${b}.${k}.revComp.sorted.md.rel.pg.bed
#
#			    # sort BEd file
#				bedtools sort -i ${TMPDIR}/bed/${b}.${k}.revComp.sorted.md.rel.pg.bed > \
#				    ${TMPDIR}/bed/${b}.${k}.revComp.sorted.md.rel.pg.sorted.bed
#
#				### TODO: script to update product end (not read end!!) => product end is the R2 start. The program that we need to develop #updates the end field of the bed file with suffix *.cigarmd.bed. Then the program if DB import must be updated!
#
#				# echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Import BED data into DB (my script - running only for MY USER! so far)"
#				TAG="${b}.${k}"
#
#				echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Create the BAM file of this final output with ONLY the valid ISs "
#				# extract unique headers from the BED file
#				cut -f 4 ${TMPDIR}/bed/${b}.${k}.revComp.sorted.md.rel.pg.bed | sed 's/\/1//g' | sed 's/\/2//g' | sort | uniq > ${TMPDIR}/bed/${b}.${k}.revComp.sorted.md.rel.pg.sorted.list
#				BEDNROWS=`wc -l ${TMPDIR}/bed/${b}.${k}.revComp.sorted.md.rel.pg.sorted.list | cut -d' ' -f1`
#				# if ROWS > 0 => do the filtering, else the BAM is already at the final version (: it means that no simple repeats are in the file)
#				if [ $BEDNROWS -gt 0 ]
#				then
#					#extract the ISs of the BED file into a BAM file
#					picard-tools FilterSamReads INPUT=${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.rel.bam FILTER=includeReadList RLF=${TMPDIR}/bed/${b}.${k}.revComp.sorted.md.rel.pg.sorted.list SO=coordinate O=${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.rel.pg.iss.bam
#
#				else
#					cp ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.rel.bam ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.rel.pg.iss.bam
#
#				fi
#				rm ${TMPDIR}/bam/${b}.${k}*.reads;
 #               /home/carlo/workspace/filtered/filter_by_cigar_bam_noprint_onlyrep_stats --bam \
  #               ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.rel.pg.iss.bam -o ${TMPDIR}/bam/${b}.${k}.revComp.toremovebycigar.list \
   #              --minStartingMatches ${minSM_rp} -p ${MAXTHREADS} -g ${CIGARGENOMEID}\
#                       --minStartingBasesNoIndels 50 --compareSubOptimal --suboptimalThreshold ${SUBOPTIMALTHRESHOLD_rp}\
#                        --SAalnAction ignore --XSlikeTag XS --ASlikeTag AS --endClipThreshold 12 --pruneByPair read1\
#                         --minStartingMatches_fromNbp ${minSM_fromN_rp}
#                if [ ! -s ${TMPDIR}/bam/${b}.${k}.revComp.toremovebycigar.list ]
#                then
#	                echo "No reads to be filtered, going on"
#                else
#   	            picard-tools FilterSamReads I=${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.rel.pg.iss.bam\
#    	             O=${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.rel.pg.iss.filtered.bam FILTER=excludeReadList\
#    	              RLF=${TMPDIR}/bam/${b}.${k}.revComp.toremovebycigar.list
#    	            mv ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.rel.pg.iss.filtered.bam ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.rel.pg.iss.bam
#                fi
#               samtools index ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.rel.pg.iss.bam
#                rm ${TMPDIR}/bam/${b}.${k}.revComp.toremovebycigar.list.;
		#bamPEFragmentSize -b ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.pg.iss.bam -hist ${OUTDIR_POOL_QUAL}/FragmentSize/${b}.${k}.sorted.md.rel.pg.iss.histLengths.png -p ${MAXTHREADS} -T ${TAG} --samplesLabel ${TAG}
#
#				### import data using the new DB structure
#				echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Import data into DB from BAM (new table different from the standard one)"
#				TAG="${b}.${k}"
#				import_iss --bam ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.rel.pg.iss.bam --locusRead read1 -p ${MAXTHREADS} -g ${CIGARGENOMEID} --dbhost ${DBHOSTIP} --dbschema ${DBSCHEMA} --dbtable ${DBTABLE}_refactored --dbuser andrea --dbpassword andrea -a ${ASSOCIATIONFILE} --tag ${TAG} --outfile ${OUTDIR_POOL_ISS}/${DBSCHEMA}_${DBTABLE}_refactored.tsv --bypassDB;
#				import_iss --bam ${TMPDIR}/bam/${b}.${k}.revComp.sorted.md.rel.pg.iss.bam --locusRead read2 -p ${MAXTHREADS} -g ${CIGARGENOMEID} --dbhost ${DBHOSTIP} --dbschema ${DBSCHEMA} --dbtable ${DBTABLE}_refactored --dbuser andrea --dbpassword andrea -a ${ASSOCIATIONFILE} --tag ${TAG} --outfile ${OUTDIR_POOL_ISS}/${DBSCHEMA}_${DBTABLE}_refactored.read2.tsv --bypassDB;
#				rm ${TMPDIR}/bam/${b}.${k}*.reads;
#
 #           fi
	#	done
	#done
 #done




#### ======================================= STATS ================================================
for b in ${BCLTR[@]}; do
  for k in ${BCLC[@]}; do
    TAG="${b}.${k}"
    # echo $TAG
    for a in ${ASSOBCLIST[@]}; do
    	if [ ${TAG} == ${a} ]; then echo "--- This tag exists in the association file -> analyze it ---";
    		echo $TAG
			# checking VARS
			echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Checking Variables ";
			LV_MAPPING_RESULTS=`samtools flagstat ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.sorted.md.bam ` ; # ${TMPDIR}/bcmuxall/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.noLTRLC.sorted.md.bam
			LV_MAPPING_PP=$((`echo ${LV_MAPPING_RESULTS} | grep ' properly paired ' | cut -d' ' -f34`/2)) ;
			LV_MAPPING_ST=$((`echo ${LV_MAPPING_RESULTS} | grep ' singletons ' | cut -d' ' -f48`)) ;
			LV_MAPPING_OVERALL=$((${LV_MAPPING_ST}+${LV_MAPPING_PP})) ;
			LV_MAPPING_MAPPED=$((`echo ${LV_MAPPING_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;
			LV_MAPPED=`zcat ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${b}.${k}.sorted.md.vector.list.gz | wc -l | cut -d' ' -f1 `

			PLASMID_MAPPED_BYPOOL=`zcat ${OUTDIR_VECTOR_POOL}/${DISEASE}_${PATIENT}_${POOL}.plasmids.list.gz | wc -l | cut -d' ' -f1 `

			BARCODE_MUX=$((`zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.fastq.gz | wc -l `/4)) ;
			LTR_IDENTIFIED="0"; #$((`zcat ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTR_1.fastq.gz | wc -l | cut -d' ' -f1 `/4)) ;
			LTR_SEQAFTERTRIMMING="0"; #$((`zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.fastq.gz | wc -l | cut -d' ' -f1 `/4)) ;
			LC_SEQAFTERTRIMMING="0"; #$((`zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.fastq.gz | wc -l | cut -d' ' -f1 `/4)) ;
			TRIMMING_LTRR1="0"; #$((`zcat ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTR_1.fastq.gz | wc -l | cut -d' ' -f1 `/4)) ;
			TRIMMING_LTRR1R2="0"; #$((`zcat ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTRR2LTR_1.fastq.gz | wc -l | cut -d' ' -f1 `/4)) ;
			TRIMMING_LTRR1R2_LCR1="0"; #$((`zcat ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LTRR2LTRR1LC_1.fastq.gz | wc -l | cut -d' ' -f1 `/4)) ;
			TRIMMING_FINAL_RESCUED="0"; #$((`zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.rescued.fastq.gz | wc -l | cut -d' ' -f1 `/4)) ;
			TRIMMING_FINAL_LTRLC="0"; #$((`zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.fastq.gz | wc -l | cut -d' ' -f1 `/4)) ;

			BWA_RESULTS=`samtools flagstat ${TMPDIR}/bam/${b}.${k}.sorted.md.cleaned.bam `
			BWA_INPUT_PAIRED=$((`zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.Vector.fastq.gz  | wc -l | cut -d' ' -f1 `/4)) ;
			BWA_INPUT_SINGLE="0"; #$((`zcat ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.noVector.single.fastq.gz  | wc -l | cut -d' ' -f1 `/4)) ;
			BWA_INPUT=$((${BWA_INPUT_PAIRED})) ;
			BWA_MAPPED=$((`echo ${BWA_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;
			BWA_MAPPED_PP=$((`echo ${BWA_RESULTS} | grep ' properly paired ' | cut -d' ' -f34`/2)) ;
			BWA_MAPPED_ST=$((`echo ${BWA_RESULTS} | grep ' singletons (' | cut -d' ' -f48`)) ;
			BWA_MAPPED_OVERALL=$((${BWA_MAPPED_PP})) ;
			BWA_ALIGNED_R1=$((`echo ${BWA_RESULTS} | grep ' read1' | cut -d' ' -f 26`)) ;

			# RECALIB_RESULTS=`samtools flagstat ${TMPDIR}/bam/${b}.${k}.sorted.md.pg.bam `
			# RECALIB_MAPPED=$((`echo ${RECALIB_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;
			# RECALIB_MAPPED_PP=$((`echo ${RECALIB_RESULTS} | grep ' properly paired ' | cut -d' ' -f34`/2)) ;
			# RECALIB_MAPPED_ST=$((`echo ${RECALIB_RESULTS} | grep ' singletons (' | cut -d' ' -f48`)) ;
			# RECALIB_MAPPED_OVERALL=$((${RECALIB_MAPPED_PP}+${RECALIB_MAPPED_ST})) ;
			# RECALIB_ALIGNED_R1=$((`echo ${RECALIB_RESULTS} | grep ' read1' | cut -d' ' -f 26`)) ;
			# RECALIB_SOFCLIPPED_READS=` grep 'Number of clipped reads' ${TMPDIR}/bam/${b}.${k}.sorted.md.pg.contaminamnts.clipstats | cut -d' ' -f6-100 | sed 's/ //g' ` # number of soft clipped reads

			RECALIB_RESULTS="0"
			RECALIB_MAPPED="0"
			RECALIB_MAPPED_PP="0"
			RECALIB_MAPPED_ST="0"
			RECALIB_MAPPED_OVERALL="0"
			RECALIB_ALIGNED_R1="0"
			RECALIB_SOFCLIPPED_READS="0"

			FILTER_MATE_TO_REMOVE="0"; #`wc -l ${TMPDIR}/bam/${b}.${k}.sorted.md.toremovebymate.list | cut -d' ' -f1 `
			FILTER_CIGARMD_TO_REMOVE="0"; #`wc -l ${TMPDIR}/bam/${b}.${k}.sorted.md.toremovebycigar.list | cut -d' ' -f1 `
			FILTER_JOINT_MC_TO_REMOVE="0"; #`wc -l ${TMPDIR}/bam/${b}.${k}.sorted.md.toremove.list | cut -d' ' -f1 `
			FILTER_JOINT_MC_RESULTS="0"; #`samtools flagstat ${TMPDIR}/bam/${b}.${k}.sorted.md.filter.bam ` ;
			FILTER_JOINT_MC_PP="0"; #$((`echo ${FILTER_JOINT_MC_RESULTS} | grep ' properly paired ' | cut -d' ' -f34`/2)) ;
			FILTER_JOINT_MC_ST="0"; #$((`echo ${FILTER_JOINT_MC_RESULTS} | grep ' singletons ' | cut -d' ' -f48`)) ;
			FILTER_JOINT_MC_OVERALL="0"; #$((${FILTER_JOINT_MC_ST}+${FILTER_JOINT_MC_PP})) ;
			FILTER_JOINT_ALIGNED_R1="0"; #$((`echo ${FILTER_JOINT_MC_RESULTS} | grep ' read1' | cut -d' ' -f 26`)) ;

			FILTER_ALMQUAL_RESULTS=`samtools flagstat ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.bam ` ;
			FILTER_ALMQUAL_PP=$((`echo ${FILTER_ALMQUAL_RESULTS} | grep ' properly paired ' | cut -d' ' -f34`/2)) ;
			FILTER_ALMQUAL_ST=$((`echo ${FILTER_ALMQUAL_RESULTS} | grep ' singletons ' | cut -d' ' -f48`)) ;
			FILTER_ALMQUAL_OVERALL=$((${FILTER_ALMQUAL_ST}+${FILTER_ALMQUAL_PP})) ;
			FILTER_ALMQUAL_ALIGNED_R1=$((`echo ${FILTER_ALMQUAL_RESULTS} | grep ' read1' | cut -d' ' -f 26`)) ;

			ISS_RESULTS=`samtools flagstat ${TMPDIR}/bam/${b}.${k}.sorted.md.rel.pg.iss.bam `
			ISS_MAPPED=$((`echo ${ISS_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;
			ISS_MAPPED_PP=$((`echo ${ISS_RESULTS} | grep ' properly paired ' | cut -d' ' -f34`/2)) ;
			ISS_MAPPED_ST=$((`echo ${ISS_RESULTS} | grep ' singletons (' | cut -d' ' -f48`)) ;
			ISS_MAPPED_OVERALL=$((${ISS_MAPPED_PP}+${ISS_MAPPED_ST})) ;
			ISS_ALIGNED_R1=$((`echo ${ISS_RESULTS} | grep ' read1' | cut -d' ' -f 26`)) ;

			ISS_FINAL=`wc -l ${TMPDIR}/bed/${b}.${k}.sorted.md.rel.pg.bed | cut -d' ' -f1 ` ;

			LTR_ID="${b}"
			LC_ID="${k}"

			############## DB STATS - insert - start #########################################
			echo "
			::VARIABLE SUMMARY:: dlimiters<>
				RUN_ID=<${RUN_ID}>
				RUN_NAME=<${RUN_NAME}>
				DISEASE=<${DISEASE}>
				PATIENT=<${PATIENT}>
				POOL=<${POOL}>
				TAG=<${TAG}>
				LTR_ID=<${LTR_ID}>
				LC_ID=<${LC_ID}>
				RAW_READS=<${RAW_READS}>
				QUALITY_PASSED=<${QUALITY_PASSED}>
				PHIX_MAPPING=<${PHIX_MAPPING}>
				BARCODE_MUX=<${BARCODE_MUX}>
				LTR_IDENTIFIED=<${LTR_IDENTIFIED}>
				TRIMMING_LTRR1=<${TRIMMING_LTRR1}>
				TRIMMING_LTRR1R2=<${TRIMMING_LTRR1R2}>
				TRIMMING_LTRR1R2_LCR1=<${TRIMMING_LTRR1R2_LCR1}>
				TRIMMING_FINAL_RESCUED=<${TRIMMING_FINAL_RESCUED}>
				TRIMMING_FINAL_LTRLC=<${TRIMMING_FINAL_LTRLC}>
				LV_MAPPED=<${LV_MAPPED}>
				PLASMID_MAPPED_BYPOOL=<${PLASMID_MAPPED_BYPOOL}>
				RAW_NO_PLASMID=<${RAW_NO_PLASMID}>
				BWA_INPUT=<${BWA_INPUT}>
				BWA_MAPPED=<${BWA_MAPPED}>
				BWA_MAPPED_PP=<${BWA_MAPPED_PP}>
				BWA_MAPPED_ST=<${BWA_MAPPED_ST}>
				BWA_MAPPED_OVERALL=<${BWA_MAPPED_OVERALL}>
				BWA_ALIGNED_R1=<${BWA_ALIGNED_R1}>
				RECALIB_MAPPED=<${RECALIB_MAPPED}>
				RECALIB_MAPPED_PP=<${RECALIB_MAPPED_PP}>
				RECALIB_MAPPED_ST=<${RECALIB_MAPPED_ST}>
				RECALIB_MAPPED_OVERALL=<${RECALIB_MAPPED_OVERALL}>
				RECALIB_ALIGNED_R1=<${RECALIB_ALIGNED_R1}>
				RECALIB_SOFCLIPPED_READS=<${RECALIB_SOFCLIPPED_READS}>
				FILTER_MATE_TO_REMOVE=<${FILTER_MATE_TO_REMOVE}>
				FILTER_CIGARMD_TO_REMOVE=<${FILTER_CIGARMD_TO_REMOVE}>
				FILTER_JOINT_MC_TO_REMOVE=<${FILTER_JOINT_MC_TO_REMOVE}>
				FILTER_JOINT_MC_PP=<${FILTER_JOINT_MC_PP}>
				FILTER_JOINT_MC_ST=<${FILTER_JOINT_MC_ST}>
				FILTER_JOINT_MC_OVERALL=<${FILTER_JOINT_MC_OVERALL}>
				FILTER_JOINT_ALIGNED_R1=<${FILTER_JOINT_ALIGNED_R1}>
				FILTER_ALMQUAL_PP=<${FILTER_ALMQUAL_PP}>
				FILTER_ALMQUAL_ST=<${FILTER_ALMQUAL_ST}>
				FILTER_ALMQUAL_OVERALL=<${FILTER_ALMQUAL_OVERALL}>
				FILTER_ALMQUAL_ALIGNED_R1=<${FILTER_ALMQUAL_ALIGNED_R1}>
				ISS_FINAL=<${ISS_FINAL}>
				ISS_MAPPED=<${ISS_MAPPED}>
				ISS_MAPPED_PP=<${ISS_MAPPED_PP}>
				ISS_MAPPED_ST=<${ISS_MAPPED_ST}>
				ISS_MAPPED_OVERALL=<${ISS_MAPPED_OVERALL}>
				ISS_ALIGNED_R1=<${ISS_ALIGNED_R1}>
			" ;

		   	echo "[TIGET] Import STATS into SUMMARY table"
		   	TAG="${b}.${k}" ;
			mysql -u andrea --password=andrea -e "INSERT INTO ${DBSCHEMA}.${STATSDBNAME}
			(RUN_ID, RUN_NAME, DISEASE, PATIENT, POOL, TAG, LTR_ID, LC_ID, RAW_READS, QUALITY_PASSED, PHIX_MAPPING, PLASMID_MAPPED_BYPOOL, RAW_NO_PLASMID, BARCODE_MUX, LTR_IDENTIFIED, TRIMMING_LTRR1, TRIMMING_LTRR1R2, TRIMMING_LTRR1R2_LCR1, TRIMMING_FINAL_RESCUED, TRIMMING_FINAL_LTRLC, LV_MAPPED, BWA_INPUT, BWA_MAPPED, BWA_MAPPED_PP, BWA_MAPPED_ST, BWA_MAPPED_OVERALL, BWA_ALIGNED_R1, RECALIB_MAPPED, RECALIB_MAPPED_PP, RECALIB_MAPPED_ST, RECALIB_MAPPED_OVERALL, RECALIB_ALIGNED_R1, RECALIB_SOFCLIPPED_READS, FILTER_MATE_TO_REMOVE, FILTER_CIGARMD_TO_REMOVE, FILTER_JOINT_MC_TO_REMOVE, FILTER_JOINT_MC_PP, FILTER_JOINT_MC_ST, FILTER_JOINT_MC_OVERALL, FILTER_JOINT_ALIGNED_R1, FILTER_ALMQUAL_PP, FILTER_ALMQUAL_ST, FILTER_ALMQUAL_OVERALL, FILTER_ALMQUAL_ALIGNED_R1, ISS_FINAL, ISS_MAPPED, ISS_MAPPED_PP, ISS_MAPPED_ST, ISS_MAPPED_OVERALL, ISS_ALIGNED_R1)
			VALUES
			('${RUN_ID}', '${RUN_NAME}', '${DISEASE}', '${PATIENT}', '${POOL}', '${TAG}', '${LTR_ID}', '${LC_ID}', '${RAW_READS}', '${QUALITY_PASSED}', '${PHIX_MAPPING}', '${PLASMID_MAPPED_BYPOOL}', '${RAW_NO_PLASMID}', '${BARCODE_MUX}', '${LTR_IDENTIFIED}', '${TRIMMING_LTRR1}', '${TRIMMING_LTRR1R2}', '${TRIMMING_LTRR1R2_LCR1}', '${TRIMMING_FINAL_RESCUED}', '${TRIMMING_FINAL_LTRLC}', '${LV_MAPPED}', '${BWA_INPUT}', '${BWA_MAPPED}', '${BWA_MAPPED_PP}', '${BWA_MAPPED_ST}', '${BWA_MAPPED_OVERALL}', '${BWA_ALIGNED_R1}', '${RECALIB_MAPPED}', '${RECALIB_MAPPED_PP}', '${RECALIB_MAPPED_ST}', '${RECALIB_MAPPED_OVERALL}', '${RECALIB_ALIGNED_R1}', '${RECALIB_SOFCLIPPED_READS}', '${FILTER_MATE_TO_REMOVE}', '${FILTER_CIGARMD_TO_REMOVE}', '${FILTER_JOINT_MC_TO_REMOVE}', '${FILTER_JOINT_MC_PP}', '${FILTER_JOINT_MC_ST}', '${FILTER_JOINT_MC_OVERALL}', '${FILTER_JOINT_ALIGNED_R1}', '${FILTER_ALMQUAL_PP}', '${FILTER_ALMQUAL_ST}', '${FILTER_ALMQUAL_OVERALL}', '${FILTER_ALMQUAL_ALIGNED_R1}', '${ISS_FINAL}', '${ISS_MAPPED}', '${ISS_MAPPED_PP}', '${ISS_MAPPED_ST}', '${ISS_MAPPED_OVERALL}', '${ISS_ALIGNED_R1}')
			" ${DBSCHEMA} ;
			############## DB STATS - insert - end #########################################
		# else
		# 	echo "Tag not in  the association file, go to the next TAG please.";
		fi
	done
  done;
done

# add comment to DB
NOWIS=`date +'%y-%m-%d %H:%M:%S'`
mysql -u andrea --password=andrea -e "ALTER TABLE ${DBSCHEMA}.${DBTABLE} COMMENT = 'Pipe: v.${PIPEVERSION}, Author: AC GS, from ${STARTTIME} to ${NOWIS}' ; " ;

# RUN Stats
mysql -u andrea --password=andrea --column-names=true -e "SELECT RUN_NAME, POOL, TAG, RAW_READS, QUALITY_PASSED, PHIX_MAPPING, PLASMID_MAPPED_BYPOOL, BARCODE_MUX, LTR_IDENTIFIED, TRIMMING_FINAL_LTRLC, LV_MAPPED, BWA_MAPPED_OVERALL, ISS_MAPPED_OVERALL, ISS_MAPPED_PP FROM ${STATSDBNAME} WHERE POOL LIKE '${POOL}'" ${DBSCHEMA} > ${OUTDIR_POOL_ISS}/stats.${DBSCHEMA}.${POOL}.tsv

# merge all together
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Merge BAM with PICARD"
INPUT="INPUT";
INSTR="";
for k in $( ls ${TMPDIR}/bam/*.sorted.md.rel.pg.iss.bam ); do
	INSTR+=" ${INPUT}=${k}";
done
picard-tools MergeSamFiles ${INSTR} OUTPUT="${TMPDIR}/bam/${PATIENT}_${POOL}.wholepool.rel.pg.iss.bam" SORT_ORDER=coordinate ;
# index grouped bam
samtools index ${TMPDIR}/bam/${PATIENT}_${POOL}.wholepool.rel.pg.iss.bam
bamCoverage -p ${MAXTHREADS} -b ${TMPDIR}/bam/${PATIENT}_${POOL}.wholepool.rel.pg.iss.bam --outFileFormat bigwig -o ${TMPDIR}/bam/${PATIENT}_${POOL}.wholepool.rel.pg.iss.bw
bamCoverage -p ${MAXTHREADS} -b ${TMPDIR}/bam/${PATIENT}_${POOL}.wholepool.rel.pg.iss.bam --outFileFormat bedgraph -o ${TMPDIR}/bam/${PATIENT}_${POOL}.wholepool.rel.pg.iss.bedgraph


# # create WHOLE POOL BED with ONLY starting read (r1)
# echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Whole pool analysis: Extract R1 only reads -> IS markers"
# # bamToBed -i ${TMPDIR}/bam/${PATIENT}_${POOL}.wholepool.rel.bam | grep '/1' > ${TMPDIR}/bed/${PATIENT}_${POOL}.wholepool.rel.bed
# # echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Whole pool analysis: Filter by CIGAR and MD tag"
# # filter_by_cigar_tagmd -bed ${TMPDIR}/bed/${PATIENT}_${POOL}.wholepool.rel.bed --bam ${TMPDIR}/bam/${PATIENT}_${POOL}.wholepool.rel.bam -o ${TMPDIR}/bed/${PATIENT}_${POOL}.wholepool.rel.cigarmd.bed
# for b in ${BCLTR[@]}; do
#   for k in ${BCLC[@]}; do
#     echo "

#     TAG:: ${b} ${k} [for whole BED analysis append]

#     where LTR=<${b}> and LC=<${k}> [marginal signs <> are delimiters]
#     " ;
#     cat ${TMPDIR}/bed/${b}.${k}.sorted.md.rel.pg.bed >> ${TMPDIR}/bed/${PATIENT}_${POOL}.wholepool.rel.iss.bed
#   done;
# done
# echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Sort merged BED file"
# bedtools sort -i ${TMPDIR}/bed/${PATIENT}_${POOL}.wholepool.rel.cigarmd.bed > ${TMPDIR}/bed/${PATIENT}_${POOL}.wholepool.rel.cigarmd.sorted.bed

# # determina quality in output
# echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Statistics for reliable merged BAM (both fastqc and samstat)"
# # fastqc -o ${OUTDIR_POOL_QUAL} --contaminants ${SEQCONTAM} -t ${MAXTHREADS} -f bam ${TMPDIR}/bam/${PATIENT}_${POOL}.wholepool.rel.bam
# # samstat ${TMPDIR}/bam/${PATIENT}_${POOL}.wholepool.rel.bam
# for k in $( ls ${TMPDIR}/bam/LTR*bam ); do
# 	samstat $k;
# done


# sync folder from TMP to results BAM files
#mv ${TMPDIR}/bam/LTR* ${OUTDIR_POOL_BAM}/;
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Syncing folder from TMP to results BAM files..."
for k in $( ls ${TMPDIR}/bam/* ); do
	mv $k ${OUTDIR_POOL_BAM}/ ;
done
mv ${TMPDIR}/bam/${PATIENT}_${POOL}.wholepool* ${OUTDIR_MERGE_BAM}/;


# sync folder from TMP to results BED files
pigz ${TMPDIR}/bed/*
mv ${TMPDIR}/bed/* ${OUTDIR_POOL_BED}/;
#mv ${TMPDIR}/bed/${PATIENT}_${POOL}.wholepool* ${OUTDIR_MERGE_BED}/;
#rm ${OUTDIR_MERGE_BED}/${PATIENT}_${POOL}.wholepool.rel*bed
#mv ${TMPDIR}/bed/${PATIENT}_${POOL}.wholepool.rel.iss.bed ${OUTDIR_MERGE_BED}/${PATIENT}_${POOL}.wholepool.rel.iss.bed

# sync folder from TMP to results SAM files
#mv ${TMPDIR}/sam/* ${OUTDIR_POOL_SAM}/;

# sync folder from TMP to results BCMUXALL (FASTQ) files
#mv ${TMPDIR}/bcmuxall/*.noLTRLC.noVector.fastq.gz ${OUTDIR_POOL_BCMUXALL}/;


# echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Collect insert size stats of reliable merged BAM"
# CollectInsertSizeMetrics HISTOGRAM_FILE=${OUTDIR_MERGE_BAM}/${PATIENT}_${POOL}.wholepool.rel.hist.pdf INPUT=${OUTDIR_MERGE_BAM}/${PATIENT}_${POOL}.wholepool.rel.bam OUTPUT=${OUTDIR_MERGE_BAM}/${PATIENT}_${POOL}.wholepool.rel.summary

if [ ${FASTQ_QF} = "slim" ]; then
	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filtering Random TAGs from R2 with final ISs"
	cat ${OUTDIR_POOL_ISS}/${DBSCHEMA}_${DBTABLE}_refactored.tsv | cut -f1 | tail -n +2 > ${TMPDIR}/tmp_header_r2.txt
	zcat ${OUTDIR_POOL_QUAL}/r2.${POOL}.qf.noPlasmids.noPhiX.TAGs.fa.gz | faextract_pureheader ${TMPDIR}/tmp_header_r2.txt | pigz -f -c > ${OUTDIR_POOL_QUAL}/r2.${POOL}.qf.noPlasmids.noPhiX.TAGs.filtered.fa.gz
	mv ${OUTDIR_POOL_QUAL}/r2.${POOL}.qf.noPlasmids.noPhiX.TAGs.fa.gz ${OUTDIR_POOL_BCMUXALL}
	mv ${OUTDIR_POOL_QUAL}/r2.${POOL}.qf.noPlasmids.noPhiX.TAGs.filtered.fa.gz ${OUTDIR_POOL_QUAL}/r2.${POOL}.qf.noPlasmids.noPhiX.TAGs.fa.gz
fi

#### Compressing refactored iss file
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Compressing refactored iss file..."
pigz -f ${OUTDIR_POOL_ISS}/${DBSCHEMA}_${DBTABLE}_refactored.tsv;
pigz -f ${OUTDIR_POOL_ISS}/${DBSCHEMA}_${DBTABLE}_refactored.read2.tsv;


# sync folder from TMP to results BCMUXALL (FASTA) files
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Syncing folder from TMP to results BCMUXALL (FASTA) files..."
for k in $( ls ${TMPDIR}/bcmuxall/*.noLTRLC.Vector.fastq.gz ); do
	NAME=`basename $k .fastq.gz`;
	zcat $k | fastq_to_fasta -n -Q 33 -z > ${OUTDIR_POOL_BCMUXALL}/$NAME.fa.gz
done

# rename files with complete name in association file
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Renaming FASTA (only genomic) files with complete name in association file..."
ASSOBARCODES=(`cut -f2 ${ASSOCIATIONFILE}`);
COMPLETENAME=(`cut -f7 ${ASSOCIATIONFILE}`);
for k in $(ls ${OUTDIR_POOL_BCMUXALL}/r1.*.Vector.fa.gz); do
	for ((j=0;j<${#ASSOBARCODES[@]};++j)); do
		i=`basename ${k}`;
		if [[ `basename ${k} .noLTRLC.Vector.fa.gz` == r1.no12.${ASSOBARCODES[j]} ]]; then mv $k ${OUTDIR_POOL_BCMUXALL}/r1.${COMPLETENAME[j]}.fa.gz; mv ${OUTDIR_POOL_BCMUXALL}/r2.no12.${i:8} ${OUTDIR_POOL_BCMUXALL}/r2.${COMPLETENAME[j]}.fa.gz; fi;
	done
done
if [ ${REPEATS} = "repeats_yes" ]
	then
	for k in $(ls ${OUTDIR_POOL_ISS_REPEATS}/*.sorted.md.lowMQ.R1.fa.out); do
		for ((j=0;j<${#ASSOBARCODES[@]};++j)); do
			i=`basename ${k}`;
			if [[ `basename ${k} .sorted.md.lowMQ.R1.fa.out` == ${ASSOBARCODES[j]} ]]; then mv $k ${OUTDIR_POOL_ISS_REPEATS}/${COMPLETENAME[j]}.lowMQ.R1.fa.out; fi;
		done
	done
	for k in $(ls ${OUTDIR_POOL_ISS_REPEATS}/*.sorted.md.lowMQ.R1.fa.tbl); do
		for ((j=0;j<${#ASSOBARCODES[@]};++j)); do
			i=`basename ${k}`;
			if [[ `basename ${k} .sorted.md.lowMQ.R1.fa.tbl` == ${ASSOBARCODES[j]} ]]; then mv $k ${OUTDIR_POOL_ISS_REPEATS}/${COMPLETENAME[j]}.lowMQ.R1.fa.tbl; fi;
		done
	done
	for k in $(ls ${OUTDIR_POOL_ISS_REPEATS}/*.sorted.md.lowMQ.R1.bam); do
		for ((j=0;j<${#ASSOBARCODES[@]};++j)); do
			i=`basename ${k}`;
			if [[ `basename ${k} .sorted.md.lowMQ.R1.bam` == ${ASSOBARCODES[j]} ]]; then mv $k ${OUTDIR_POOL_ISS_REPEATS}/${COMPLETENAME[j]}.lowMQ.R1.bam; fi;
		done
	done
fi


#### clean tmp dir
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Cleaning TMP Directory... "
if [ ${REMOVE_TMP_DIR} = "remove_tmp_yes" ]
	then
	rm -fr ${TMPDIR};
fi
