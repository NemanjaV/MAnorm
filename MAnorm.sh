#--------------------------------------------------------------------------------------------------------------------------------------------------#
# NOTE: This script was originally downloaded from the official MAnorm website (http://bcb.dfci.harvard.edu/~gcyuan/MAnorm/MAnorm.htm#Downloads),  #
# but it has been modified by Roger Mulet in order to allow its execution from another directory - normally the program would expect the input     #
# files and the other scripts to be in the same directory																						   #
# UPDATE: 26/07/2017																															   #
# LAST UPDATE: Add option to filter narrowPeak files  									  											               # 
#--------------------------------------------------------------------------------------------------------------------------------------------------#

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ $# -lt 4 ]
then
	echo "Usage: `basename $0` peak1.bed peak2.bed (read1.bed read2.bed | read1.bam read2.bam) [bp_shift_1 bp_shift_2] [-o OUTPUT -f FILTER]"

	echo "Optional arguments:"
	echo "-o|--output		Name of the output file"
	echo "-f|--filter		Significance threshold for MACS FDR [0.05]"
	exit
fi

#--------------------------------------------------------------------------------------------------
#Check input arguments
#--------------------------------------------------------------------------------------------------

check_file() {

	file=$1; ext=$2

	if [[ -z $file ]]; then

		echo 'Error: Please specify an input file (-h for help)'
		exit -1

	elif [[ ! -f $file ]]; then

		echo "Error: Specified file $file is not a valid file or does not exist"
		exit -1

	fi

	if [[ ! $file =~ (${ext})$ ]]; then

		echo "Error: Input file $file does not have $ext extension as requested"
		exit -1
	fi
}

check_file $1 "bed|narrowPeak"; PEAKS1=$1
check_file $2 "bed|narrowPeak"; PEAKS2=$2
check_file $3 "bed|bam"; READS1=$3
check_file $4 "bed|bam"; READS2=$4

FILTER=1.30103

while [[ $# -gt 0 ]]
do
        case "$1" in
			-o|--output)
			OUTPUT=$2
			shift
			;;
			-f|--filter)
			FILTER=$2
			shift
			;;
			-h|--help)
			display_usage
			exit 0
			;;
        esac
shift
done 


if [[ $PEAKS1 =~ ".narrowPeak" ]]; then

	echo "Info: $PEAKS1 will be converted to BED using filter $FILTER"
	awk -v FILTER=$FILTER '$7>FILTER{print $1,$2,$3,$4,$7}' $PEAKS1 > ${PEAKS1%%.narrowPeak}.bed

fi

if [[ $PEAKS2 =~ ".narrowPeak" ]]; then

	echo "Info: $PEAKS2 will be converted to BED using filter $FILTER"
	awk -v FILTER=$FILTER '$7>FILTER{print $1,$2,$3,$4,$7}' $PEAKS2 > ${PEAKS2%%.narrowPeak}.bed

fi
	

#if [[ $3 == *.bam || $4 == *.bam ]]
#then
#	echo "Warning: BAM files provided instead of BED files. Reads will be automatically extracted"
#	READS1=${3%%.bam}_reads.bed
#	READS2=${4%%.bam}_reads.bed
#	if [[ ! -f $READS1 ]]
#	then 
#		samtools view $3 | awk -F'\t' '{if ($2==0) {print $3,$4,($4+length($10)-1),"+"} else if ($2==16) {print $3,$4,($4+length($10)-1),"-"}}' > ${3%%.bam}_reads.bed
#	fi
#	if [[ ! -f $READS2 ]]
#	then
#		samtools view $4 | awk -F'\t' '{if ($2==0) {print $3,$4,($4+length($10)-1),"+"} else if ($2==16) {print $3,$4,($4+length($10)-1),"-"}}' > ${4%%.bam}_reads.bed
#	fi
#	if [[ $? -ne 0 ]]
#	then
#		echo 'Error: Extract reads to bed format failed'
#		exit -1
#	fi
#else 
#	READS1=$3
#	READS2=$4
#
#fi

# Added feature to count the number of bp from XLS file (Roger)

if [ -z $5 -o -z $6 ]
then
  echo "Warning: You have not specified a number of bp to be shifted. They will be automatically retrieved from XLS file"

  if [ -e ${PEAKS1%%.bed}.xls -a -e ${PEAKS2%%.bed}.xls ]
  then
    BP1=$(grep -oP '(?<=#\sd\s=\s)\d+' ${PEAKS1%%.bed}.xls)
    BP2=$(grep -oP '(?<=#\sd\s=\s)\d+' ${PEAKS2%%.bed}.xls)
	echo "BP to be used: $BP1 $BP2"
  else 
	echo "Error: BP could not be recovered! Please provide either BP or make sure that XLS file exists"
	exit -1
  fi

else
  BP1=$5
  BP2=$6

fi

if [[ -z $OUTPUT ]]; then

	# Define output name. By default, a combination of the two file names (Roger)
	OUTPUT=$(echo ${READS1##*/} | grep -oP '.*(?=_[12][0-9][012][0-9][0123][0-9])')_vs_$(echo ${READS2#*/} | grep -oP '.*(?=_[12][0-9][012][0-9][0123][0-9])')
	
	echo $OUTPUT  # CHANGE, only for our server
	# echo  $OUTPUT $PEAKS1 $PEAKS2 $READS1 $READS2

fi

# echo  $OUTPUT $PEAKS1 $PEAKS2 $READS1 $READS2
echo $OUTPUT

#--------------------------------------------------------------------------------------------------
#Script execution
#--------------------------------------------------------------------------------------------------

echo "StepI: clean input"

sed 's/\s$//g' $PEAKS1 | awk 'BEGIN {OFS="\t"}
     {if ($1~/chr/ && $1 !="chrM" && $1 !~/random/ && $3>$2 && $2>0 && $3>0)
          print $1,$2,$3>"peak1.bed";
      else 
          print $0 > "peak1_dump.bed"}'

sed 's/\s$//g' $PEAKS2 | awk 'BEGIN {OFS="\t"}
     {if ($$1~/chr/ && 1 !="chrM"  && $1 !~/random/ && $3>$2  && $2>0 && $3>0)
          print $1,$2,$3>"peak2.bed";
      else
          print $0 > "peak2_dump.bed"}'

if [[ $READS1 == *bam ]]; then
	echo "The provided file $READS1 is a BAM file. It will be converted to BED"
	# Simply replace 4 (strand) with 6
	bedtools bamtobed -i $READS1 | awk -v var=$BP1 'BEGIN {OFS="\t"}
	     {if ($1~/chr/ && $1 !="chrM" && $6=="+" && $1 !~/random/ && $3>$2 && $2>0 && $3>0)
		  print $1,$2+var,$3+var>"read1.bed";
	      else if ($1~/chr/  && $1 !="chrM" && $6=="-" && $1 !~/random/ && $3>$2  && $2>var && $3>var)
		  print $1,$2-var,$3-var>"read1.bed";
	      else 
		  print $0 > "read1_dump.bed"}'
else

	sed 's/\s$//g' $READS1 | awk -v var=$BP1 'BEGIN {OFS="\t"}
	     {if ($1~/chr/ && $1 !="chrM" && $4=="+" && $1 !~/random/ && $3>$2 && $2>0 && $3>0)
		  print $1,$2+var,$3+var>"read1.bed";
	      else if ($1~/chr/  && $1 !="chrM" && $4=="-" && $1 !~/random/ && $3>$2  && $2>var && $3>var)
		  print $1,$2-var,$3-var>"read1.bed";
	      else 
		  print $0 > "read1_dump.bed"}'
fi

if [[ $READS2 == *bam ]]; then
	echo "The provided file $READS2 a BAM file. It will be converted to BED"
	# Simply replace 4 (strand) with 6
	bedtools bamtobed -i $READS2 | awk -v var=$BP2 'BEGIN {OFS="\t"}
	     {if ($1~/chr/ && $1 !="chrM" && $6=="+" && $1 !~/random/ && $3>$2 && $2>0 && $3>0)
		  print $1,$2+var,$3+var>"read2.bed";
	      else if ($1~/chr/  && $1 !="chrM" && $6=="-" && $1 !~/random/ && $3>$2  && $2>var && $3>var)
		  print $1,$2-var,$3-var>"read2.bed";
	      else 
		  print $0 > "read2_dump.bed"}'
else

	sed 's/\s$//g' $READS2 | awk -v var=$BP2 'BEGIN {OFS="\t"}
	     {if ($1~/chr/ && $1 !="chrM" && $4=="+" && $1 !~/random/ && $3>$2 && $2>0 && $3>0)
		  print $1,$2+var,$3+var>"read2.bed";
	      else if ($1~/chr/  && $1 !="chrM" && $4=="-" && $1 !~/random/ && $3>$2 && $2>var && $3>var)
		  print $1,$2-var,$3-var>"read2.bed";
	      else 
		  print $0 > "read2_dump.bed"}'
fi

#--------------------------------------------------------------------------------------------------
#Classification of peaks and read counting
#----------------------------------------------------------------------------------------------

echo "StepII: classify common or unique peaks"
intersectBed -a peak1.bed -b peak2.bed -u | sort -k1,1 -k2,2n -k3,3n > common_peak1.bed
intersectBed -a peak2.bed -b peak1.bed -u | sort -k1,1 -k2,2n -k3,3n > common_peak2.bed
intersectBed -a peak1.bed -b peak2.bed -v | sort -k1,1 -k2,2n -k3,3n > unique_peak1.bed
intersectBed -a peak2.bed -b peak1.bed -v | sort -k1,1 -k2,2n -k3,3n > unique_peak2.bed

#cat common_peak1.bed common_peak2.bed | mergeBed -i - > common_peak.bed
# Added sorting to suppress warning messages (RM)
cat common_peak1.bed common_peak2.bed | sort -k1,1 -k2,2n > temp_common_peak.bed
mergeBed -i temp_common_peak.bed > common_peak.bed

echo "StepIII: count peak read"
if [ -f MAnorm.bed ];
then
rm MAnorm.bed
fi
# Modified order of coverageBed (RM)
coverageBed -b read1.bed -a unique_peak1.bed | sort -k1,1 -k2,2n -k3,3n  | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"unique_peak1" >> "MAnorm.bed"; print $4 > "unique_peak1_count_read1"}'
coverageBed -b read2.bed -a unique_peak1.bed | sort -k1,1 -k2,2n -k3,3n  | awk '{print $4 > "unique_peak1_count_read2"}'
coverageBed -b read1.bed -a common_peak1.bed | sort -k1,1 -k2,2n -k3,3n  | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"common_peak1" >> "MAnorm.bed";print $4 > "common_peak1_count_read1"}'
coverageBed -b read2.bed -a common_peak1.bed | sort -k1,1 -k2,2n -k3,3n  | awk '{print $4 > "common_peak1_count_read2"}'
coverageBed -b read1.bed -a common_peak2.bed | sort -k1,1 -k2,2n -k3,3n  | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"common_peak2"  >> "MAnorm.bed";print $4 > "common_peak2_count_read1"}'
coverageBed -b read2.bed -a common_peak2.bed | sort -k1,1 -k2,2n -k3,3n  | awk '{print $4 > "common_peak2_count_read2"}'
coverageBed -b read1.bed -a unique_peak2.bed | sort -k1,1 -k2,2n -k3,3n  | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"unique_peak2">> "MAnorm.bed";print $4 > "unique_peak2_count_read1"}'
coverageBed -b read2.bed -a unique_peak2.bed | sort -k1,1 -k2,2n -k3,3n  | awk '{print $4 > "unique_peak2_count_read2"}'


cat common_peak1_count_read1 common_peak2_count_read1 > common_peak_count_read1
cat common_peak1_count_read2 common_peak2_count_read2 > common_peak_count_read2
cat unique_peak1_count_read1 common_peak1_count_read1 common_peak2_count_read1 unique_peak2_count_read1 > peak_count_read1
cat unique_peak1_count_read2 common_peak1_count_read2 common_peak2_count_read2 unique_peak2_count_read2 > peak_count_read2

if [ -f MAnorm_merge.bed ];
then
rm MAnorm_merge.bed
fi

cat  unique_peak1.bed | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"unique_peak1" >> "MAnorm_merge.bed"}'
coverageBed -b read1.bed -a common_peak.bed | sort -k1,1 -k2,2n -k3,3n  | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"merged_common_peak" >> "MAnorm_merge.bed"; print $4 > "merge_common_read1"}'
coverageBed -b read2.bed -a common_peak.bed | sort -k1,1 -k2,2n -k3,3n  | awk '{print $4 > "merge_common_read2"}'
cat  unique_peak2.bed | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"unique_peak2" >> "MAnorm_merge.bed"}'

cat unique_peak1_count_read1 merge_common_read1  unique_peak2_count_read1 > merge_common_peak_count_read1
cat unique_peak1_count_read2 merge_common_read2  unique_peak2_count_read2 > merge_common_peak_count_read2


echo "StepIV: normalize using common peaks"
#R --vanilla MAnorm.r >Rcommand.out
R CMD BATCH $DIR/MAnorm.r Rcommand.out

if [ $? -ne 0 ]
then
	echo 'Error: MAnorm.r failed! Execution halted.'
	exit -1
fi

# Removed unnecessary sample peaks (Roger Mulet)
#awk 'BEGIN{OFS="\t"}{if($4~/1/) print $1,$2,$3,$7>"sample1_peaks.wig"}' MAnorm_result.xls
#awk 'BEGIN{OFS="\t"}{if($4~/2/) print $1,$2,$3,$7>"sample2_peaks.wig"}' MAnorm_result.xls

rm temp_common_peak.bed
rm *count*
rm *read1*
rm *read2*
rm *peak1*
rm *peak2*
rm MAnorm.bed
rm MAnorm_merge.bed
rm common_peak.bed

echo "INFO: Renaming files as specified"
rename "s/MAnorm_result/MAnorm_${OUTPUT}_result/" MAnorm*
rename "s/MAplot_after/MAplot_${OUTPUT}_after/" MAplot*
rename "s/MAplot_before/MAplot_${OUTPUT}_before/" MAplot*

echo -e "MAnorm.sh execution finished. Files saved as ${OUTPUT}\n"
