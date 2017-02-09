DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ $# -ne 6 ]
then
  echo "Usage: `basename $0` peak1.bed peak2.bed read1.bed read2.bed bp_shift_1 bp_shift_2"
  exit
fi

if [ -z $5 -o -z $6 ]
then
  echo "Warning: You have not specified a number of bp to be shifted. They will be automatically retrieved from XLS file"

  if [ -e $(basename $1 bed).xls -a -e $(basename $2 bed).xls ]
  then
    BP1=grep -oP '(?<=#\sd\s=\s)\d+' $(basename $1 bed).xls
    BP2=grep -oP '(?<=#\sd\s=\s)\d+' $(basename $2 bed).xls
  fi
else
  BP1=$5
  BP2=$6
fi

echo "StepI: clean input"
sed 's/\s$//g' $1 | awk 'BEGIN {OFS="\t"}
     {if ($1~/chr/ && $1 !="chrM" && $1 !~/random/ && $3>$2 && $2>0 && $3>0)
          print $1,$2,$3>"peak1.bed";
      else
          print $0 > "peak1_dump.bed"}'
sed 's/\s$//g' $2 | awk 'BEGIN {OFS="\t"}
     {if ($$1~/chr/ && 1 !="chrM"  && $1 !~/random/ && $3>$2  && $2>0 && $3>0)
          print $1,$2,$3>"peak2.bed";
      else
          print $0 > "peak2_dump.bed"}'
sed 's/\s$//g' $3 | awk -v var=$BP1 'BEGIN {OFS="\t"}
     {if ($1~/chr/ && $1 !="chrM" && $4=="+" && $1 !~/random/ && $3>$2  && $2>0 && $3>0)
          print $1,$2+var,$3+var>"read1.bed";
      else if ($1~/chr/  && $1 !="chrM" && $4=="-" && $1 !~/random/ && $3>$2  && $2>var && $3>var)
          print $1,$2-var,$3-var>"read1.bed";
      else
          print $0 > "read1_dump.bed"}'
sed 's/\s$//g' $4 | awk -v var=$BP2 'BEGIN {OFS="\t"}
     {if ($1~/chr/ && $1 !="chrM" && $4=="+" && $1 !~/random/ && $3>$2  && $2>0 && $3>0)
          print $1,$2+var,$3+var>"read2.bed";
      else if ($1~/chr/  && $1 !="chrM" && $4=="-" && $1 !~/random/ && $3>$2  && $2>var && $3>var)
          print $1,$2-var,$3-var>"read2.bed";
      else
          print $0 > "read2_dump.bed"}'


echo "StepII: classify common or unique peaks"
intersectBed -a peak1.bed -b peak2.bed -u | sort -k1,1 -k2,2n -k3,3n > common_peak1.bed
intersectBed -a peak2.bed -b peak1.bed -u | sort -k1,1 -k2,2n -k3,3n > common_peak2.bed
intersectBed -a peak1.bed -b peak2.bed -v | sort -k1,1 -k2,2n -k3,3n > unique_peak1.bed
intersectBed -a peak2.bed -b peak1.bed -v | sort -k1,1 -k2,2n -k3,3n > unique_peak2.bed

#cat common_peak1.bed common_peak2.bed | mergeBed -i - > common_peak.bed
cat common_peak1.bed common_peak2.bed > temp_common_peak.bed
mergeBed -i temp_common_peak.bed > common_peak.bed



echo "StepIII: count peak read"
if [ -f MAnorm.bed ];
then
rm MAnorm.bed
fi
coverageBed -b read1.bed -a unique_peak1.bed | sort -k1,1 -k2,2n -k3,3n  | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"unique_peak1" >> "MAnorm.bed"; print $4 > "unique_peak1_count_read1"}'
coverageBed -b read2.bed -a unique_peak1.bed | sort -k1,1 -k2,2n -k3,3n  | awk '{print $4 > "unique_peak1_count_read2"}'
coverageBed -b read1.bed -a common_peak1.bed | sort -k1,1 -k2,2n -k3,3n  | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"common_peak1" >> "MAnorm.bed";print $4 > "common_peak1_count_read1"}'
coverageBed -b read2.bed -a common_peak1.bed | sort -k1,1 -k2,2n -k3,3n  | awk '{print $4 > "common_peak1_count_read2"}'
coverageBed -b read1.bed -a common_peak2.bed | sort -k1,1 -k2,2n -k3,3n  | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"common_peak2"  >> "MAnorm.bed";print $4 > "common_peak2_count_read1"}'
coverageBed -b read2.bed -a common_peak2.bed |sort -k1,1 -k2,2n -k3,3n  |  awk '{print $4 > "common_peak2_count_read2"}'
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




echo "SetpIV: normalize using common peaks"
#R --vanilla MAnorm.r >Rcommand.out
R CMD BATCH $DIR/MAnorm.r Rcommand.out

if [ $? -ne 0 ]
then
	echo 'Error: MAnorm.r failed! Execution halted.'
	exit -1
fi

awk 'BEGIN{OFS="\t"}{if($4~/1/) print $1,$2,$3,$7>"sample1_peaks.wig"}' MAnorm_result.xls
awk 'BEGIN{OFS="\t"}{if($4~/2/) print $1,$2,$3,$7>"sample2_peaks.wig"}' MAnorm_result.xls


rm temp_common_peak.bed
rm *count*
rm *read1*
rm *read2*
rm *peak1*
rm *peak2*
rm MAnorm.bed
rm MAnorm_merge.bed
rm common_peak.bed
