#!/bin/bash



Help()
{
   # Display Help
   echo "Generating a .csv file of all sequences included in each MOTU from the output of SWARM"
   echo
   echo "Syntax: bash MOTUs_from_SWARM.sh [-h] [-i] [motu_list_file] [-o] [output_file_from_SWARM] [-r] [-t] [sample_abundances] [-d] [output_directory] [-l] [output_file_from_lulu]"
   echo "options:"
   echo "h     Print this Help."
   echo "i     .txt containing MOTU ids for which to create .csv files"
   echo "      example: cat motus_to_denoise.txt"
   echo "      		seq1_id"
   echo "      		seq2_id"
   echo "      		seq3_id"
   echo "o     swarm output file"
   echo "r     remove databases that will be created during the process in the output directory"
   echo "t     .tab file containing sample information of original sequences (from obitab)"
   echo "d     output directory"
   echo "l     lulu corrected_sequences file, an output file from lulu (optional) "
   echo
}

while getopts hri:o:t:d:l: flag
do
    case "${flag}" in
	h) Help
		exit;;
	i) motu_list="$( cd -P "$( dirname "${OPTARG}" )" >/dev/null 2>&1 && pwd )/${OPTARG}";;
	o) output_SWARM="$( cd -P "$( dirname "${OPTARG}" )" >/dev/null 2>&1 && pwd )/${OPTARG}";;
	r) remove_db=true;;
	t) tab_file="$( cd -P "$( dirname "${OPTARG}" )" >/dev/null 2>&1 && pwd )/${OPTARG}";;
	d) output_dir="$( cd -P "$( dirname "${OPTARG}" )" >/dev/null 2>&1 && pwd )/";;
	l) lulu_deleted=${OPTARG};;
	\?) echo "usage: bash MOTUs_from_SWARM.sh [-h|i|o|r|t|d|l]"
		exit;;
    esac
done

if [ -z "${motu_list}" ]
 then
 echo 'ERROR! motu_list (-i) needed'
 Help
 exit
 fi
if [ -z "${output_SWARM}" ]
 then
 echo 'ERROR! SWARM_output_file (-o) needed'
 Help
 exit
 fi
if [ -z "${tab_file}" ]
 then
 echo 'ERROR! tab_file (-t) needed'
 Help
 exit
 fi
if [ -z "${output_dir}" ]
 then
 echo 'ERROR! output_dir (-d) needed'
 Help
 exit
 fi
if [ -z "${remove_db}" ]
 then
 remove_db=false
else
 echo database will be removed at the end
fi


SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  TARGET="$(readlink "$SOURCE")"
  if [[ $TARGET == /* ]]; then
    SOURCE="$TARGET"
  else
    DIR="$( dirname "$SOURCE" )"
    SOURCE="$DIR/$TARGET" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
  fi
done

script_dir="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )/"




mkdir ${output_dir}database
db_dir=${output_dir}database/

mkdir ${db_dir}motus
motus_dir=${db_dir}motus/

mkdir ${output_dir}motus_csv
motus_csv_dir=${output_dir}motus_csv/

 
linies=$(wc -l ${motu_list} | cut -f1 -d ' ')

# creating a version of output_swarm with just the wanted MOTUs

motus_pipeline=$(sed -e ':a;N;$!ba;s/\n/\\|/g' ${motu_list}) #converts line feeds to "\|"
grep ${motus_pipeline} ${output_SWARM} >${db_dir}output_swarm

# for i in $(seq 1 ${linies}) 
#  do 
#  var=$(sed "${i}q;d" ${motu_list})
#  if [ i = 1 ]
#   then
#   grep ${var} ${output_SWARM} >${db_dir}output_swarm
#  else
#   grep ${var} ${output_SWARM} >>${db_dir}output_swarm
#  fi
# done
 
output_SWARM=${db_dir}output_swarm


sed -i -e "s/;size=[0-9]*;/;/g" ${output_SWARM}

if [ -n "${lulu_deleted}" ]
 then
 for i in $(seq 1 ${linies}) 
  do 
  var=$(sed "${i}q;d" ${motu_list})
  grep ${var} ${output_SWARM} >${motus_dir}${var}_nolulu
  sed -i -e "s/; /\\n/g" ${motus_dir}${var}_nolulu
  sed -i -e "s/;//g" ${motus_dir}${var}_nolulu
  Rscript ${scripts_dir}add_lulu2MOTUfiles.R ${var} ${lulu_deleted} ${motus_dir}${var}_lulu
  cat ${motus_dir}${var}_nolulu ${motus_dir}${var}_lulu > ${motus_dir}${var} & 
  done
 else
  for i in $(seq 1 ${linies}) 
  do 
  var=$(sed "${i}q;d" ${motu_list})
  grep ${var} ${output_SWARM} >${motus_dir}${var}
  sed -i -e "s/; /\\n/g" ${motus_dir}${var}
  sed -i -e "s/;//g" ${motus_dir}${var} &
  done
fi # [ -n "${lulu_deleted}" ]
wait
echo motu composition done



for i in $(seq 1 ${linies}) 
 do
 var=$(sed "${i}q;d" ${motu_list})
 sed "1q;d" ${tab_file} >${motus_csv_dir}${var}
 mkdir ${motus_dir}${var}_dir
 rasputin=${motus_dir}${var}_dir/
 cd ${rasputin}
 split -l200 -d ${motus_dir}${var}
 for z in *
  do
  liniesmotu=$(sed -e ':a;N;$!ba;s/\n/\\|/g' ${rasputin}${z}) #converts line feeds to "\|"
  grep ${liniesmotu} ${tab_file} >${rasputin}${z}_grep&
  done
 wait
 cat ${rasputin}*_grep >>${motus_csv_dir}${var}
 echo MOTU ${i} finished | rm -r ${rasputin} & 
 done
echo now just wait please
wait
echo motu_csv finished


if ${remove_db}
 then
 echo removing database
 rm -r ${db_dir}
fi



