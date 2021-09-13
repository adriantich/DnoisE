#!/bin/bash



Help()
{
   # Display Help
   echo "Converting multiple line sequence fasta to unique line"
   echo
   echo "Syntax: bash to_uniq_line_fasta.sh [-h] [-i] [input_file] [-o] [output_file] "
   echo "options:"
   echo "h     Print this Help."
   echo "i     .fasta containing all sequences"
   echo "o     output file"
   echo
}

while getopts hi:o: flag
do
    case "${flag}" in
	h) Help
		exit;;
	i) input_file="$( cd -P "$( dirname "${OPTARG}" )" >/dev/null 2>&1 && pwd )/${OPTARG}";;
	o) output_file="$( cd -P "$( dirname "${OPTARG}" )" >/dev/null 2>&1 && pwd )/${OPTARG}";;
	\?) echo "usage: bash to_uniq_line_fasta.sh [-h|i|o|r|t|d|l]"
		exit;;
    esac
done

if [ -z "${input_file}" ]
 then
 echo 'ERROR! input_file (-i) needed'
 Help
 exit
 fi
if [ -z "${output_file}" ]
 then
 echo 'ERROR! output_file (-o) needed'
 Help
 exit
 fi

sed -E ':a;N;$!ba;s/([ACGTacgt])\n([ACGTacgt])/\1\2/g' ${input_file} > ${output_file}