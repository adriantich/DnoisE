#!/bin/bash

set -e

# Test both native python and nuitka binary with small input data set

for DNOISE in dnoise DnoisE.bin; do
  echo $DNOISE
  ## csv input file

  # 1. with entropy correction
  $DNOISE --csv_input sample1000.csv --csv_output output -n count -s 4 -z 79 -y
  # 2. without entropy correction
  $DNOISE --csv_input sample1000.csv --csv_output output -n count -s 4 -z 79
  # 3. parallel
  $DNOISE --csv_input sample1000.csv --csv_output output -n count -s 4 -z 79 -c 2
  # 4. parallel within MOTU
  $DNOISE --csv_input sample1000.csv --csv_output output -n count -s 4 -z 79 -c 2 -w MOTU -p 2


  ## fasta input file

  # 1. with entropy correction
  $DNOISE --fasta_input sample1000.fasta --fasta_output output -y
  # 2. without entropy correction
  $DNOISE --fasta_input sample1000.fasta --fasta_output output
  # 3. parallel
  $DNOISE --fasta_input sample1000.fasta --fasta_output output -c 2



  ## fastq input file

  # 1. with entropy correction
  $DNOISE --fastq_input sample1000.fastq --fasta_output output -y
  # 2. without entropy correction
  $DNOISE --fastq_input sample1000.fastq --fasta_output output
  # 3. parallel
  $DNOISE --fastq_input sample1000.fastq --fasta_output output -c 2


  # ## merging from info
  # dnoise --fasta_input ./test-DnoisE/DnoisE_example.fasta --fasta_output output --joining_file DnoisE_output_example_info.csv -j 4
  # # using modal_length and unique_length
  # dnoise --fasta_input ./test-DnoisE/DnoisE_example.fasta --fasta_output output -y -c 3 -m 313 -u
done
