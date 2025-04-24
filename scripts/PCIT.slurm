#!/bin/bash
#PBS -N PCIT
#PBS -q medium
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=2
#PBS -l mem=16gb
#PBS -o PCIT_Pipeline.out
#PBS -e PCIT_Pipeline.err
#PBS -V


cd /home/aubclsd0330/PCIT

echo input_PCIT_ok.txt | ./PCIT.exe

awk '$4 != 0.00000 {print $1, $2, $3, ($3 < 0 ? "NEG" : "POS")}' input_PCIT_ok.txt_out.txt > PCIT_final_signif.txt
