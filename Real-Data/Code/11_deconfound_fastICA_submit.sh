#!/bin/bash

cd /home/zwan873/Real-Data/Data/deconfound_fastICA/bash_outputs

for seedID in $(seq 1 400)
do

sed -e s:seedID:"${seedID}":g </home/zwan873/Real-Data/Code/11_deconfound_singleseed_fastICA.R >/home/zwan873/Real-Data/Data/deconfound_fastICA/bash_outputs/deconfound_singleseed_fastICA_${seedID}.R

sed -e s:seedID:"${seedID}":g </home/zwan873/Real-Data/Code/11_deconfound_singleseed_fastICA.sh >/home/zwan873/Real-Data/Data/deconfound_fastICA/bash_outputs/deconfound_singleseed_fastICA_${seedID}.sh

qsub  -cwd -N seed${seedID} deconfound_singleseed_fastICA_${seedID}.sh
#qsub  -cwd -q R4.q -N seed${seedID} deconfound_singleseed_fastICA_${seedID}.sh

done

