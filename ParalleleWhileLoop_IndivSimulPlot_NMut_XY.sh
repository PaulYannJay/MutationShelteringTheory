#!/bin/bash

slimFc(){
#Init=`ls ../InitialState/slim_g15000_N=1000_r=1e-05_u=$1_s=$3_h=$2_*`
Init=`ls ../InitialState/slim_g15000_MidSDR_10MChrom_XY_N=$1_r=1e-06_u=$2_s=$4_h=$3_*[0-9].txt`
position=$6
size=$7
end=$((position+size/2+1))
start=$((position-size/2+1))
echo $end
echo $start
~/Software/SLiM/build/slim -d N=$1 -d mu=$2 -d h=$3 -d s=$4 -d r=1e-6 -d rep=$5 -d start=$start -d end=$end -d "Init='$Init'" -d "SexChrom='$8'" Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_NMut_XY.slim
}



