Collection of scripts to perform the new simulations and analyses found in the 2014 bioRxiv paper: Deleterious mutations can contribute to the evolution of recombination suppression between sex chromosomes


# Usage:

## Step 1: Create the necessary directories:
>mkdir ../InitialState

>mkdir ../Output

## Step 2: Creare the initial population:
Make a simulation of a population with only an autosome or a sex chromosome during 200.000 generations:
>slim -d N=12500 -d s=-0.001 -d Rep=1 ScriptNeutralInversion_DefineInitialStateSexChrom_NoFullRecess.slim #Here with sex chromosome, N=12500, mean of gamma distribution of s=-0.001 et repetition n°1, and a scenario without fully recessive mutation

## Step 3: Introduce inversions in the population
>source Parallele_IntroInv_XY_Sex_NoFullRecess.sh

>export -f slimFc

>parallel -j 30 slimFc ::: 12500 ::: -0.001 ::: 2500000 ::: 2000000 ::: 0 ::: {1..100000} #Here, introduce 100000 independant inversions of 2000000bp, centred on position 2500000 (middle of the chromosome), with the same parameter values as above. These inversions have no ("0") initial selective advantage, except the one provided by the putative lower load


The output file is: IntroInv_g200000_XYsyst_5M_Sex_N12500_r1e-8_u1e-08_s-0.001_sInv0_hmean0.2_Stat.txt
The column names are:
"N", "mu", "r", "h", "s", "sInv", "sim_cycle" (generation), "Inversion_start", "Inversion_end", "Repetition", "Chromosome", "mean number of mutations in the inverted segments", "minimum number of mutations in the inverted segments", "maximum number of mutations in the inverted segments", "standard deviation of the number of mutations in the inverted segments", "mean frequency of mutations in the inverted segments", "mean number of mutations in the non-inverted segments", "minimum number of mutations in the non-inverted segments", "maximum number of mutations in the non-inverted segments", "standard deviation of the number of mutations in the non-inverted segments", "mean frequency of mutations in the non-inverted segments", "Mean Fitness of inversion carriers", "mean Fitness of non-inverted segment carriers", "Frequency of the inversion"

A summary of the number of fixed inversions can be obtained using the script:

>python Summarise_NoChunk.py -i IntroInv_g200000_XYsyst_5M_Sex_N12500_r1e-8_u1e-08_s-0.001_sInv0_hmean0.2_Stat.txt -o IntroInv_g200000_XYsyst_5M_Sex_N12500_r1e-8_u1e-08_s-0.001_sInv0_hmean0.2_Stat.summary.txt

These summary files can be used to plot the data using the Plot_SexChromSheltering.R script

Please contact me paul.yann.jay[at]gmail.com for any further information or help running these scripts
