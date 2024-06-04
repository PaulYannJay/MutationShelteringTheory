library(tidyverse)
library(cowplot)
## Load the file containing the simulation results (The same as in Jay et al. (2022), to be downloaded here: dx.doi.org/10.6084/m9.figshare.19961033)
Simul=read.table(paste("~/Paper/ModelSexChrom/V3/CleanDataset/Compressed/InversionTrajectories_N=1000_Fig3-S13-S15-S19.txt",sep=""), stringsAsFactors = F) #File containing all simulation with N=1000

#Set the column names
colnames(Simul)=c("N", "u", "r", "h", "s", "Gen", "StartInv", "EndInv", "Rep", 
                  "MeanMutInv","MinMutInv","MaxMutInv","sdMutInv","FreqMutInv",
                  "MeanMutNoInv","MinMutNoInv","MaxMutNoInv","sdMutNoInv","FreqMutNoInv",
                  "InvFit", "NoInvFit","Freq","Chromosome")

# Inversions between position 1-10Mb are on sex chromosome, Inversions between position 10-20Mb are on Autosome. Modify the "Chromosome" column to indicate this
Simul[Simul$StartInv>10000000,]$Chromosome="Autosome"

# Create a new column indicating the colums size
Simul$InvSize=Simul$EndInv - Simul$StartInv 

# For each simulation, determine the maximum frequency the inversion reached, and the maximum generation it segregated (I.e. the generation were it becomes lost of fixed)
SimulSum=Simul %>% group_by(u,r,h,s,StartInv, InvSize, Chromosome, Rep) %>% summarise(maxGen=max(Gen), Freq=max(Freq))

#Save the resulting table, to avoid regenerating it
write.table(SimulSum, "~/Projects/MutationShelteringV2/Output/OldDataJayEtAl2022/Reanalyses_Jay_etAl2022/InversionTrajectories_N=1000_Fig3-S13-S15-S19_Summary.txt", quote=F, row.names = F)

# Load this Table
SimulSum=read.table("~/Projects/MutationShelteringV2/Output/OldDataJayEtAl2022/Reanalyses_Jay_etAl2022/InversionTrajectories_N=1000_Fig3-S13-S15-S19_Summary.txt", header=T)

SimulSum$Fixed=0 #Create a new column that indicate if the inversion became fixed. 
SimulSum[SimulSum$Freq>0.95,]$Fixed=1 # If maximum frequency reached by the inversion is above 0.95, it is considered fixed (For computation purpose, we stopped the simulation when the inversion get fixed, and so we do not observe inversion at frequency 1.0)
SimulSum[(SimulSum$Chromosome=="Y" & SimulSum$Freq==0.25),]$Fixed=1 # If the inversion is on chromosome Y and reach a frequency of 0.25, it is fixed
SimulSumSum=SimulSum %>% group_by(u,r,h,s,StartInv, InvSize, Chromosome) %>% summarise(nFixed=sum(Fixed), n=n()) #For each value of u, r, h, s, Inversion size, and inversion position (Y or autosome), count the number of inversion fixed, and of inversion analysed
SimulSumSum$FracFixed=SimulSumSum$nFixed/SimulSumSum$n #Define the proportion of inversion fixed. 
Col=scales::viridis_pal(begin=0.1, end=0.8, option="A")(5) #Define the color palette to be used

SimulSumSum$InvSize_f=factor(paste0("Inv. size=", SimulSumSum$InvSize, " bp"),
                             levels=c("Inv. size=500000 bp", "Inv. size=1000000 bp","Inv. size=2000000 bp","Inv. size=5000000 bp")) #Refactor the the Inv. size, so it appear in the right order)

for (u in unique(SimulSumSum$u)) #For each value of u (1e-08 or 1e-09), create a plot
{
    Base=ggplot(SimulSumSum[(SimulSumSum$u==u & SimulSumSum$s<0),], aes( y=FracFixed))  # Define the plot basics # If you want to plot the simulation corresponding to s=0, just remove the "& SimulSumSum$s<0", and set above : "Col=scales::viridis_pal(begin=0.1, end=0.8, option="A")(6)" (But again, it is not informative, it does not allow to measure the strength of the sheltering effect or of the lower-load advantage, nor allow to measure the effect of different effective population size between autosomes and sex chromosomes on the fixation of loaded inversions) 
    PlotFracFixed=Base+geom_line(aes(x=h, color=as.factor(s)), size=1)+ #Plot the line
        geom_point(aes(x=h, color=as.factor(s)), size=5)+
        scale_color_manual("Selection coefficient \nof mutations (s)", values=Col)+
        ylab("Fraction of inversions fixed after 10,000 generations")+
        xlab("Dominance coefficient of mutations (h)")+
        theme( #Improve plot design
            panel.grid.major = element_line(color="grey95"),
            panel.background = element_blank(),
            text = element_text(size=18),
            axis.line = element_line(colour = "grey"),
            axis.title = element_text(size = 20, face="bold"),
            legend.spacing.y= unit(0, 'cm'),
            legend.title = element_text(size = 20, face="bold"),
            legend.text = element_text(size = 18),
            strip.text.x = element_text(size = 20, face="bold"),
            strip.text.y = element_text(size = 18, face="bold")
        )+
        facet_grid(InvSize_f ~ Chromosome, scales = "free_y")
    
    save_plot(paste0("~/Projects/MutationShelteringV2/Output/OldDataJayEtAl2022/Reanalyses_Jay_etAl2022/Figure_i_u",u,".pdf"), PlotFracFixed, nrow=4, ncol=2, base_aspect_ratio = 2)
    save_plot(paste0("~/Projects/MutationShelteringV2/Output/OldDataJayEtAl2022/Reanalyses_Jay_etAl2022/Figure_i_u",u,".png"), PlotFracFixed, nrow=4, ncol=2, base_aspect_ratio = 2)
}

SimulSumSumWide_FracFixed=SimulSumSum %>%
    pivot_wider(id_cols=-c(nFixed, n, StartInv), values_from = FracFixed, names_from = Chromosome) %>%
    mutate(RatioYvsAuto=Y/Autosome) #Reshape the dataframe to compare the rate of fixation on sex chromosomes and autosomes. The column "Y" contains de rate of inversion fixation on Y chromosomes, while "Autosome" contains de rate of inversion fixation on autosomes. The Column "RatioYvsAuto" contains the ratio of these value (>1 if the rate is higher on sex chromosomes, <1 otherwise)
sum(SimulSumSumWide_FracFixed$RatioYvsAuto>1)/nrow(SimulSumSumWide_FracFixed) #Inversion are more likely to fix on Y chromosomes than on Autosomes in 94.9 % of the parameter space.
sum(SimulSumSumWide_FracFixed$Y)/sum(SimulSumSumWide_FracFixed$Autosome) #Inversion were 5.81 more likely to fix on Y chromosomes than on Autosomes in the parameter space studied.
