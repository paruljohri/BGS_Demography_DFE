#!/bin/env Rscript
#This is the final list of stats:
#chmod +x get_final_statistics.R
#Rscript ./get_final_statistics_dpgp3.R 5p derived numbp50 0.8
#setwd("/home/pjohri/eqm_disc_5/")

#options("winSize" = ()$win_size)
options(scipen=999)
args = commandArgs(trailingOnly=TRUE)
side <- args[1]
state <- args[2]
cutoffLink <- args[3]
cons_cutoff <- args[4]

t <- read.table(paste("/Users/paruljohri/BgsDfeDemo/Drosophila/dpgp3_statistics/InsectCons_", cons_cutoff, "/SingExon_", side, "_", state, "_bigwindow.stats", sep=""),h=T,fill=T)


#To calculate mean and variance of all statistics per window:
t_windows <- c("functional", "linked", "neutral")
for (win in t_windows) {
        wpii <- t$thetapi[which(t$WinType==win)]/t$WinSize[which(t$WinType==win)]
        wwi <- t$thetaw[which(t$WinType==win)]/t$WinSize[which(t$WinType==win)]
        whi <- t$thetah[which(t$WinType==win)]/t$WinSize[which(t$WinType==win)]
        whpi <- t$hprime[which(t$WinType==win)]
        #wtdi <- t$tajimasd[which(t$WinType==win)]
        wtdi <- t$tajimasd[which(t$WinType==win & is.finite(t$tajimasd))]
        wsingi <- t$numSing[which(t$WinType==win)]/t$WinSize[which(t$WinType==win)]
        whapdivi <- t$hapdiv[which(t$WinType==win)]
        wrsqi <- as.numeric(as.character(t$rsq[which(t$WinType==win & t$rsq!="<NA>")]))
        wDi <- t$D[which(t$WinType==win)]
        wDpri <- t$Dprime[which(t$WinType==win)]
        wdivi <- t$div[which(t$WinType==win)]/(t$WinSize[which(t$WinType==win)] - t$S[which(t$WinType==win)])
        v_win_m <- c(mean(wpii, na.rm=T), mean(wwi, na.rm=T), mean(whi, na.rm=T), mean(whpi, na.rm=T), mean(wtdi, na.rm=T), mean(wsingi, na.rm=T), mean(whapdivi, na.rm=T), mean(wrsqi, na.rm=T), mean(wDi, na.rm=T), mean(wDpri, na.rm=T), mean(wdivi, na.rm=T))
        v_win_sd <- c(sd(wpii, na.rm=T), sd(wwi, na.rm=T), sd(whi, na.rm=T), sd(whpi, na.rm=T), sd(wtdi, na.rm=T), sd(wsingi, na.rm=T), sd(whapdivi, na.rm=T), sd(wrsqi, na.rm=T), sd(wDi, na.rm=T), sd(wDpri, na.rm=T), sd(wdivi, na.rm=T))
        if (win=="functional"){
                stats <- matrix(c(v_win_m, v_win_sd))
                colnames(stats) <- "win_name"
                rownames(stats) <- c("thetapi_m", "thetaw_m", "thetah_m", "hprime_m", "tajimasd_m", "numSing_m", "hapdiv_m", "rsq_m", "D_m", "Dprime_m", "div_m", "thetapi_sd", "thetaw_sd", "thetah_sd", "hprime_sd", "tajimasd_sd", "numSing_sd", "hapdiv_sd", "rsq_sd", "D_sd", "Dprime_sd", "div_sd")
                }
        else{

                stats <- cbind(stats, c(v_win_m, v_win_sd))
                }
        }
colnames(stats, do.NULL = FALSE)
colnames(stats) <- t_windows

write.table(stats, file=paste("/Users/paruljohri/BgsDfeDemo/Drosophila/dpgp3_statistics/InsectCons_", cons_cutoff, "/SingExon_", side, "_", state, ".bigwinsummary", sep=""), sep="\t")

print ("done")


