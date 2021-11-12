library(PHONEMeS)
library(readr)

CPTAC_phospho_ttop <- as.data.frame(
  read_csv("data/CPTAC_phospho_ttop.csv"))

CPTAC_Kinase_Activities_omnipath <- as.data.frame(
  read_csv("data/CPTAC_Kinase_Activities_omnipath.csv"))

# data("phonemesPKN")

# phonemesPKN <- PHONEMeSPKN

CPTAC_phospho_ttop <- CPTAC_phospho_ttop[
  CPTAC_phospho_ttop$ID %in% phonemesPKN$target,]
CPTAC_Kinase_Activities_omnipath <- CPTAC_Kinase_Activities_omnipath[
  CPTAC_Kinase_Activities_omnipath$X1 %in% phonemesPKN$source,]

kinase_to_exclude <- CPTAC_Kinase_Activities_omnipath[
  abs(CPTAC_Kinase_Activities_omnipath$V1) < 0.5,"X1"]

CPTAC_Kinase_Activities_omnipath <- CPTAC_Kinase_Activities_omnipath[
  abs(CPTAC_Kinase_Activities_omnipath$V1) >= 2,] #2

CPTAC_phospho_ttop <- CPTAC_phospho_ttop[1:300, c(1,4)] #300

kinase_input <- as.numeric(CPTAC_Kinase_Activities_omnipath[,2])
phospho_meas <- as.numeric(CPTAC_phospho_ttop[,2])

names(kinase_input) <- CPTAC_Kinase_Activities_omnipath$X1
names(phospho_meas) <- CPTAC_phospho_ttop$ID
kinase_input <- sign(kinase_input)

PHONEMeS_res <- run_phonemes(inputObj = kinase_input, 
                               measObj = phospho_meas, 
                               rmNodes = kinase_to_exclude, 
                               pruning = T, 
                               n_steps_pruning = 10, 
                               solverPath = "~/Documents/cplex", ##Put whatever you need !!
                               timelimit = 120)

View(PHONEMeS_res$res$nodesAttributes)

PHONEMeS_res <- reattach_psites(PHONEMeS_res)

sif <- as.data.frame(PHONEMeS_res$res$weightedSIF)
att <- as.data.frame(PHONEMeS_res$res$nodesAttributes)

write_csv(sif, file = "sif.csv")
write_csv(att, file = "att.csv")

save.image("session_result.Rdata")