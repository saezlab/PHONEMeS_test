library(readr)
library(PHONEMeS)

load("session_result.Rdata")

pred_reg_sites <- regulatory_psites(PHONEMeS_res)
pred_reg_sites <- pred_reg_sites$regulatory_psites

benchmark <- benchmark_regulatory_psites(regulatory_psites = pred_reg_sites)
