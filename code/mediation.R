library(vroom)
library(data.table)
# library(openxlsx)
library(vroom)
library(lavaan)
library(dplyr)


# exposure<-read.csv("/home/cmet-standard/data/OneDrive/T2DM-IDPS-CF-research/Summary_final_v1/Data/Prepared_data/HbA1c.csv")

exposure<-read.csv("/home/cmet-standard/data/OneDrive/T2DM-IDPS-CF-research/Summary_final_v1/Data/Prepared_data/Diabetes_diagnosis.csv")
exposure <- exposure[!grepl("T1D", exposure$diagnosis_result), ]
# table(exposure$diagnosis_result)
exposure$diagnosis_result <- replace(exposure$diagnosis_result, exposure$diagnosis_result == "T2D", 1)
colnames(exposure)[5]<-"diagnosis_result_T2D"
Covariance<-vroom("/home/cmet-standard/data/OneDrive/T2DM-IDPS-CF-research/Summary_final_v1/Data/Prepared_data/Covariance.csv")
outcome<-vroom("/home/cmet-standard/data/OneDrive/T2DM-IDPS-CF-research/Summary_final_v1/Data/Prepared_data/Sleep_diagnosis.csv")
IDP<-vroom("/home/cmet-standard/data/OneDrive/T2DM-IDPS-CF-research/Summary_final_v1/Data/Prepared_data/Imaging_addending.csv")


DF<-merge(IDP, outcome, by="eid")
DF<-merge(DF, Covariance, by="eid")
DF<-merge(DF, exposure, by="eid")

table(DF$diagnosis_result)

med_func <- function(DATA, M_IDX){ 
  # M_IDX<-7
  # DATA = T2D_AD_data
  # DATA$X <- T2D_AD_data$Diabetes
  # DATA$Y <- T2D_AD_data$AD_PRS
  
  DATA$M <- scale(DATA[M_IDX], center = TRUE, scale = TRUE)
  
  mod <- "# a path
     M ~ a * X

     # b path
     Y ~ b * M + age + sex + ethic + BMI + alcohol + smoke + Townsend + edu

     # c prime path 
     Y ~ cp * X

     # indirect and total effects
     ab := a * b
     total := cp + ab"
  
  set.seed(1234)
  fsem <- sem(mod, data = DATA, se = "bootstrap", bootstrap = 1000)
  # pval<-summary(fsem, standardized = TRUE)
  # temp_res <- c(as.numeric(colnames(DATA)[M_IDX]), pval$pe$est[c(1, 2, 5, 14, 15)], pval$pe$pvalue[c(1, 2, 5, 14, 15)])
  # print(temp_res)
  # return(temp_res)
  return(fsem)
  
}

res <- matrix(1, 625*1, 11)
colnames(res) <- c('M', 'est_a', 'est_b', 'est_cp', 'est_ab', 'est_total', 'p_a', 'p_b', 'p_cp', 'p_ab', 'p_total')
# colnames(diabetes_sample)
DATA <- DF
DATA$X <- DF$diagnosis_result_T2D
DATA$Y <- DF$diagnosis_result
# colnames(Xx)

for (i in 2: 626) {
  i=2
  test <- med_func(DATA, i)
  re<-summary(test)
  re
  
  res[i-1, ]<-c(colnames(DATA)[i], re$pe$est[1], re$pe$est[2], re$pe$est[11], re$pe$est[59], re$pe$est[60], 
                re$pe$pvalue[1], re$pe$pvalue[2], re$pe$pvalue[11], re$pe$pvalue[59], re$pe$pvalue[60])
  print(res[i-1, ])
  print(i)
}

fwrite(res, "/home/cmet-standard/data/OneDrive/T2DM-IDPS-CF-research/Summary/summary-new/Summary/mediation analysis/HbA1c_IDP_meantal_health_depression.csv", sep = ",", row.names = F)
