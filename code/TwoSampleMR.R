library(vroom)
library(dplyr)
library(vroom)
library(TwoSampleMR)
library(LDlinkR)
library(plinkbinr)
library(dplyr)
library(ieugwasr)
library(genetics.binaRies)


# Function 
Radial_data = function(dat){
  library(dplyr)
  library(RadialMR)
  RadialMR_r<- ivw_radial(r_input = dat, alpha = 0.05,
                          weights = 1, tol = 0.0001, summary = T)
  RadialMR_r1<- egger_radial(r_input = dat, alpha = 0.05,
                             weights = 1, summary = T)
  plot_radial(c(RadialMR_r,RadialMR_r1))
  which(RadialMR_r$data$Outliers == "Outlier")
  out_data<-RadialMR_r$data[RadialMR_r$data$Outliers != "Outlier", ]
  variant<-out_data$SNP
  variant<-as.data.frame(variant)
  colnames(variant)<-"SNP"
  out<-left_join(variant, dat, by="SNP")
  return(out)
}

het_test = function(dat, exp_name="exposure", out_name="outcome"){
  out<-mr_heterogeneity(dat)[ ,c("exposure", "outcome", "method", "Q", "Q_df", "Q_pval")]
  out$exposure<- exp_name
  out$outcome<-out_name
  return(out)
} #异质性检验,就是看不同工具什么的会不会影响

ple_test = function(dat, exp_name="exposure", out_name="outcome"){
  out<-mr_pleiotropy_test(dat)[ ,c("exposure", "outcome", "egger_intercept", "se", "pval")]
  out$exposure<- exp_name
  out$outcome<-out_name
  return(out)
} #多效性检验,排除工具变量直接影响结局和混杂因素

steiger_test<-function(dat, exp_name="exposure", out_name="outcome", exposure_type = "quantitative", outcome_type = "quantitative"){
  
  if(exposure_type == "quantitative" && outcome_type == "quantitative"){
    out<-directionality_test(dat)[, c("exposure", "outcome", "snp_r2.exposure", "snp_r2.outcome", "correct_causal_direction", "steiger_pval")]
    out$exposure<- exp_name
    out$outcome<-out_name
    
  }else if(exposure_type == "binary" && outcome_type == "quantitative"){
    dat$r.exposure<-get_r_from_bsen(dat$beta.exposure, dat$se.exposure, dat$samplesize.exposure)
    out<-directionality_test(dat)[, c("exposure", "outcome", "snp_r2.exposure", "snp_r2.outcome", "correct_causal_direction", "steiger_pval")]
    out$exposure<- exp_name
    out$outcome<-out_name
    
  }else if(exposure_type == "quantitative" && outcome_type == "binary"){
    dat$r.outcome<-get_r_from_bsen(dat$beta.outcome, dat$se.outcome, dat$samplesize.outcome)
    out<-directionality_test(dat)[, c("exposure", "outcome", "snp_r2.exposure", "snp_r2.outcome", "correct_causal_direction", "steiger_pval")]
    out$exposure<- exp_name
    out$outcome<-out_name
    
  }else{
    dat$r.exposure<-get_r_from_bsen(dat$beta.exposure, dat$se.exposure, dat$samplesize.exposure)
    dat$r.outcome<-get_r_from_bsen(dat$beta.outcome, dat$se.outcome, dat$samplesize.outcome)
    out<-directionality_test(dat)[, c("exposure", "outcome", "snp_r2.exposure", "snp_r2.outcome", "correct_causal_direction", "steiger_pval")]
    out$exposure<- exp_name
    out$outcome<-out_name
  }
  return(out)
} #暴露变量是否比结局变量解释更多的遗传方差,就是谁是因谁是果

Robust_analysis = function(dat, exp_name="exposure", out_name="outcome"){
  library(mr.raps)
  a<-mr.raps(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)
  x<-matrix(1, 1, 7)
  x<-as.data.frame(x)
  colnames(x)<-c("exposure", "outcome", "beta.hat", "beta.se", "beta.p.value", "naive.se", "chi.sq.test")
  x$exposure<-exp_name
  x$outcome<-out_name
  x$beta.hat<-a$beta.hat
  x$beta.se<-a$beta.se
  x$beta.p.value<-a$beta.p.value
  x$naive.se<-a$naive.se
  x$chi.sq.test<-a$chi.sq.test
  out<-x
  return(out)
}

# write_csv = function(name){
#   write.table(main_res, "", sep = ",", row.names = F)
# }

#提取G6和F5
x<-vroom("/media/cmet-standard/MyBook/Yunfan/finngen_R12_manifest.tsv")

y<-x[grep("F5|G6", x$category), ]
y1 <-x[ x$phenocode %in% c("C_STROKE"),]
y <- rbind(y, y1)

z<-subset(y, y$num_cases>1000)
#处理manifest.tsv
# thattsv <- vroom("D:/Brain Geo/R/finngen_R12_manifest.tsv")
# thattsv <-thattsv[ thattsv$phenocode %in% c("C_STROKE","F5_ANXIETY","F5_DEPRESSIO","F5_EATING"
#                                             ,"G6_EPLEPSY","G6_MIGRAINE","G6_PARKINSON","G6_SLEEPDISOTH","G6_TIA"),]
# thattsv <-thattsv[ thattsv$phenocode %in% c("C_STROKE","G6_EPLEPSY","G6_MIGRAINE","G6_PARKINSON","G6_SLEEPDISOTH","I9_TIA"),]
# write.table(thattsv,"D:/Brain Geo/R/myset.csv",sep = ",", row.names = F)





# file_name<-vroom("/media/cmet-standard/MyBook/5820_data/Summary_level_GWAS_data/FinnGEN_dementia_R5R9/R5/dementia.csv")
# file_name<-vroom("D:/Brain Geo/R/myset.csv")
file_name<-z
# 这个到底是什么数据,这个是多种病症的一个汇总表
# file_name<-vroom("D:/Brain Geo/R/finngen_R5_F5_DEMENTIA.gz.csv")
# file_name<-file_name[file_name$phenocode %in% c("finngen_R5_C_STROKE", "finngen_R5_F5_ANXIETY", "finngen_R5_F5_DEPRESSIO", "finngen_R5_F5_EATING", "finngen_R5_G6_EPLEPSY","finngen_R5_G6_MIGRAINE","finngen_R5_G6_PARKINSON","finngen_R5_G6_SLEEPDISOTH","finngen_R5_G6_TIA","F5_DEMNAS"), ]
# file_name<-file_name[file_name$phenocode %in% c("C_STROKE","F5_ANXIETY","F5_DEPRESSIO","F5_EATING"
#                                                 ,"G6_EPLEPSY","G6_MIGRAINE","G6_PARKINSON","G6_SLEEPDISOTH","G6_TIA"), ]
# file_name<-file_name[file_name$phenocode %in% c("C_STROKE","G6_EPLEPSY","G6_MIGRAINE","G6_PARKINSON","G6_SLEEPDISOTH","I9_TIA"), ]
# path_on_mylaptop <- file_name$path_on_mylaptop

#read数据
# for (i in 1:128){
# url <- file_name$path_https[i]
# download.file(url, destfile = paste0("/media/cmet-standard/MyBook/Yunfan/R12/", file_name$phenocode[i]))
# data <- vroom(gzfile(url))
# }
# data <- read.csv(gzfile(url))
# path1 <- paste0("/media/cmet-standard/MyBook/Yunfan/R12/", file_name$phenocode[j])
# RT: 330069, occupational: 248847, CP: 257828, intelligence: 269867

het_analysis_result<-NULL
ple_analysis_result<-NULL
steiger_analysis_result<-NULL
robust_analysis_result<-NULL
presso_analysis_result<-NULL
radial_analysis_result<-NULL
exp_name<-"T2D"
# exp_name<-"HbA1c"
# exp_name<-"FI_EUR"
# exp_name<-"FG_EUR"
# exp_name<-"2hGlu"



# exp<-vroom("/media/cmet-standard/MyBook/5820_data/Summary_level_GWAS_data/T2D_related_traits/summary_sig/T2D_sig.csv")
# exp<-vroom("/media/cmet-standard/MyBook/Yunfan/Glycemic_GWAS_exposure/CHEN_MAGIC1000G_2hGlu_EUR_sig.csv")
# exp<-vroom("/media/cmet-standard/MyBook/Yunfan/Glycemic_GWAS_exposure/CHEN_MAGIC1000G_FG_EUR_sig.csv")
# exp<-vroom("/media/cmet-standard/MyBook/Yunfan/Glycemic_GWAS_exposure/CHEN_MAGIC1000G_FI_EUR_sig.csv")
# exp<-vroom("/media/cmet-standard/MyBook/Yunfan/Glycemic_GWAS_exposure/CHEN_MAGIC1000G_HbA1c_EUR_sig.csv")
exp<-vroom("/media/cmet-standard/MyBook/Yunfan/Glycemic_GWAS_exposure/T2D_sig.csv")


exposure_dat<-exp

exposure_dat$samplesize_col<-exposure_dat$sample_size

#在mr分析之前要生成两个标准文件作为输入,一个是标准格式的exposure文件,一个是标准格式的outcome文件
exposure_dat <- format_data(as.data.frame(exposure_dat), 
                            type = "exposure", snps = NULL, 
                            snp_col = "snp", beta_col = "beta", se_col = "se", 
                            effect_allele_col = "effect_allele", 
                            other_allele_col = "other_allele", pval_col = "pvalue", 
                            samplesize_col = "samplesize_col")

exp_dat<-exposure_dat



exp_dat <- exp_dat %>% 
  rename(
    rsid = SNP,
    pval = pval.exposure
  )

#去掉连锁不平衡,筛选显著位点
exp_dat_clumped <- ld_clump(
  dat = exp_dat,
  # clump_kb = 10000,
  clump_kb = 1000000000000,
  clump_r2 = 0.001, 
  clump_p = 5e-8,
  # plink_bin = get_plink_binary(),
  plink_bin = get_plink_exe(),
  # bfile = "/media/cmet-standard/MyBook/5820_data/local_ldsc_reference/LD_reference_panel/EUR/EUR" #path to LD reference dataset
  bfile = "/media/cmet-standard/MyBook/Yunfan/EUR/EUR/EUR" #path to LD reference dataset
)

exp_dat_clumped <- exp_dat_clumped %>% 
  rename(
    SNP = rsid,
    pval.exposure = pval
  )

print(paste0("Number of IVs: ", as.character(length(exp_dat_clumped$SNP))))


for (i in 1:128) {
  
  # out_path=paste("/media/cmet-standard/MyBook/5820_data/Summary_level_GWAS_data/FinnGEN_dementia_R5R9/R5/summary_all/", "finngen_R5_",file_name$phenocode[i], ".gz.csv", sep = "")
  # out_path=paste("D:/Brain Geo/R/finngen_R5_F5_DEMENTIA.gz.csv")
  out_path=paste0("/media/cmet-standard/MyBook/Yunfan/R12/", file_name$phenocode[i]) #换自己的
  # out_path=path_on_mylaptop[i] #换自己的
  
  # out_path<-paste0("/home/cmet-standard/data/Summary_level_GWAS_data/FinnGEN_dementia_R5R9/R5/summary_all/", "finngen_R5_", file_name$phenocode[i], ".gz.csv")
  outcome_dat<-vroom(out_path)
  outcome_dat$samplesize_col<-file_name$num_cases[i]+file_name$num_controls[i]
  # outcome_dat$samplesize_col<-file_name$N_case+file_name$N_control #新创建一列samplesize_col
  # outcome_dat$samplesize_col<- sum(outcome_dat$n_hom_cases, outcome_dat$n_het_cases, outcome_dat$n_hom_controls, outcome_dat$n_het_controls, na.rm = TRUE) #新创建一列samplesize_col
  outcome_dat$se<-sqrt((outcome_dat$beta^2) / qchisq(outcome_dat$pval, df = 1, lower.tail = FALSE)) #算se
  # outcome_dat$samplesize_col<-IDP$`N(fullscan)`[i]
  colnames(outcome_dat)
  
  # out_dat <- format_data(as.data.frame(outcome_dat), type = "outcome", 
  #                        snps = exp_dat_clumped$SNP, snp_col = "SNP", 
  #                        beta_col = "BETA", se_col = "SE", 
  #                        effect_allele_col = "A1", other_allele_col = "A2", 
  #                        pval_col = "P", samplesize_col = "samplesize_col")
  #在mr分析之前要生成两个标准文件作为输入,一个是标准格式的exposure文件,一个是标准格式的outcome文件
  out_dat <- format_data(as.data.frame(outcome_dat), type = "outcome",
                         snps = exp_dat_clumped$SNP, snp_col = "rsids",
                         beta_col = "beta", se_col = "se",
                         effect_allele_col = "alt", other_allele_col = "ref",
                         pval_col = "pval", samplesize_col = "samplesize_col")
  
  
  # pro_out<-proxy_snp(exp_dat_clumped, out_dat, pop="EUR", RR = 0.8, token = "282d875159e2",out_beta_col="BETA", out_se_col="SE", out_effect_allele_col="A1",
  # out_other_allele_col="A2", out_pval_col="P", out_samplesize = "samplesize_col", out_path)
  
  dat<-harmonise_data(exp_dat_clumped, out_dat, action = 3) #对齐效应等位基因（Effect Alleles）,剔除那些无法对齐或存在问题的 SNP 数据
  
  radial_data<-Radial_data(dat)
  
  out_name<-file_name$phenocode[i]
  # out_name<-"dementia"
  
  radial_analysis<-mr(radial_data) #RadialMR 包进行因果推断分析，并移除数据中的离群点
  radial_analysis$exposure<-exp_name
  radial_analysis$outcome<-out_name
  radial_analysis_result<-rbind(radial_analysis_result, radial_analysis)
  
  # het test
  het_analysis<-het_test(radial_data, exp_name = exp_name, out_name = out_name)
  het_analysis_result<-rbind(het_analysis_result, het_analysis)
  
  # ple test
  ple_analysis<-ple_test(radial_data, exp_name = exp_name, out_name = out_name)
  ple_analysis_result<-rbind(ple_analysis_result, ple_analysis)
  
  # steiger test
  steiger_analysis<-steiger_test(radial_data, exp_name = exp_name , out_name = out_name, exposure_type = "binary", outcome_type = "binary")
  steiger_analysis_result<-rbind(steiger_analysis_result, steiger_analysis)
  
  # Robust analysis
  robust_analysis<-Robust_analysis(radial_data, exp_name = exp_name, out_name = out_name)
  robust_analysis_result<-rbind(robust_analysis_result, robust_analysis)
  
  
  print(i)
}
#最重要的就是看Inverse variance weighted这个方法的P值
main_res<-generate_odds_ratios(radial_analysis_result)

# write.table(main_res, "/home/cmet-standard/data/OneDrive/T2DM-IDPS-CF-research/Summary/summary-new/radial_main_result/2hGlu-IDP/res/main_analysis_result.csv", sep = ",", row.names = F)
# write.table(het_analysis_result, "/home/cmet-standard/data/OneDrive/T2DM-IDPS-CF-research/Summary/summary-new/radial_main_result/2hGlu-IDP/res/het_analysis_result.csv", sep = ",", row.names = F)
# write.table(ple_analysis_result, "/home/cmet-standard/data/OneDrive/T2DM-IDPS-CF-research/Summary/summary-new/radial_main_result/2hGlu-IDP/res/ple_analysis_result.csv", sep = ",", row.names = F)
# write.table(steiger_analysis_result, "/home/cmet-standard/data/OneDrive/T2DM-IDPS-CF-research/Summary/summary-new/radial_main_result/2hGlu-IDP/res/steiger_analysis_result.csv", sep = ",", row.names = F)
# write.table(robust_analysis_result, "/home/cmet-standard/data/OneDrive/T2DM-IDPS-CF-research/Summary/summary-new/radial_main_result/2hGlu-IDP/res/robust_analysis_result.csv", sep = ",", row.names = F)

               write.table(main_res, "/media/cmet-standard/MyBook/Yunfan/R12_output_T2D_kb12_0/main_analysis_result.csv", sep = ",", row.names = F)
    write.table(het_analysis_result, "/media/cmet-standard/MyBook/Yunfan/R12_output_T2D_kb12_0/het_analysis_result.csv", sep = ",", row.names = F)
    write.table(ple_analysis_result, "/media/cmet-standard/MyBook/Yunfan/R12_output_T2D_kb12_0/ple_analysis_result.csv", sep = ",", row.names = F)
write.table(steiger_analysis_result, "/media/cmet-standard/MyBook/Yunfan/R12_output_T2D_kb12_0/steiger_analysis_result.csv", sep = ",", row.names = F)
 write.table(robust_analysis_result, "/media/cmet-standard/MyBook/Yunfan/R12_output_T2D_kb12_0/robust_analysis_result.csv", sep = ",", row.names = F)
