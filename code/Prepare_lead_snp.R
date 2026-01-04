library(TwoSampleMR)
library(LDlinkR)
library(plinkbinr)
library(ieugwasr)
library(pbmcapply)




# file_path <- "/media/cmet-admin/36fd75e3-6af2-4beb-8b69-54232857271c/merged_glycemic/coloc"
#
# file_list <- list.files(file_path)
#
# for (i in 1:length(file_list)) {
#
#   # i=1
#
#   dt_diseases <- fread(paste0(file_path, "/"  ,file_list[i]))
#
#   # colnames(dt_diseases)<-c("chr_38","pos_38","ref_allele","alt_allele","snp","nearest_genes","pval","mlogp",
#   #                          "beta","sebeta","af_alt","af_alt_cases","af_alt_controls", "se","Z","maf",
#   #                          "samplesize","chr","pos")
#
#   dt_diseases    <- drop_palindromic(dt_diseases)
#
#   # colnames(dt_eqtl)
#   dt_diseases_sig <- dt_diseases[dt_diseases$pval<5e-8, ]
#
#   dat <- data.frame(rsid = dt_diseases_sig$snp,
#                     beta = dt_diseases_sig$beta,
#                     se = dt_diseases_sig$se,
#                     effect_allele = dt_diseases_sig$alt_allele,
#                     other_allele = dt_diseases_sig$ref_allele,
#                     pval = dt_diseases_sig$pval,
#                     chr = dt_diseases_sig$chr,
#                     pos = dt_diseases_sig$pos)
#
#
#   dat_clumped <- ld_clump( dat = dat,
#                            clump_kb = 10000,
#                            clump_r2 = 0.001,
#                            clump_p = 5e-8,
#                            plink_bin = get_plink_exe(),
#                            bfile = "/media/cmet-admin/d087db17-fa04-4471-b171-8138fae36357/Onedriver_boxing/T2DM-IDPS-CF-research/Summary_final_v1/EUR/EUR/EUR" )
#
#   write.table(dat_clumped, paste0("/media/cmet-admin/36fd75e3-6af2-4beb-8b69-54232857271c/Multi_coloc_data/Diseases/", file_list[i]), sep = ",", row.names = F)
# }





drop_palindromic <- function(dt, maf_low=0.45, maf_high=0.55) {
  pal <- (dt$alt_allele=="A" & dt$ref_allele=="T") |
    (dt$alt_allele=="T" & dt$ref_allele=="A") |
    (dt$alt_allele=="C" & dt$ref_allele=="G") |
    (dt$alt_allele=="G" & dt$ref_allele=="C")
  ok <- !(pal & !is.na(dt$maf) & dt$maf >= maf_low & dt$maf <= maf_high)
  dt[ok]
}



# 必要包
library(data.table)
library(pbmcapply)
# ld_clump 所在包（示例 TwoSampleMR；按你的实际来源加载）
# library(TwoSampleMR)



file_path <- "/media/cmet-admin/36fd75e3-6af2-4beb-8b69-54232857271c/eQTL/Brain_cell/Prepared"
out_dir   <- "/media/cmet-admin/d087db17-fa04-4471-b171-8138fae36357/Onedriver_boxing/T2DM-IDPS-CF-research/Summary_final_v1/Data/Multi_coloc_data/Brain_cell_eQTL"
plink_bfile <- "/media/cmet-admin/d087db17-fa04-4471-b171-8138fae36357/Onedriver_boxing/T2DM-IDPS-CF-research/Summary_final_v1/EUR/EUR/EUR"
# file_sig_path <- "/media/cmet-admin/d087db17-fa04-4471-b171-8138fae36357/Onedriver_boxing/T2DM-IDPS-CF-research/Summary_final_v1/Multi_coloc/Inverse_IDP"

# # 列出所有 CSV 文件
# file_list_IDP <- list.files(path = file_sig_path, pattern = "\\.csv$", full.names = TRUE)
#
# file_list <- list.files(path = file_path, pattern = "\\.csv$", full.names = TRUE)
# # 读取并合并
# data_all <- do.call(rbind, lapply(file_list_IDP, read.csv))
#
# ID<-unique(data_all$Pheno)
# ID<-na.omit(ID)
# ID <- sprintf("%04d", ID)
#
# pattern <- paste(ID, collapse = "|")
# matched_files <- file_list[grepl(pattern, file_list)]

file_list <- list.files("/media/cmet-admin/36fd75e3-6af2-4beb-8b69-54232857271c/eQTL/Brain_cell/Prepared", full.names = T)

# 文件列表（直接存绝对路径更稳）
# file_list <- matched_files

# 并行核数（预留 1 核给系统）
n_cores <- max(1L, parallel::detectCores() - 5)

# 单文件处理函数
process_one <- function(infile) {
  # 让每个任务用自己的临时目录，避免并发读写冲突
  old_tmp <- tempdir()
  new_tmp <- tempfile(pattern = "pbmc_")
  dir.create(new_tmp, showWarnings = FALSE, recursive = TRUE)
  withr::local_tempdir(new_tmp)  # 需要 withr 包；如果没装，也可用 Sys.setenv(TMPDIR=...)

  # 输出文件名与路径
  outfile <- file.path(out_dir, basename(infile))

  tryCatch({
    # 读取
    dt_diseases <- fread(infile)

    # 统一列名
    # setnames(dt_diseases, c("chr_38","pos_38","ref_allele","alt_allele","snp","nearest_genes","pval","mlogp",
    #                         "beta","sebeta","af_alt","af_alt_cases","af_alt_controls","se","Z","maf",
    #                         "samplesize","chr","pos"))

    # 去除回文位点（假设你已有该函数）
    dt_diseases <- drop_palindromic(dt_diseases)

    # 筛选显著
    dt_sig <- dt_diseases[pval < 5e-8]
    if (nrow(dt_sig) == 0L) {
      # 没有显著位点也写个空表，保持可追踪
      fwrite(dt_sig, outfile, sep = ",")
      return(list(file = basename(infile), n_sig = 0L, n_out = 0L, status = "no hits"))
    }

    # 准备 clump 输入
    dat <- data.frame(
      rsid = dt_sig$snp,
      beta = dt_sig$beta,
      se   = dt_sig$se,
      effect_allele = dt_sig$alt_allele,
      other_allele  = dt_sig$ref_allele,
      pval = dt_sig$pval,
      chr  = dt_sig$chr,
      pos  = dt_sig$pos,
      stringsAsFactors = FALSE
    )

    # 运行 LD clump
    dat_clumped <- ld_clump(
      dat        = dat,
      clump_kb   = 10000,
      clump_r2   = 0.001,
      clump_p    = 5e-8,
      plink_bin  = get_plink_exe(),
      bfile      = plink_bfile
    )

    # 写出
    # 用 fwrite 更快；sep="," 保持你原来的 csv
    fwrite(dat_clumped, outfile, sep = ",")

    list(file = basename(infile), n_sig = nrow(dt_sig), n_out = nrow(dat_clumped), status = "ok")
  }, error = function(e) {
    # 出错也返回简要信息
    list(file = basename(infile), n_sig = NA_integer_, n_out = NA_integer_,
         status = paste0("error: ", conditionMessage(e)))
  })
}

# 创建输出目录
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 并行执行（带进度条）
res_list <- pbmclapply(file_list, process_one, mc.cores = n_cores)

# 汇总结果
res <- rbindlist(lapply(res_list, as.data.table), fill = TRUE)
print(res)
