
# ##############################################################################

#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# 0. 前期设置 -------------------------------------------------------------
library(data.table)
library(ldscr)
library(parallel)

# 1. 定义 LDSC 分析函数 ---------------------------------------------------
ldscrg_analysis <- function(trait1, trait2, trait1_name, trait2_name) {
  res <- ldsc_rg(
    munged_sumstats = list(
      trait1_name = trait1,
      trait2_name = trait2
    ),
    ancestry = "EUR"
  )
  res$rg[1] <- trait1_name
  res$rg[2] <- trait2_name
  as.data.frame(res$rg)
}

# 2. 路径配置 ---------------------------------------------------------------
trait1_path   <- "/media/cmet-admin/36fd75e3-6af2-4beb-8b69-54232857271c/GWAS_PA_data/MVPA.csv"
trait2_path   <- "/media/cmet-admin/36fd75e3-6af2-4beb-8b69-54232857271c/UKB_IDP_MAF2126_v2/"
output_dir   <- "/media/cmet-admin/d087db17-fa04-4471-b171-8138fae36357/Onedriver_boxing/ShanDong Menglu/PA/results/LDSC/"

# 确保输出目录存在
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 3. 枚举所有 trait2 文件 --------------------------------------------------
trait2_files <- list.files(trait2_path, pattern="\\.csv$", full.names=TRUE)
# 只跑前 10 个做测试
# trait2_files <- trait2_files[1:10]

# 4. 并行参数配置 -----------------------------------------------------------
num_cores <- max(1, detectCores() - 70)

# 5. 启动 PSOCK 集群 --------------------------------------------------------
cl <- makeCluster(num_cores, type = "PSOCK")
clusterEvalQ(cl, {
  library(data.table)
  library(ldscr)
})
# 把函数和路径都传给 worker
clusterExport(cl, c("ldscrg_analysis", "trait1_path", "output_dir"), envir = environment())

# 6. 并行处理：每跑一个就写出一个 CSV ---------------------------------------
parLapply(cl, trait2_files, function(fn) {
  # 6.1 读主性状
  trait1 <- fread(
    trait1_path,
    select = c("snp","alt_allele","ref_allele","samplesize","zscore"),
    col.names = c("SNP","A1","A2","N","Z")
  )

  # 6.2 读次性状
  trait2 <- tryCatch({
    fread(
      fn,
      select = c("snp","alt_allele","ref_allele","samplesize","zscore"),
      col.names = c("SNP","A1","A2","N","Z")
    )
  }, error = function(e){
    message("❌ 读文件失败：", fn, " -> ", e$message)
    return(NULL)
  })
  if (is.null(trait2)) return(NULL)

  # 6.3 运行遗传相关性
  trait2_name <- gsub("\\.csv$", "", basename(fn))
  res <- tryCatch({
    ldscrg_analysis(trait1, trait2, "MVPA", trait2_name)
  }, error = function(e){
    message("❌ LDSC 计算失败：", fn, " -> ", e$message)
    return(NULL)
  })
  if (is.null(res)) return(NULL)

  # 6.4 把结果写成单个 CSV
  out_file <- file.path(output_dir, paste0(trait2_name, "_rg.csv"))
  fwrite(res, out_file, sep = ",")
  message("✅ 写出结果：", out_file)

  # 6.5 清理内存
  rm(trait1, trait2, res)
  gc()

  invisible(TRUE)
})

# 7. 收尾 --------------------------------------------------------------
stopCluster(cl)
message("🎉 所有子任务完成，单文件结果保存在：", output_dir)
