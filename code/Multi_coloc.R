# ====== 加载包 ======
suppressPackageStartupMessages({
  library(data.table)
  library(hyprcoloc)
  library(pbmcapply)   # 并行 + 进度条 (Linux/macOS 真并行; Windows 仅进度)
})

normalize_dt <- function(dt) {
  req <- c("snp","chr","pos","alt_allele","ref_allele","beta")
  miss <- setdiff(req, names(dt))
  if (length(miss)) stop("缺少必要列: ", paste(miss, collapse=", "))

  setDT(dt)
  dt[, chr := as.integer(chr)]
  dt[, pos := as.integer(pos)]
  dt[, alt_allele   := toupper(as.character(alt_allele))]
  dt[, ref_allele := toupper(as.character(ref_allele))]

  if (!"pval" %in% names(dt)) dt[, pval := NA_real_]
  if (!"se"   %in% names(dt) || all(is.na(dt$se))) {
    if (!all(is.na(dt$pval))) {
      z <- qnorm(dt$pval/2, lower.tail = FALSE)                         # 双尾 p 的 z
      dt[, se := ifelse(is.na(se), abs(beta)/pmax(z, .Machine$double.eps), se)]
    } else {
      stop("既没有 se，也没有 pval 以回推 se。")
    }
  }

  if (!"maf" %in% names(dt)) dt[, maf := NA_real_]                      # 方便后续过滤
  dt[]
}

# 1.2 去除“易混”的回文位点（A/T 或 C/G）且 MAF≈0.5（默认开启更稳）
drop_palindromic <- function(dt, maf_low=0.45, maf_high=0.55) {
  pal <- (dt$alt_allele=="A" & dt$ref_allele=="T") |
    (dt$alt_allele=="T" & dt$ref_allele=="A") |
    (dt$alt_allele=="C" & dt$ref_allele=="G") |
    (dt$alt_allele=="G" & dt$ref_allele=="C")
  ok <- !(pal & !is.na(dt$maf) & dt$maf >= maf_low & dt$maf <= maf_high)
  dt[ok]
}

# 1.3 用“chr:pos:alt:other”作为稳定键去重；p 更小（更显著）的优先
dedup_full_variant <- function(dt) {
  dt[, variant_id := sprintf("%s:%d:%s:%s", chr, pos, alt_allele, ref_allele)]
  setorder(dt, pval)                     # p 小在前（NA 自动排后）
  dt <- dt[!duplicated(variant_id)]
  dt[]
}

# 1.3 用“chr:pos:alt:other”作为稳定键去重；p 更
harmonize_four <- function(..., remove_pal = TRUE, trait_names = NULL, ref_index = 1) {
  dfs <- list(...)
  K   <- length(dfs)
  if (K < 2) stop("需要至少 2 个数据表。")

  # 预处理（与原逻辑一致）
  dfs <- lapply(dfs, normalize_dt)
  if (remove_pal) dfs <- lapply(dfs, drop_palindromic)
  dfs <- lapply(dfs, dedup_full_variant)

  # 仅按位置对齐（chr:pos）
  for (i in seq_len(K)) dfs[[i]][, locus_id := sprintf("%s:%d", chr, pos)]

  # 参考方向：第 ref_index 个数据
  ref <- dfs[[ref_index]][, .(locus_id, alt_ref = alt_allele, ref_ref = ref_allele, snp_ref = snp)]

  # 方向校正（保持原来只处理“互换”而不做互补的逻辑；并做 NA 安全）
  fix_dir <- function(x, ref) {
    x <- merge(x, ref, by = "locus_id", all = FALSE)
    # 丢掉等位缺失，避免比较产生 NA
    x <- x[!is.na(alt_allele) & !is.na(ref_allele) &
             !is.na(alt_ref)    & !is.na(ref_ref)]

    flip <- (x$alt_allele != x$alt_ref) & (x$alt_allele == x$ref_ref)
    flip[is.na(flip)] <- FALSE

    if (any(flip)) x[flip, beta := -beta]           # 翻转 β
    x[, `:=`(alt_allele = alt_ref, ref_allele = ref_ref)]  # 强制写成参考方向
    x
  }
  out <- lapply(dfs, fix_dir, ref = ref)

  # 取所有表的交集，并按参考顺序对齐
  keep <- Reduce(intersect, lapply(out, `[[`, "locus_id"))
  out  <- lapply(out, function(x) x[locus_id %in% keep])

  setkey(out[[ref_index]], locus_id)
  for (i in seq_len(K)) if (i != ref_index) {
    setkey(out[[i]], locus_id)
    out[[i]] <- out[[i]][out[[ref_index]]]
  }

  # 构造 B/S 矩阵
  B <- do.call(cbind, lapply(out, `[[`, "beta"))
  S <- do.call(cbind, lapply(out, `[[`, "se"))

  # 列名：保持兼容性（4 个输入时仍用原来的默认名）
  if (is.null(trait_names)) {
    if (K == 4) {
      trait_names <- c("eQTL","HbA1c","CortVol","AD")
    } else {
      trait_names <- paste0("Trait", seq_len(K))
    }
  }
  if (length(trait_names) != K) stop("trait_names 长度需与输入数据个数一致。")
  colnames(B) <- trait_names
  colnames(S) <- trait_names

  # SNP 标识：优先用参考数据的 rsID，否则用 locus:alt
  snp_id <- ifelse(!is.na(out[[ref_index]]$snp_ref) & out[[ref_index]]$snp_ref != "",
                   out[[ref_index]]$snp_ref,
                   sprintf("%s:%s", out[[ref_index]]$locus_id, out[[ref_index]]$alt_allele))

  list(B = B, S = S, snps = snp_id)
}

# ====== 基本配置 ======
glycemic_name <- "HbA1c"

glycemic_path     <- "/media/cmet-admin/36fd75e3-6af2-4beb-8b69-54232857271c/Multi_coloc_data/Multi_coloc_leadsnp/Glycemic_traits"
IDP_ID_path       <- "/media/cmet-admin/d087db17-fa04-4471-b171-8138fae36357/Onedriver_boxing/T2DM-IDPS-CF-research/Summary_final_v1/Multi_coloc/Forward_IDP"
IDP_forward_path  <- "/media/cmet-admin/36fd75e3-6af2-4beb-8b69-54232857271c/Multi_coloc_data/Multi_coloc_leadsnp/IDP/Forward"
disease_ID_path   <- "/media/cmet-admin/d087db17-fa04-4471-b171-8138fae36357/Onedriver_boxing/T2DM-IDPS-CF-research/Summary_final_v1/Figure/Data/heatmap_11disea_forward"
disease_path      <- "/media/cmet-admin/36fd75e3-6af2-4beb-8b69-54232857271c/Multi_coloc_data/Multi_coloc_leadsnp/Diseases"
QTL_path          <- "/media/cmet-admin/36fd75e3-6af2-4beb-8b69-54232857271c/Multi_coloc_data/Multi_coloc_leadsnp/Brain_cell_eQTL"

glycemic_region_path <- "/media/cmet-admin/36fd75e3-6af2-4beb-8b69-54232857271c/Multi_coloc_data/Multi_coloc_region/glycemic"
IDP_region_path      <- "/media/cmet-admin/36fd75e3-6af2-4beb-8b69-54232857271c/Multi_coloc_data/Multi_coloc_region/IDP"
disease_region_path  <- "/media/cmet-admin/36fd75e3-6af2-4beb-8b69-54232857271c/Multi_coloc_data/Multi_coloc_region/disease"
QTL_region_path      <- "/media/cmet-admin/36fd75e3-6af2-4beb-8b69-54232857271c/Multi_coloc_data/Multi_coloc_region/brain_cell_all"

out_path <- "/media/cmet-admin/36fd75e3-6af2-4beb-8b69-54232857271c/Multi_colco_result/brain_cell"
dir.create(file.path(out_path, glycemic_name), recursive = TRUE, showWarnings = FALSE)

# ====== 读入 IDP / disease / QTL 列表 ======
IDP <- fread(file.path(IDP_ID_path, paste0(glycemic_name, "_IDP.csv")))
IDP_pheno <- sprintf("%04d", IDP$Pheno)

disease <- fread(file.path(disease_ID_path, paste0(glycemic_name, ".csv")))
disease_name <- disease[pval < 0.05, outcome]

QTL_name <- sub("\\.csv$", "", list.files(QTL_path, pattern = "\\.csv$", full.names = FALSE))

# ====== 组合表 ======
trait_table <- expand.grid(
  glycemic = glycemic_name,
  IDP      = IDP_pheno,
  disease  = disease_name,
  QTL      = QTL_name,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# 为避免 data.frame 子集变原子向量，把每行拆成列表元素
rows <- split(trait_table, seq_len(nrow(trait_table)))

# ====== 并行核心数 ======
cores <- max(1, parallel::detectCores() - 50)

# ====== 工作函数：处理一行组合 ======
run_one_row <- function(traits) {
  # traits 是一行 data.frame
  g   <- as.character(traits[["glycemic"]])
  idp <- as.character(traits[["IDP"]])
  dis <- as.character(traits[["disease"]])
  qtl <- as.character(traits[["QTL"]])

  safe_read <- function(path) tryCatch(fread(path), error = function(e) NULL)

  # 读取四类 leadSNP 列表
  glycemic_leadsnp <- safe_read(file.path(glycemic_path,    paste0(g,   ".csv")))
  IDP_leadsnp      <- safe_read(file.path(IDP_forward_path, paste0(idp, ".csv")))
  disease_leadsnp  <- safe_read(file.path(disease_path,     paste0(dis, ".csv")))
  QTL_leadsnp      <- safe_read(file.path(QTL_path,         paste0(qtl, ".csv")))

  # 任一缺失则跳过
  if (is.null(glycemic_leadsnp) || is.null(IDP_leadsnp) || is.null(disease_leadsnp) || is.null(QTL_leadsnp)) {
    message("[skip] leadsnp list missing for: ", paste(c(g,idp,dis,qtl), collapse = " | "))
    return(invisible(NULL))
  }

  leadsnp_all <- unique(c(glycemic_leadsnp$rsid,
                          IDP_leadsnp$rsid,
                          disease_leadsnp$rsid,
                          QTL_leadsnp$rsid))

  trait_cor      <- diag(4)
  sample_overlap <- diag(4)

  # 逐 SNP 处理（顺序即可，避免并发写文件冲突）
  for (j in seq_along(leadsnp_all)) {
    tag <- leadsnp_all[j]

    sub_glycemic <- safe_read(file.path(glycemic_region_path, paste0(tag, "_", g,   ".csv")))
    sub_IDP      <- safe_read(file.path(IDP_region_path,      paste0(tag, "_", idp, ".csv")))
    sub_disease  <- safe_read(file.path(disease_region_path,  paste0(tag, "_", dis, ".csv")))
    sub_QTL      <- safe_read(file.path(QTL_region_path,      paste0(tag, "_", qtl, ".csv")))

    if (is.null(sub_glycemic) || is.null(sub_IDP) || is.null(sub_disease) || is.null(sub_QTL)) next
    if (min(nrow(sub_glycemic), nrow(sub_IDP), nrow(sub_disease), nrow(sub_QTL)) < 50) next

    # 协调
    mats <- tryCatch(
      harmonize_four(sub_glycemic, sub_IDP, sub_disease, sub_QTL,
                     trait_names = c(g, idp, dis, qtl),
                     ref_index   = 2, remove_pal = FALSE),
      error = function(e) NULL
    )
    if (is.null(mats)) next

    # HyPrColoc
    res <- tryCatch(
      hyprcoloc(
        effect.est     = mats$B,
        effect.se      = mats$S,
        trait.names    = c(g, idp, dis, qtl),
        snp.id         = mats$snps,
        trait.cor      = trait_cor,
        sample.overlap = sample_overlap,
        uniform.priors = FALSE,
        prior.1        = 1e-4,
        prior.c        = 0.02,
        reg.thresh     = "default",
        align.thresh   = "default",
        snpscores      = TRUE
      ),
      error = function(e) NULL
    )
    if (is.null(res) || is.null(res$results)) next

    out <- as.data.table(res$results)

    # 区间信息
    idx  <- which(sub_glycemic$snp == tag)
    chr  <- if (length(idx)) sub_glycemic$chr[idx[1]] else NA
    w0   <- sub_glycemic$pos[1]
    w1   <- sub_glycemic$pos[nrow(sub_glycemic)]

    out[, `:=`(lead_id = tag, chr = chr, w_start = w0, w_end = w1,
               n_snps = nrow(mats$B),
               glycemic = g, IDP = idp, disease = dis, QTL = qtl)]

    out_file <- file.path(out_path, glycemic_name,
                          paste0(tag, "_", g, "_", idp, "_", dis, "_", qtl, ".csv"))
    fwrite(out, out_file)
  }

  invisible(NULL)
}

# ====== 并行执行 ======
pbmcapply::pbmclapply(
  rows,
  run_one_row,
  mc.cores = cores,
  mc.preschedule = FALSE
)

message("Done.")
