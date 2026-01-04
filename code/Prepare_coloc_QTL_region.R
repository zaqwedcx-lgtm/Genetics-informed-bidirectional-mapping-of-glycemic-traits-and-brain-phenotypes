#### =========================
#### 0) Packages & parameters
#### =========================
library(data.table)
library(hyprcoloc)
library(plinkbinr)     # get_plink_exe()
library(vroom)
library(pbmcapply)     # Linux 并行 + 进度条
library(LDlinkR)
library(ieugwasr)

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



#

QTL_path <- "/media/cmet-admin/36fd75e3-6af2-4beb-8b69-54232857271c/eQTL/Brain_cell/Prepared"
