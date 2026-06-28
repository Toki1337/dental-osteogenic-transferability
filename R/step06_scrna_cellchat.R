# R/step06_scrna_cellchat.R
# Step 6 — Single-cell "success vs failure" contrast + cell-cell communication rewiring.
# 失败侧主轴 = 鼠下颌特异 scRNA(对照 vs 病损骨髓/基质):
#   GSE295106  (control vs BRONJ mandibular marrow)  -> PRIMARY (有对照)
#   GSE269255  (irradiated mandibular stromal cells, ORNJ) -> 辅
#   GSE303003  (human MRONJ granulation, n=2 scRNA)   -> 仅示意 (n 极小,不做定量推断)
# 流程:read_sc -> QC -> Normalize/HVG/Scale/PCA -> Harmony(按 sample/condition 整合) ->
#       FindClusters(res) -> UMAP -> canonical-marker 注释 -> 群比例 ctrl vs failure ->
#       CellChat 分组建对象 -> compareInteractions / 关键通路(RANKL/TNF/SPP1)成败差异。
#
# Output:
#   06_scrna/sc_celltype_props.tsv        群比例 ctrl vs failure(标注样本量限制)
#   06_scrna/cellchat_pathway_diff.tsv    通路级 ctrl vs failure 通讯强度差
#   figures/sc_umap_*.pdf, figures/cellchat_*.pdf
#
# !! 重要声明:scRNA 在本项目仅作机制假设(hypothesis-generating),NOT 队列级推断。
#    样本/动物数极小(每组少数 GSM/动物),群比例差异仅为示意,不声称统计学因果。
#    GSE303003 n=2 仅作人源示意,绝不做定量队列结论。跨物种比较一律先 to_human_symbols
#    做鼠->人同源,再在人源 symbol 空间上比对/取交集。

source("R/utils.R"); init_pipeline()
suppressMessages({
  library(Seurat)
  library(ggplot2)
})
set.seed(PARAMS$seed)

# 失败侧 scRNA 计划(主轴 + 辅 + 示意)。condition 取 control / failure。
# 一个 GSE 可能既含对照又含病损(混合在 supplementary 内),由 read_sc 解析样本元数据时切分;
# 若无法切分,则按 dataset 级 default_condition 处理(并在比例分析中标注)。
SC_PLAN <- data.frame(
  accession        = c("GSE295106",   "GSE269255",   "GSE303003"),
  role             = c("failure",     "failure",     "failure"),
  organism         = c("mouse",       "mouse",       "human"),
  tier             = c("primary",     "secondary",   "illustrative"),
  default_condition= c(NA,            "failure",     "failure"),  # GSE269255/303003 多为病损侧
  stringsAsFactors = FALSE
)

# canonical markers(人源 symbol 空间;鼠数据先同源化后用同一套)
MARKERS <- list(
  Osteo_MSC   = c("RUNX2", "SP7", "COL1A1", "LEPR", "BGLAP", "PDGFRA"),
  Osteoclast  = c("CTSK", "ACP5", "MMP9", "OSCAR"),
  Macrophage  = c("CD68", "LYZ", "ADGRE1", "ITGAM"),
  Endothelial = c("PECAM1", "CDH5", "EMCN"),
  Immune_Tcell= c("PTPRC", "CD3E", "CD8A")
)

# 关注的通讯通路(成败重构重点;CellChat pathway 名)
FOCUS_PATHWAYS <- c("RANKL", "TNF", "SPP1", "TGFb", "CXCL", "WNT", "BMP", "PDGF")

# ============================================================================
# read_sc — 读取单个 scRNA 数据集为 Seurat 对象(或 NULL 跳过,绝不假装有数据)
#   优先级:10x-style supplementary (matrix.mtx[.gz] + barcodes + features/genes)
#           或 *.h5 (10x HDF5) ；否则处理后的表达矩阵 (cell x gene 或 gene x cell 的 txt/csv/tsv)。
#   找不到任何可读 supplementary -> .log WARN 并返回 NULL。
# ============================================================================
read_sc <- function(accession, role) {
  dir <- file.path(PARAMS$dir_data, role)
  if (!dir.exists(dir)) {
    .log(accession, ": data dir ", dir, " absent; download supplementary first -> SKIP", level = "WARN")
    return(NULL)
  }
  # --- 1) 10x HDF5 (*.h5) ---
  h5 <- list.files(dir, pattern = paste0("^", accession, ".*\\.h5$"), full.names = TRUE, ignore.case = TRUE)
  if (length(h5)) {
    obj <- tryCatch({
      counts <- Read10X_h5(h5[1])
      if (is.list(counts)) counts <- counts[["Gene Expression"]] %||% counts[[1]]
      CreateSeuratObject(counts, project = accession, min.cells = 3, min.features = PARAMS$sc_min_genes)
    }, error = function(e) { .log(accession, " .h5 read failed: ", conditionMessage(e), level = "WARN"); NULL })
    if (!is.null(obj)) { .log(accession, ": loaded from 10x .h5 (", ncol(obj), " cells)"); return(obj) }
  }
  # --- 2) 10x-style triplet (matrix.mtx + barcodes + features/genes) ---
  mtx <- list.files(dir, pattern = "matrix\\.mtx(\\.gz)?$", full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
  mtx <- mtx[grepl(accession, mtx, ignore.case = TRUE)]
  if (length(mtx)) {
    obj <- tryCatch({
      subdir <- dirname(mtx[1])
      counts <- Read10X(subdir)                         # 自动匹配 barcodes/features
      if (is.list(counts)) counts <- counts[["Gene Expression"]] %||% counts[[1]]
      CreateSeuratObject(counts, project = accession, min.cells = 3, min.features = PARAMS$sc_min_genes)
    }, error = function(e) { .log(accession, " 10x triplet read failed: ", conditionMessage(e), level = "WARN"); NULL })
    if (!is.null(obj)) { .log(accession, ": loaded from 10x triplet (", ncol(obj), " cells)"); return(obj) }
  }
  # --- 3) processed dense matrix (gene x cell 或 cell x gene 的 txt/csv/tsv[.gz]) ---
  flat <- list.files(dir, pattern = paste0(accession, ".*(matrix|expr|counts|umi|tpm).*\\.(txt|tsv|csv)(\\.gz)?$"),
                     full.names = TRUE, ignore.case = TRUE)
  if (length(flat)) {
    obj <- tryCatch({
      m <- as.matrix(data.table::fread(flat[1]), rownames = 1)
      storage.mode(m) <- "numeric"
      # 启发式定向:基因数通常 > 细胞数;若行更像细胞则转置
      if (nrow(m) < ncol(m) && nrow(m) < 5000 && ncol(m) > 5000) m <- t(m)
      CreateSeuratObject(m, project = accession, min.cells = 3, min.features = PARAMS$sc_min_genes)
    }, error = function(e) { .log(accession, " flat-matrix read failed: ", conditionMessage(e), level = "WARN"); NULL })
    if (!is.null(obj)) { .log(accession, ": loaded from processed matrix (", ncol(obj), " cells)"); return(obj) }
  }
  .log(accession, ": no readable scRNA supplementary found in ", dir, " -> SKIP (not fabricating data)", level = "WARN")
  NULL
}

# ============================================================================
# annotate_condition — 从 barcode 前缀 / sample 标签推断 control vs failure
#   GEO 10x supplementary 常以 GSM/sample 前缀编码 barcode(e.g. "Ctrl_AAACG-1")。
#   能解析则按关键词分组;不能解析则回退 dataset 级 default_condition。
# ============================================================================
annotate_condition <- function(obj, default_condition = NA) {
  cn <- colnames(obj)
  ctrl_kw <- "ctrl|control|sham|veh|wt|wild|normal|healthy|pre|day0|d0|untreat"
  fail_kw <- "bronj|mronj|ornj|onj|bp|zol|alendr|irradiat|radiat|disease|lesion|fail|necro"
  lab <- rep(NA_character_, length(cn)); lc <- tolower(cn)
  lab[grepl(ctrl_kw, lc)] <- "control"
  lab[grepl(fail_kw, lc)] <- "failure"
  if (all(is.na(lab))) lab[] <- default_condition %||% NA
  # sample id:barcode 前缀(去掉末尾 -1 等),用于 Harmony 批次/比例统计
  samp <- sub("[_.-][ACGTN]{6,}.*$", "", cn); samp[samp == cn] <- obj$orig.ident[samp == cn]
  obj$condition <- lab
  obj$sample    <- ifelse(is.na(samp) | samp == "", as.character(obj$orig.ident), samp)
  obj
}

# ============================================================================
# preprocess_sc — QC + Normalize + HVG + Scale + PCA(单个对象)
# ============================================================================
preprocess_sc <- function(obj, organism = "mouse") {
  mt_pat <- if (organism == "mouse") "^mt-" else "^MT-"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pat)
  obj <- subset(obj, subset = nFeature_RNA >= PARAMS$sc_min_genes & percent.mt < PARAMS$sc_max_mt_pct)
  if (ncol(obj) < 50) { .log(obj@project.name, ": <50 cells after QC -> SKIP", level = "WARN"); return(NULL) }
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  npcs <- min(PARAMS$sc_n_pcs, ncol(obj) - 1)
  obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
  obj
}

# ============================================================================
# integrate_cluster — Harmony 跨样本整合 + 邻域图 + 聚类 + UMAP
#   Harmony 不可用则回退 PCA(标注 WARN)。
# ============================================================================
integrate_cluster <- function(obj) {
  npcs <- min(PARAMS$sc_n_pcs, ncol(Embeddings(obj, "pca")))
  reduction <- "pca"
  if (length(unique(obj$sample)) > 1 && requireNamespace("harmony", quietly = TRUE)) {
    obj <- tryCatch({
      harmony::RunHarmony(obj, group.by.vars = "sample", reduction.use = "pca",
                          dims.use = 1:npcs, verbose = FALSE)
    }, error = function(e) { .log("Harmony failed -> PCA: ", conditionMessage(e), level = "WARN"); obj })
    if ("harmony" %in% Reductions(obj)) reduction <- "harmony"
  } else {
    .log("harmony unavailable or single sample -> cluster on PCA", level = "WARN")
  }
  obj <- FindNeighbors(obj, reduction = reduction, dims = 1:npcs, verbose = FALSE)
  obj <- FindClusters(obj, resolution = PARAMS$sc_cluster_res, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = reduction, dims = 1:npcs, verbose = FALSE)
  obj
}

# ============================================================================
# annotate_celltypes — 用 canonical marker 平均表达给 cluster 打标签
#   鼠数据先 to_human_symbols 同源化(在人源 symbol 空间上评分)。
# ============================================================================
annotate_celltypes <- function(obj, organism = "mouse") {
  # 取归一化数据矩阵;鼠 -> 人同源,使行名与 MARKERS 同空间
  expr <- GetAssayData(obj, slot = "data")
  rn <- rownames(expr)
  if (organism == "mouse") {
    hs <- to_human_symbols(rn, from = "mouse")
    keep <- !is.na(hs) & hs != ""
    expr <- expr[keep, , drop = FALSE]
    rn2 <- hs[keep]
    expr <- expr[!duplicated(rn2), , drop = FALSE]   # 同源后去重
    rownames(expr) <- rn2[!duplicated(rn2)]
  }
  clusters <- Idents(obj)
  score_mat <- sapply(names(MARKERS), function(ct) {
    g <- intersect(MARKERS[[ct]], rownames(expr))
    if (!length(g)) return(rep(0, nlevels(clusters)))
    sub <- expr[g, , drop = FALSE]
    tapply(colMeans(as.matrix(sub)), clusters, mean)[levels(clusters)]
  })
  rownames(score_mat) <- levels(clusters)
  # 每个 cluster 取得分最高的细胞类型;全 0 则 Unassigned
  lab <- apply(score_mat, 1, function(r) if (all(r == 0)) "Unassigned" else names(MARKERS)[which.max(r)])
  obj$cell_type <- lab[as.character(Idents(obj))]
  obj
}

# ============================================================================
# 主循环:逐数据集处理,合表;同时保留 per-dataset Seurat 供 CellChat
# ============================================================================
sc_objs   <- list()        # 处理后的 Seurat 对象(用于 CellChat / 群比例)
prop_rows <- list()         # 群比例长表

for (i in seq_len(nrow(SC_PLAN))) {
  acc  <- SC_PLAN$accession[i]; role <- SC_PLAN$role[i]
  org  <- SC_PLAN$organism[i];  tier <- SC_PLAN$tier[i]
  log_step(6, sprintf("%s [%s, %s]", acc, tier, org))

  obj <- read_sc(acc, role)
  if (is.null(obj)) next                                   # 无数据则跳过,绝不补造
  obj <- annotate_condition(obj, default_condition = SC_PLAN$default_condition[i])
  obj <- preprocess_sc(obj, organism = org)
  if (is.null(obj)) next
  obj <- integrate_cluster(obj)
  obj <- annotate_celltypes(obj, organism = org)
  obj$dataset <- acc; obj$tier <- tier; obj$organism <- org
  sc_objs[[acc]] <- obj

  # UMAP 图:cell_type / condition
  p_ct <- DimPlot(obj, group.by = "cell_type", label = TRUE, repel = TRUE) +
    ggtitle(sprintf("%s — cell types (%s)", acc, tier))
  p_cd <- DimPlot(obj, group.by = "condition") + ggtitle(sprintf("%s — condition", acc))
  ggsave(file.path(PARAMS$dir_fig, paste0("sc_umap_", acc, ".pdf")), p_ct + p_cd, width = 12, height = 5)

  # 群比例(per sample x condition x cell_type),保留样本量以便标注限制
  md <- obj@meta.data
  tab <- as.data.frame(table(sample = md$sample, condition = md$condition, cell_type = md$cell_type))
  tab <- tab[tab$Freq > 0, ]
  if (nrow(tab)) {
    denom <- tapply(tab$Freq, tab$sample, sum)
    tab$prop <- tab$Freq / denom[as.character(tab$sample)]
    tab$dataset <- acc; tab$tier <- tier; tab$organism <- org
    prop_rows[[acc]] <- tab
  }
}

# ---- 群比例汇总 + 差异检验(propeller 优先,否则卡方;标注样本量限制) ----
if (length(prop_rows)) {
  props <- do.call(rbind, prop_rows)

  # 仅在 PRIMARY(GSE295106,含真对照)上做 control vs failure 差异;辅/示意只列比例
  diff_tab <- data.frame()
  prim <- props[props$tier == "primary" & props$condition %in% c("control", "failure"), ]
  if (nrow(prim) && length(unique(prim$condition)) == 2) {
    diff_tab <- tryCatch({
      if (requireNamespace("speckle", quietly = TRUE)) {
        # propeller 需 cell-level cluster/sample/group;从 PRIMARY 对象重建
        po <- sc_objs[["GSE295106"]]
        res <- speckle::propeller(clusters = po$cell_type, sample = po$sample, group = po$condition)
        res <- as.data.frame(res); res$cell_type <- rownames(res); res$test <- "propeller"; res
      } else stop("speckle absent")
    }, error = function(e) {
      .log("propeller unavailable -> per-celltype chi-square (n very small, illustrative): ",
           conditionMessage(e), level = "WARN")
      cts <- unique(prim$cell_type)
      do.call(rbind, lapply(cts, function(ct) {
        cnt_ct  <- tapply(prim$Freq[prim$cell_type == ct], prim$condition[prim$cell_type == ct], sum)
        cnt_all <- tapply(prim$Freq, prim$condition, sum)
        m <- rbind(ct = cnt_ct[c("control","failure")], other = (cnt_all - cnt_ct)[c("control","failure")])
        m[is.na(m)] <- 0
        p <- tryCatch(suppressWarnings(stats::chisq.test(m)$p.value), error = function(e) NA_real_)
        data.frame(cell_type = ct, p_chisq = p, test = "chisq",
                   n_control = cnt_all["control"], n_failure = cnt_all["failure"])
      }))
    })
  } else {
    .log("PRIMARY GSE295106 not available with both arms -> proportions only, no diff test", level = "WARN")
  }

  # 合并:每个数据集每群每组的平均比例 + (若有)PRIMARY 差异检验
  summ <- aggregate(prop ~ dataset + tier + organism + condition + cell_type, data = props, FUN = mean)
  if (nrow(diff_tab)) summ <- merge(summ, diff_tab, by = "cell_type", all.x = TRUE)
  summ$caveat <- "hypothesis-generating; tiny n (few GSM/animals); NOT cohort-level inference"
  save_tsv(summ, file.path(PARAMS$dir_scrna, "sc_celltype_props.tsv"))
  .log("sc_celltype_props.tsv written (", nrow(summ), " rows across ", length(prop_rows), " datasets)")
} else {
  .log("no scRNA dataset yielded cells -> sc_celltype_props.tsv NOT written (no fabricated data)", level = "WARN")
}

# ============================================================================
# CellChat — 分别建对照组 / 失败组对象,比较通讯重构
#   只在 PRIMARY(GSE295106,鼠;同源到人配体-受体)上做 ctrl-vs-failure 对比,
#   因其唯一同时含 control 与 failure 真臂。其余数据集 n 太小,仅作单组示意(可选)。
# ============================================================================
run_cellchat_group <- function(obj, cells, organism = "mouse") {
  sub <- subset(obj, cells = cells)
  data_input <- GetAssayData(sub, slot = "data")        # 归一化表达
  if (organism == "mouse") {
    hs <- to_human_symbols(rownames(data_input), from = "mouse")  # 鼠->人,与 human DB 对齐
    keep <- !is.na(hs) & hs != ""
    data_input <- data_input[keep, , drop = FALSE]
    hs <- hs[keep]; data_input <- data_input[!duplicated(hs), , drop = FALSE]
    rownames(data_input) <- hs[!duplicated(hs)]
  }
  meta <- data.frame(group = sub$cell_type, row.names = colnames(sub))
  cc <- CellChat::createCellChat(object = data_input, meta = meta, group.by = "group")
  cc@DB <- CellChat::CellChatDB.human                   # 统一用 human DB(已同源)
  cc <- CellChat::subsetData(cc)
  cc <- CellChat::identifyOverExpressedGenes(cc)
  cc <- CellChat::identifyOverExpressedInteractions(cc)
  cc <- CellChat::computeCommunProb(cc, type = "triMean")
  cc <- CellChat::filterCommunication(cc, min.cells = 10)
  cc <- CellChat::computeCommunProbPathway(cc)
  cc <- CellChat::aggregateNet(cc)
  cc
}

# ============================================================================
# scrna_dysregulated_genes.tsv — gene-level control-vs-failure DE per cell_type.
#   Consumed by Step10 (scrna_support dimension). Canonical schema (Step10 reads
#   these exact cols): gene, cell_type, is_dysregulated(logical), abs_avg_log2FC, dataset.
#   失败侧主轴 = PRIMARY GSE295106(唯一同时含 control 与 failure 真臂)。
#   逐 cell_type 用 Seurat::FindMarkers(ident.1="failure", ident.2="control") 做
#   control-vs-failure 差异;鼠 symbol 先 to_human_symbols 同源到人,再与人源空间对齐。
#   细胞数太少的群(任一臂 < SC_DE_MIN_CELLS)跳过并 WARN。
#   is_dysregulated = (p_val_adj < 0.05 & |avg_log2FC| > 0.25);abs_avg_log2FC=|avg_log2FC|。
# ============================================================================
SCRNA_DYS_COLS <- c("gene", "cell_type", "is_dysregulated", "abs_avg_log2FC", "dataset")
SC_DE_MIN_CELLS <- 30   # per-arm minimum cells within a cell_type to attempt DE

# write_empty_scrna_dys — 空 schema(0 行)+ caveat 列,绝不假装有数据
write_empty_scrna_dys <- function() {
  empty <- data.frame(
    gene = character(0), cell_type = character(0),
    is_dysregulated = logical(0), abs_avg_log2FC = numeric(0),
    dataset = character(0), caveat = character(0),
    stringsAsFactors = FALSE
  )
  save_tsv(empty, file.path(PARAMS$dir_scrna, "scrna_dysregulated_genes.tsv"))
}

cellchat_done <- tryCatch({
  suppressMessages(library(CellChat))
  po <- sc_objs[["GSE295106"]]
  stopifnot(!is.null(po))
  cells_ctrl <- colnames(po)[which(po$condition == "control")]
  cells_fail <- colnames(po)[which(po$condition == "failure")]
  if (length(cells_ctrl) < 30 || length(cells_fail) < 30)
    stop("too few cells per arm for CellChat (need both control & failure)")

  cc_ctrl <- run_cellchat_group(po, cells_ctrl, organism = "mouse")
  cc_fail <- run_cellchat_group(po, cells_fail, organism = "mouse")
  cc_list <- list(control = cc_ctrl, failure = cc_fail)
  merged  <- mergeCellChat(cc_list, add.names = names(cc_list))
  save_rds(cc_list, file.path(PARAMS$dir_scrna, "cellchat_objs.rds"))

  # 整体交互数/强度对比图
  pdf(file.path(PARAMS$dir_fig, "cellchat_compareInteractions.pdf"), width = 8, height = 4)
  print(compareInteractions(merged, show.legend = FALSE, group = c(1, 2)) +
        compareInteractions(merged, show.legend = FALSE, group = c(1, 2), measure = "weight"))
  dev.off()

  # 网络重构(差异网络:失败 vs 对照)
  pdf(file.path(PARAMS$dir_fig, "cellchat_diff_network.pdf"), width = 9, height = 5)
  par(mfrow = c(1, 2))
  netVisual_diffInteraction(merged, weight.scale = TRUE)
  netVisual_diffInteraction(merged, weight.scale = TRUE, measure = "weight")
  dev.off()

  # ---- 通路级 ctrl vs failure 差异(RANKL/TNF/SPP1 等)----
  pw_strength <- function(cc) {
    if (is.null(cc@netP$pathways) || !length(cc@netP$pathways)) return(stats::setNames(numeric(0), character(0)))
    s <- sapply(seq_along(cc@netP$pathways), function(k) sum(cc@netP$prob[, , k], na.rm = TRUE))
    stats::setNames(s, cc@netP$pathways)
  }
  s_ctrl <- pw_strength(cc_ctrl); s_fail <- pw_strength(cc_fail)
  all_pw <- union(names(s_ctrl), names(s_fail))
  pw_diff <- data.frame(
    pathway          = all_pw,
    strength_control = as.numeric(s_ctrl[all_pw] %||% 0),
    strength_failure = as.numeric(s_fail[all_pw] %||% 0),
    stringsAsFactors = FALSE
  )
  pw_diff[is.na(pw_diff)] <- 0
  pw_diff$delta_failure_minus_control <- pw_diff$strength_failure - pw_diff$strength_control
  pw_diff$is_focus <- pw_diff$pathway %in% FOCUS_PATHWAYS
  pw_diff <- pw_diff[order(-abs(pw_diff$delta_failure_minus_control)), ]
  pw_diff$caveat <- "PRIMARY GSE295106 only (mouse, orthologized); hypothesis-generating, tiny n"
  save_tsv(pw_diff, file.path(PARAMS$dir_scrna, "cellchat_pathway_diff.tsv"))

  # focus 通路若存在,画其信号网络对比
  shared_focus <- intersect(FOCUS_PATHWAYS, intersect(cc_ctrl@netP$pathways, cc_fail@netP$pathways))
  for (pw in shared_focus) {
    pdf(file.path(PARAMS$dir_fig, paste0("cellchat_pathway_", pw, ".pdf")), width = 9, height = 5)
    par(mfrow = c(1, 2))
    tryCatch({
      netVisual_aggregate(cc_ctrl, signaling = pw, layout = "circle", title.name = paste(pw, "control"))
      netVisual_aggregate(cc_fail, signaling = pw, layout = "circle", title.name = paste(pw, "failure"))
    }, error = function(e) .log("netVisual ", pw, " failed: ", conditionMessage(e), level = "WARN"))
    dev.off()
  }
  .log("cellchat_pathway_diff.tsv written (", nrow(pw_diff), " pathways; focus present: ",
       paste(shared_focus, collapse = ","), ")")
  TRUE
}, error = function(e) {
  .log("CellChat comparison skipped: ", conditionMessage(e), level = "WARN")
  .log("-> cellchat_pathway_diff.tsv NOT written (need PRIMARY control+failure arms with enough cells)", level = "WARN")
  FALSE
})

# ============================================================================
# scrna_dysregulated_genes.tsv — per-cell_type control-vs-failure DE on PRIMARY.
#   优雅降级:PRIMARY 不存在 / 缺臂 / 各群细胞不足 -> 写空 schema + WARN,不 quit、不捏造。
# ============================================================================
scrna_dys_done <- tryCatch({
  po <- sc_objs[["GSE295106"]]
  if (is.null(po)) stop("PRIMARY GSE295106 object not available")
  if (!all(c("control", "failure") %in% unique(po$condition)))
    stop("PRIMARY GSE295106 lacks both control & failure arms")

  Idents(po) <- "condition"
  org_prim <- (po$organism[1] %||% "mouse")
  cell_types <- sort(unique(po$cell_type))
  dys_rows <- list()

  for (ct in cell_types) {
    cells_ct  <- colnames(po)[which(po$cell_type == ct)]
    sub_ct    <- subset(po, cells = cells_ct)
    n_ctrl    <- sum(sub_ct$condition == "control")
    n_fail    <- sum(sub_ct$condition == "failure")
    if (n_ctrl < SC_DE_MIN_CELLS || n_fail < SC_DE_MIN_CELLS) {
      .log("cell_type '", ct, "': too few cells per arm (control=", n_ctrl,
           ", failure=", n_fail, "; need >=", SC_DE_MIN_CELLS, ") -> SKIP", level = "WARN")
      next
    }
    mk <- tryCatch(
      Seurat::FindMarkers(sub_ct, ident.1 = "failure", ident.2 = "control",
                          group.by = "condition", logfc.threshold = 0,
                          min.pct = 0.1, verbose = FALSE),
      error = function(e) { .log("FindMarkers '", ct, "' failed: ",
                                 conditionMessage(e), level = "WARN"); NULL })
    if (is.null(mk) || !nrow(mk)) next

    genes_mouse <- rownames(mk)
    genes_human <- if (org_prim == "mouse") to_human_symbols(genes_mouse, from = "mouse")
                   else genes_mouse
    keep <- !is.na(genes_human) & genes_human != ""
    if (!any(keep)) {
      .log("cell_type '", ct, "': no genes mapped to human symbols -> SKIP", level = "WARN")
      next
    }
    df <- data.frame(
      gene            = as.character(genes_human[keep]),
      cell_type       = ct,
      abs_avg_log2FC  = abs(mk$avg_log2FC[keep]),
      p_val_adj       = mk$p_val_adj[keep],
      stringsAsFactors = FALSE
    )
    df$is_dysregulated <- (df$p_val_adj < 0.05) & (df$abs_avg_log2FC > 0.25)
    df$dataset <- "GSE295106"
    # 同源后可能多鼠基因映射到同一人源 symbol;每 (gene,cell_type) 保留最大 |log2FC|
    df <- df[order(-df$abs_avg_log2FC), , drop = FALSE]
    df <- df[!duplicated(df$gene), , drop = FALSE]
    dys_rows[[ct]] <- df[, c(SCRNA_DYS_COLS), drop = FALSE]
  }

  if (!length(dys_rows)) stop("no cell_type yielded enough cells for control-vs-failure DE")

  dys <- do.call(rbind, dys_rows)
  rownames(dys) <- NULL
  dys$caveat <- "PRIMARY GSE295106 only (mouse, orthologized); hypothesis-generating, tiny n"
  save_tsv(dys, file.path(PARAMS$dir_scrna, "scrna_dysregulated_genes.tsv"))
  .log("scrna_dysregulated_genes.tsv written (", nrow(dys), " gene x cell_type rows; ",
       sum(dys$is_dysregulated), " dysregulated; ", length(dys_rows), " cell types)")
  TRUE
}, error = function(e) {
  .log("scRNA per-cell_type DE skipped: ", conditionMessage(e), level = "WARN")
  .log("-> writing empty scrna_dysregulated_genes.tsv schema (no fabricated data)", level = "WARN")
  write_empty_scrna_dys()
  FALSE
})

dump_session()
.log("Step6 done. scRNA = mechanistic hypothesis layer ONLY (tiny n; GSE303003 n=2 illustrative). ",
     "CellChat ctrl-vs-failure rewiring computed on PRIMARY GSE295106 (orthologized).",
     level = if (length(sc_objs)) "INFO" else "WARN")
