# R/real_run_step8_position.R — REAL Step8 jaw position-specific context.
# Does the regeneration-competent program sit lower in jaw bone (osteogenic
# disadvantage)? GSE58474 (human mandible vs iliac osteoblasts), GSE30167 (mouse
# jaw/alveolar vs long bone). Output 05_projection/jaw_position_context.tsv.
.libPaths(c("F:/Rlib", .libPaths()))
suppressMessages({ library(GEOquery); library(limma); library(babelgene) })
options(timeout = 1800); dir.create("05_projection", showWarnings = FALSE)

collapse_symbol <- function(m, sym){ keep <- !is.na(sym)&sym!=""&sym!="---"; m<-m[keep,,drop=FALSE]; sym<-sym[keep]
  o<-order(rowMeans(m,na.rm=TRUE),decreasing=TRUE); m<-m[o,,drop=FALSE]; sym<-sym[o]
  m<-m[!duplicated(sym),,drop=FALSE]; rownames(m)<-sym[!duplicated(sym)]; m }
to_human <- function(genes){ ort<-babelgene::orthologs(genes=genes, species="mouse", human=TRUE)
  setNames(ort$human_symbol, ort$symbol)[genes] }

jaw_logfc <- function(acc, organism, jaw_kw, other_kw){
  g <- getGEO(acc, destdir="00_data/meta", GSEMatrix=TRUE, AnnotGPL=TRUE, getGPL=TRUE)
  es <- if(is.list(g)) g[[1]] else g; ex <- Biobase::exprs(es); fd <- Biobase::fData(es)
  tt <- tolower(paste(Biobase::pData(es)$title, Biobase::pData(es)$source_name_ch1))
  sc <- grep("^gene.?symbol$", colnames(fd), ignore.case=TRUE, value=TRUE)[1]
  if(is.na(sc)) sc <- grep("symbol", colnames(fd), ignore.case=TRUE, value=TRUE)[1]
  sym <- sub(" ?//.*$","", as.character(fd[[sc]]))
  if (max(ex,na.rm=TRUE) > 100) ex <- log2(ex+1)
  ex <- collapse_symbol(ex, sym)
  grp <- rep(NA_character_, length(tt)); grp[grepl(jaw_kw,tt)]<-"jaw"; grp[grepl(other_kw,tt)]<-"other"
  cat(acc, "groups:", paste(table(grp), collapse="/"), "(jaw/other among assigned)\n")
  keep <- !is.na(grp); ex<-ex[,keep,drop=FALSE]; grp<-grp[keep]
  if(length(unique(grp))<2) { cat(acc,"cannot assign both groups\n"); return(NULL) }
  design <- model.matrix(~factor(grp, levels=c("other","jaw")))
  fit <- eBayes(lmFit(ex, design), trend=TRUE, robust=TRUE)
  lfc <- topTable(fit, coef=2, number=Inf, sort.by="none")$logFC; names(lfc) <- rownames(ex)
  if(organism=="mouse"){ hs<-to_human(names(lfc)); keep<-!is.na(hs)&hs!=""; lfc<-lfc[keep]; names(lfc)<-hs[keep]; lfc<-lfc[!duplicated(names(lfc))] }
  lfc
}

h <- tryCatch(jaw_logfc("GSE58474","human","mandib|jaw|alveol","iliac|ilium|hip"), error=function(e){cat("GSE58474 err:",conditionMessage(e),"\n");NULL})
m <- tryCatch(jaw_logfc("GSE30167","mouse","jaw|mandib|alveol","long|femur|tibia|limb"), error=function(e){cat("GSE30167 err:",conditionMessage(e),"\n");NULL})

core <- read.delim("04_modules_wgcna/core_program.tsv")
sig  <- read.delim("03_rra_meta/regeneration_competent_signature.tsv")  # broader signature for context
genes <- unique(sig$gene); sd <- setNames(sig$direction, sig$gene)
hf <- if(!is.null(h)) h[genes] else setNames(rep(NA_real_,length(genes)),genes)
mf <- if(!is.null(m)) m[genes] else setNames(rep(NA_real_,length(genes)),genes)
combined <- ifelse(!is.na(hf), hf, mf)
dirv <- sd[genes]
consistent <- (dirv=="up" & combined<0) | (dirv=="down" & combined>0)
consistent[is.na(consistent)] <- FALSE
out <- data.frame(gene=genes, jaw_weak_direction_consistent=consistent,
                  abs_jaw_effect=ifelse(is.na(combined),0,abs(combined)),
                  jaw_log2FC_human=unname(hf), jaw_log2FC_mouse=unname(mf),
                  n_sources=(!is.na(hf))+(!is.na(mf)), in_core=genes %in% core$gene)
out <- out[order(-out$jaw_weak_direction_consistent, -out$abs_jaw_effect), ]
write.table(out, "05_projection/jaw_position_context.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("\nStep8: %d signature genes; %d direction-consistent with jaw osteogenic weakness (human=%s mouse=%s)\n",
            nrow(out), sum(out$jaw_weak_direction_consistent), !is.null(h), !is.null(m)))
cat("among CORE program:", sum(out$jaw_weak_direction_consistent & out$in_core), "/", sum(out$in_core), "consistent\n")
