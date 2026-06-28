# Makefile — equivalent targets for the maxillary-regeneration dry-lab pipeline.
# Each target runs one Step script with Rscript. Thresholds live in config/params.R;
# the dataset registry in config/datasets.tsv. No result is computed here — the
# Step scripts compute everything on the user's downloaded data.
#
#   make install     # install R/Bioconductor dependencies (env/install_packages.R)
#   make step0        # auto-propose comparison arms (then HAND-CURATE — see below)
#   make step1 ...    # quantitative stages; blocked until curation files exist
#   make all          # step0 .. step11 in order (stops if curation missing)
#   make clean        # remove generated outputs (keeps 00_data cache & curated files)
#   make distclean    # also remove 00_data download cache
#
# After `make step0` you MUST hand-curate:
#   01_curation/sample_arms_final.tsv  and  01_curation/rra_usable_datasets.tsv
# before any step>=1 will run.

RSCRIPT ?= Rscript
CURATION := 01_curation/sample_arms_final.tsv 01_curation/rra_usable_datasets.tsv

.PHONY: all install step0 step1 step2 step3 step4 step5 step6 step7 step8 step9 step10 step11 \
        curation-check clean distclean help

help:
	@echo "Targets: install | step0..step11 | all | clean | distclean"
	@echo "Run 'make step0', hand-curate $(CURATION), then 'make all' (or 'make step1')."

install:
	$(RSCRIPT) env/install_packages.R

# ---- Step0: auto-proposal only (human curation follows) --------------------
step0:
	$(RSCRIPT) R/step00_curation.R

# ---- curation guard: step>=1 require the two hand-curated files ------------
curation-check:
	@for f in $(CURATION); do \
	  if [ ! -f "$$f" ]; then \
	    echo "ACTION REQUIRED: missing $$f"; \
	    echo "  Run 'make step0', then hand-curate final_group/arm_clean and run validate_curation()."; \
	    exit 1; \
	  fi; \
	done

step1: curation-check
	$(RSCRIPT) R/step01_per_dataset_DE.R
step2: curation-check
	$(RSCRIPT) R/step02_rra_meta.R
step3: curation-check
	$(RSCRIPT) R/step03_wgcna_modules.R
step4: curation-check
	$(RSCRIPT) R/step04_success_projection.R
step5: curation-check
	$(RSCRIPT) R/step05_failure_projection.R
step6: curation-check
	$(RSCRIPT) R/step06_scrna_cellchat.R
step7: curation-check
	$(RSCRIPT) R/step07_mr_coloc.R
step8: curation-check
	$(RSCRIPT) R/step08_position_context.R
step9: curation-check
	$(RSCRIPT) R/step09_cmap_repurposing.R
step10: curation-check
	$(RSCRIPT) R/step10_integration.R
step11: curation-check
	$(RSCRIPT) R/step11_limitations.R

# ---- full run (mirrors run_all.R ordering) ---------------------------------
all: step0 step1 step2 step3 step4 step5 step6 step7 step8 step9 step10 step11

# ---- cleanup ---------------------------------------------------------------
# clean keeps the GEO download cache (00_data) and your curated arm files.
clean:
	rm -rf 02_per_dataset_DE/* 03_rra_meta/* 04_modules_wgcna/* \
	       05_projection/* 06_scrna/* 07_mr_coloc/* 08_cmap/* \
	       09_integration/* figures/* supp/sessionInfo.txt
	rm -f 01_curation/sample_arms_template.tsv

distclean: clean
	rm -rf 00_data/*
