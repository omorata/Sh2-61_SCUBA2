#  physplots.mk
#  O. Morata 2020
#
#  rules to produce plots out of the physical parameters calculated on
#  the JCMT SCUBA2 data of Sh2-61
#
#  To be included by the main Makefile

define MapPhysParam
# Template to plot maps of the calculated physical parameter
#
# Arguments: 1- tgt, 2- findclump id, 3- physparam calculation variant,
#            4- parameter to map
#
$(eval tgt_dir := $(RES_DIR)/analysis_maps)
$(eval out_dir := $(RES_DIR)/analysis_maps)

$(eval orig_file := $(tgt_dir)/calcs-$(1)-$(2)_$(3).log)

$(eval out_file := $(out_dir)/$(SNAME)-$(1)-$(2)_$(3)-$(4)-map.pdf)
$(eval cfg_file := $(CFG_DIR)/analysis/$(SNAME)-$(1)-$(2)_$(3)-$(4)-map.yml)


$(out_file): $$(wildcard $$(cfg_file)) $(orig_file)
	@if [ -f $(cfg_file) ]; then \
	     $(EXT_DIR)/dbxmap.py \
                 -c $(cfg_file) \
                 -o $(out_dir) \
                 -w $(tgt_dir) ;\
         else \
             echo -e "\n++ Ignoring rule:\n      $(out_file)" ;\
             echo -e "    No cfg file:\n      $(cfg_file)" ;\
         fi


map-physpar-$(1)-$(2)_$(3)-$(4): $(out_file)
.PHONY: map-physpar-$(1)-$(2)_$(3)-$(4)

maps-physpar-$(1)-$(2)_$(3): map-physpar-$(1)-$(2)_$(3)-$(4)
.PHONY: maps-physpar-$(1)-$(2)_$(3)

maps-physpar-$(1)-$(2): maps-physpar-$(1)-$(2)_$(3)
.PHONY: maps-physpar-$(1)-$(2)

maps-physpar-$(1): maps-physpar-$(1)-$(2)
.PHONY: maps-physpar-$(1)

maps-physpar: maps-physpar-$(1)
.PHONY: maps-physpar


clean-map-physpar-$(1)-$(2)_$(3)-$(4):
	@rm -fv $(out_file)
.PHONY: clean-map-physpar-$(1)-$(2)_$(3)-$(4)


clean-maps-physpar-$(1)-$(2)_$(3): clean-map-physpar-$(1)-$(2)_$(3)-$(4)
.PHONY: clean-maps-physpar-$(1)-$(2)_$(3)

clean-maps-physpar-$(1)-$(2): clean-maps-physpar-$(1)-$(2)_$(3)
.PHONY: clean-maps-physpar-$(1)-$(2)

clean-maps-physpar-$(1): clean-maps-physpar-$(1)-$(2)
.PHONY: clean-maps-physpar-$(1)

clean-maps-physpar: clean-maps-physpar-$(1)
.PHONY: clean-maps-physpar-$(1)

endef



define HistoPlots
# Template to plot histograms from the pixel maps of calculated physical
# parameters
#
# Arguments: 1- tgt, 2- findclump id, 3- physcalc tag, 4-histoplot tag
#
# The plots depend on calcs.log as a generic way  and on the wildcard of
# a .yml file with the same stem as the pdf
#
# models
#  Sh2_61-j850r0_co_mb__j450ro_mb-fw_01_ref-histo-N_all.pdf
#

$(eval out_dir := $(RES_DIR)/analysis_maps)

$(eval orig_calc := $(outdir)/calcs-$(1)-$(2)_$(3).log)

$(eval cfg_file := $(CFG_DIR)/analysis/$(SNAME)-$(1)-$(2)_$(3)-histo-$(4).yml)
$(eval histo_file := $(out_dir)/$(SNAME)-$(1)-$(2)_$(3)-histo-$(4).pdf)

$(histo_file): $$(wildcard $$(cfg_file)) $(orig_calc)
	@if [ -f $(cfg_file) ]; then \
	     $(BIN_DIR)/mapstats.py \
                 -c $(cfg_file) \
                 -o $(out_dir) \
                 -w $(RES_DIR) ;\
         else \
             echo -e "\n++ Ignoring rule\n      $(histo_file)" ;\
             echo -e "    No cfg file\n      $(cfg_file)" ;\
         fi


histo-$(1)-$(2)_$(3)-$(4): $(histo_file)

.PHONY: histo-$(1)-$(2)_$(3)-$(4)


histo-$(1)-$(2)_$(3): histo-$(1)-$(2)_$(3)-$(4)

.PHONY: histo-$(1)-$(2)_$(3)


histo-$(1)-$(2): histo-$(1)-$(2)_$(3)

.PHONY: histo-$(1)-$(2)


histo-$(1): histo-$(1)-$(2)

.PHONY: histo-$(1)


histo: histo-$(1)

.PHONY: histo


clean-histo-$(1)-$(2)_$(3)-$(4):
	@rm -fv $(histo_file)

.PHONY: clean-histo-$(1)-$(2)_$(3)-$(4)


clean-histo-$(1)-$(2)_$(3): clean-histo-$(1)-$(2)_$(3)-$(4)

.PHONY: clean-histo-$(1)-$(2)_$(3)


clean-histo-$(1)-$(2): clean-histo-$(1)-$(2)_$(3)

.PHONY: clean-histo-$(1)-$(2)


clean-histo-$(1): clean-histo-$(1)-$(2)

.PHONY: clean-histo-$(1)


clean-histo: clean-histo-$(1)

.PHONY: clean-histo

endef



define XYPlots
# Template to plot xyplots from the pixel maps of calculated physical
# parameters
#
# Arguments: 1- tgt, 2- findclump id, 3- physcalc tag, 4-xyplot tag
#
# the plots depend on calcs.log as a generic way  and on the wildcard of
# a .yml file with the same stem as the pdf
#
# The output name must have xyplot between arg 3 and arg 4
#
# models
#  Sh2_61-j850r0_co_mb__j450ro_mb-fw_01_ref-xyplot-N_vs_tdust_all.pdf
#

$(eval out_dir := $(RES_DIR)/analysis_maps)

$(eval orig_calc := $(outdir)/calcs-$(1)-$(2)_$(3).log)

$(eval cfg_file := $(CFG_DIR)/analysis/$(SNAME)-$(1)-$(2)_$(3)-xyplot-$(4).yml)
$(eval xyplot_file := $(out_dir)/$(SNAME)-$(1)-$(2)_$(3)-xyplot-$(4).pdf)

$(xyplot_file): $$(wildcard $$(cfg_file)) $(orig_calc)
	@if [ -f $(cfg_file) ]; then \
	     $(BIN_DIR)/mapstats.py \
                 -c $(cfg_file) \
                 -o $(out_dir) \
                 -w $(RES_DIR) ;\
         else \
             echo -e "\n++ Ignoring rule\n      $(xyplot_file)" ;\
             echo -e "    No cfg file\n      $(cfg_file)" ;\
         fi


xyplot-$(1)-$(2)_$(3)-$(4): $(xyplot_file)

.PHONY: xyplot-$(1)-$(2)_$(3)-$(4)


xyplot-$(1)-$(2)_$(3): xyplot-$(1)-$(2)_$(3)-$(4)

.PHONY: xyplot-$(1)-$(2)_$(3)


xyplot-$(1)-$(2): xyplot-$(1)-$(2)_$(3)

.PHONY: xyplot-$(1)-$(2)


xyplot-$(1): xyplot-$(1)-$(2)

.PHONY: xyplot-$(1)


xyplot: xyplot-$(1)

.PHONY: xyplot


clean-xyplot-$(1)-$(2)_$(3)-$(4):
	@rm -fv $(xyplot_file)

.PHONY: clean-xyplot-$(1)-$(2)_$(3)-$(4)


clean-xyplot-$(1)-$(2)_$(3): clean-xyplot-$(1)-$(2)_$(3)-$(4)

.PHONY: clean-xyplot-$(1)-$(2)_$(3)


clean-xyplot-$(1)-$(2): clean-xyplot-$(1)-$(2)_$(3)

.PHONY: clean-xyplot-$(1)-$(2)


clean-xyplot-$(1): clean-xyplot-$(1)-$(2)

.PHONY: clean-xyplot-$(1)


clean-xyplot: clean-xyplot-$(1)

.PHONY: clean-xyplot


endef
