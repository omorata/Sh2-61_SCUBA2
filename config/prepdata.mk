#  prepdata.mk
#  O. Morata 2020
#
#  function to generate the rules to prepare reduced data for further
#  processing and to plot the resulting maps of the JCMT SCUBA2 data of
#  Sh2-61
#
#  To be included by the main Makefile



define PlotMaps
# Template to plot maps for targets
#
# Parameter: 1- target
#
$(eval map_dir := $(RES_DIR)/analysis_maps)

$(eval reduc_file := $(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc.sdf)
$(eval reduc_snrfile := $(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc_snr.sdf)

$(eval fits_reducfile := $(map_dir)/$(SNAME)-$(1)-reduc.fits)
$(eval fits_reducsnrfile := $(map_dir)/$(SNAME)-$(1)-reduc_snr.fits)

$(fits_reducfile): $(reduc_file)
	@ $(BIN_DIR)/prepare_maps.sh \
            -f $(reduc_file) \
            -o $(fits_reducfile) \
            -t "tofits" 

$(fits_reducsnrfile): $(reduc_snrfile)
	@ $(BIN_DIR)/prepare_maps.sh \
            -f $(reduc_snrfile) \
            -o $(fits_reducsnrfile) \
            -t "tofits"


tofits-$(1): $(fits_reducfile) $(fits_reducsnrfile)

.PHONY: tofits-$(1)


tofits: tofits-$(1)

.PHONY: tofits


clean-fits_datasets-$(1):
	@rm -fv $(fits_origfile)
	@rm -fv $(fits_origsnrfile)

.PHONY: clean-fits_datasets-$(1)


clean-fits_datasets: clean-fits_datasets-$(1)

.PHONY: clean-fits_datasets



$(eval orig_file := $(RES_DIR)/analysis_maps/$(SNAME)-$(1).fits)

$(eval map_file := $(map_dir)/$(SNAME)-$(1)-map.pdf)
$(eval cfg_file := $(CFG_DIR)/figures/$(SNAME)-$(1)-map.yml)

$(map_file): $(orig_file) $$(wildcard $$(cfg_file))
	@if [ -f $(cfg_file) ]; then \
	    $(EXT_DIR)/dbxmap.py \
                -c $(cfg_file) \
                -o $(map_dir) \
                -w $(map_dir) ;\
        fi
        #else \
            echo -e "\n++ Ignoring rule:\n      $(out_fc)" ;\
            echo -e "    No cfg file:\n      $(cfg_file)" ;\
        fi

plotmap-$(1): $(map_file)
.PHONY: plotmap-$(1)

plotmaps: plotmap-$(1)
.PHONY: plotmaps

clean-plotmap-$(1):
	@rm -fv $(map_file)
.PHONY: clean-plotmap-$(1)

clean-plotmaps: clean-plotmap-$(1)
.PHONY: clean-plotmaps

endef




define PrepareDataset
# Template to prepare data files for further processing
#
#  Parameter: 1- target
#
$(eval out_dir := $(RES_DIR)/analysis_maps)

$(eval tgt_dir := $(RES_DIR)/$(1))

$(eval orig_file := $(tgt_dir)/$(SNAME)-$(1)-reduc.sdf)
$(eval orig_snrfile := $(tgt_dir)/$(SNAME)-$(1)-reduc_snr.sdf)

$(eval strip_file := $(out_dir)/$(SNAME)-$(1).sdf)
$(eval strip_snrfile := $(out_dir)/$(SNAME)-$(1)-snr.sdf)

$(strip_file): $(orig_file) 
	@$(BIN_DIR)/prepare_maps.sh \
            -f $(orig_file) \
            -o $(strip_file) \
            -t "strip" 

$(strip_snrfile): $(orig_snrfile) 
	@$(BIN_DIR)/prepare_maps.sh \
            -f $(orig_snrfile) \
            -o $(strip_snrfile) \
            -t "strip" 

strip-$(1): $(strip_file) $(strip_snrfile)

.PHONY: strip-$(1)


strip: strip-$(1)

.PHONY: strip




$(eval fits_stripfile := $(out_dir)/$(SNAME)-$(1).fits)
$(eval fits_stripsnrfile := $(out_dir)/$(SNAME)-$(1)-snr.fits)

$(fits_stripfile): $(strip_file)
	@$(BIN_DIR)/prepare_maps.sh \
            -f $(strip_file) \
            -o $(fits_stripfile) \
            -t "tofits" 

$(fits_stripsnrfile): $(strip_snrfile)
	@$(BIN_DIR)/prepare_maps.sh \
            -f $(strip_snrfile) \
            -o $(fits_stripsnrfile) \
            -t "tofits" 

tofits_strip-$(1): $(fits_stripfile) $(fits_stripsnrfile)

.PHONY: tofits_strip-$(1)


tofits_strip: tofits_strip-$(1)

.PHONY: tofits_strip


clean-strip-$(1):
	@rm -fv $(strip_file)
	@rm -fv $(strip_snrfile)
	@rm -fv $(fits_stripfile)
	@rm -fv $(fits_stripsnrfile)

.PHONY: clean-strip-$(1)


clean-strip: clean-strip-$(1)

.PHONY: clean-strip

endef
