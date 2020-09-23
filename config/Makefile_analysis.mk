##
##  Makefile to run scripts to reduce and analyse the JCMT SCUBA2 data
##  of Sh2-61
##
##  O. Morata
##  2020
##

##-- Info --------------------------------------------------------------
PRJ_NAME=Sh2_61-SCUBA2
HOME_DIR := .
SNAME := Sh2_61

targets := j850r0 j850r1 j450r0 j450r1
targets += j850r0_mb j850r1_mb j450r0_mb j450r1_mb
targets += j850r0_co j850r1_co
targets += j850r0_co_mb

findclump_tags := fw_01 fw_02 cf_01

combined := j850r0_co_mb__j450r0_mb j850r0_mb__j450r0_mb
combined += j850r1_mb__j450r0_mb

comb_maps := ratio tdust N mass

physcalc_tags := ref Tcalc Tfix Tffx
physcalc_tags += td1 td2 td3 beta1 beta2 beta3

histo_tags := N tdust N_superp tdust_superp

xyplots_tags := N_vs_tdust N_vs_f850
#
##-- End info ----------------------------------------------------------

# names of directories
#
BIN_DIR := $(HOME_DIR)/src
CFG_DIR := $(HOME_DIR)/config
DATA_DIR := $(HOME_DIR)/data
RES_DIR := $(HOME_DIR)/results
EXT_DIR := $(HOME_DIR)/bin

# defaults
#
SHELL := /bin/bash
.DELETE_ON_ERROR:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-varibles
MAKEFLAFS += --no-builtin_rules


export


# include function to define rules for reduction
#
include $(CFG_DIR)/reduction.mk

#include $(CFG_DIR)/????.mk
#include $(CFG_DIR)/findclumps.mk


##-- Template definition -----------------------------------------------

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
         else \
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
	@ $(BIN_DIR)/prepare_maps.sh \
                 -f $(orig_file) \
                 -o $(strip_file) \
                 -t "strip" 

$(strip_snrfile): $(orig_snrfile) 
	@$(BIN_DIR)/prepare_maps.sh \
                 -f $(orig_snrfile) \
                 -o $(strip_snrfile) \
                 -t "strip" 

.PHONY: strip-$(1)
strip-$(1): $(strip_file) $(strip_snrfile)

.PHONY: strip
strip: strip-$(1)




$(eval fits_stripfile := $(out_dir)/$(SNAME)-$(1).fits)
$(eval fits_stripsnrfile := $(out_dir)/$(SNAME)-$(1)-snr.fits)

$(fits_stripfile): $(strip_file)
	@ $(BIN_DIR)/prepare_maps.sh \
            -f $(strip_file) \
            -o $(fits_stripfile) \
            -t "tofits" 

$(fits_stripsnrfile): $(strip_snrfile)
	@ $(BIN_DIR)/prepare_maps.sh \
            -f $(strip_snrfile) \
            -o $(fits_stripsnrfile) \
            -t "tofits" 

.PHONY: tofits_strip-$(1)
tofits_strip-$(1): $(fits_stripfile) $(fits_stripsnrfile)

.PHONY: tofits_strip
tofits_strip: tofits_strip-$(1)


.PHONY: clean-strip-$(1)
clean-strip-$(1):
	@rm -fv $(strip_file)
	@rm -fv $(strip_snrfile)
	@rm -fv $(fits_stripfile)
	@rm -fv $(fits_stripsnrfile)

.PHONY: clean-strip
clean-strip: clean-strip-$(1)



endef




define Findclumps
# Template to find clumps in continuum maps
#
#  Parameters: 1- target; 2- findclump tag
#
$(eval analysis_dir := $(RES_DIR)/analysis_maps)
$(eval findclumps_dir := $(RES_DIR)/findclumps)


$(eval cfg_file := $(CFG_DIR)/analysis/$(SNAME)-$(1)-$(2).cfg)
$(eval par_file := $(CFG_DIR)/analysis/$(SNAME)-$(1)-$(2).par)

$(eval in_fc := $(analysis_dir)/$(SNAME)-$(1).sdf)
$(eval insnr_fc := $(analysis_dir)/$(SNAME)-$(1)-snr.sdf)

$(eval out_fc := $(findclumps_dir)/$(SNAME)-$(1)-$(2)-clumps.sdf)
$(eval out_fc_fits := $(findclumps_dir)/$(SNAME)-$(1)-$(2)-clumps.fits)


$(out_fc): $$(wildcard $$(cfg_file) $$(par_file)) $(in_fc) $(insnr_fc)
	@if [ -f $(cfg_file) ]; then \
	     sh $(BIN_DIR)/findclumps.sh \
                 -c $(cfg_file) \
                 -o $(findclumps_dir) \
                 -i $(analysis_dir) \
                 -d $(CFG_DIR)/analysis; \
         else \
             echo -e "\n++ Ignoring rule $(out_fc)" ;\
             echo -e "    No cfg file $(cfg_file)" ;\
             touch $(out_fc);\
         fi


$(out_fc_fits): $$(wildcard $$(cfg_file)) $(out_fc)
	@if [ -f $(cfg_file) ]; then \
             $(BIN_DIR)/prepare_maps.sh \
                 -f $(out_fc) \
                 -o $(out_fc_fits) \
                 -t "tofits" ;\
         else \
             echo -e "\n++ Ignoring rule $(out_fc)" ;\
             echo -e "    No cfg file $(cfg_file)" ;\
             touch $(out_fc_fits);\
         fi


.PHONY: findclumps_snr-$(1)-$(2)
findclumps_snr-$(1)-$(2): $(out_fc) $(out_fc_fits)

.PHONY: findclumps_snr-$(1)
findclumps_snr-$(1): findclumps_snr-$(1)-$(2)

.PHONY: findclumps_snr
findclumps_snr: findclumps_snr-$(1)


.PHONY: clean-findclumps-$(1)-$(2)
clean-findclumps-$(1)-$(2):
	@rm -fv $(out_fc)
	@rm -fv $(out_fc_fits)
	@rm -fv $(findclumps_dir)/$(SNAME)-$(1)-$(2)-catalog.fits
	@rm -fv $(findclumps_dir)/$(SNAME)-$(1)-$(2).log
	@rm -fv $(findclumps_dir)/$(SNAME)-$(1)-snr-$(2)-clumps.sdf
	@rm -fv $(findclumps_dir)/$(SNAME)-$(1)-snr-$(2)-catalog.fits
	@rm -fv $(findclumps_dir)/$(SNAME)-$(1)-snr-$(2).log

.PHONY: clean-findclumps-$(1)
clean-findclumps-$(1): clean-findclumps-$(1)-$(2)

.PHONY: clean-findclumps
clean-findclumps: clean-findclumps-$(1)


$(eval out_shapes := $(findclumps_dir)/$(SNAME)-$(1)-$(2)-shapes.dat)
$(eval out_catalog := $(findclumps_dir)/$(SNAME)-$(1)-$(2)-catalog.fits)
$(out_shapes): $$(wildcard $$(out_catalog) $$(out_fc))
	@if [ -f $(out_catalog) ];then \
	    sh $(BIN_DIR)/catalog_to_polygonfile.sh \
                $(out_catalog) $(out_shapes);\
        fi


polygonfile-$(1)-$(2): $(out_shapes)

.PHONY: polygonfile-$(1)-$(2)


polygonfiles-$(1): polygonfile-$(1)-$(2)

.PHONY: polygonfiles-$(1)


polygonfiles: polygonfiles-$(1)

.PHONY: polygonfiles


clean-polygonfile-$(1)-$(2):
	@rm -fv $(out_shapes)

.PHONY: clean-polygonfile-$(1)-$(2)


clean-polygonfiles-$(1): clean-polygonfile-$(1)-$(2)

.PHONY: clean-polygonfiles-$(1)


clean-polygonfiles: clean-polygonfiles-$(1)

.PHONY: polygonfiles



$(eval map_file := $(analysis_dir)/$(SNAME)-$(1).fits)
$(eval cfg_mapfile := $(CFG_DIR)/figures/$(SNAME)-$(1)-$(2)-clumps-map.yml)
$(eval cl_mapfile := $(analysis_dir)/$(SNAME)-$(1)-$(2)-clumps-map.pdf)

$(cl_mapfile): $(map_file) $$(wildcard $$(cfg_mapfile) $$(out_shapes))
	@if [ -f $(cfg_mapfile) ]; then \
	     $(EXT_DIR)/dbxmap.py \
                 -c $(cfg_mapfile) \
                 -o $(out_dir) \
                 -w $(out_dir) ;\
         else \
             echo -e "\n++ Ignoring rule $(out_file)" ;\
             echo -e "    No cfg file $(cfg_mapfile)" ;\
         fi

clumps-map-$(1)-$(2): $(cl_mapfile)
.PHONY: clumps-map-$(1)-$(2)

clumps-map-$(1): clumps-map-$(1)-$(2)
.PHONY: clumps-map-$(1)

clumps-map: clumps-map-$(1)
.PHONY: clumps-map


clean-clumps-map-$(1)-$(2):
	@rm -vf $(cl_mapfile)
.PHONY: clean-clumps-map-$(1)-$(2)

clean-clumps-map-$(1): clean-clumps-map-$(1)-$(2)
.PHONY: clean-clumps-map-$(1)

clean-clumps-map: clean-clumps-map-$(1)
.PHONY: clean-clumps-map

endef



include $(CFG_DIR)/physcalcs.mk


include $(CFG_DIR)/physplots.mk



##-- End of template definition ----------------------------------------


list_tasks := reduction
list_tasks += tofits plotmaps
list_tasks += findclumps_snr polygonfiles clumps-map
list_tasks += strip tofits_strip
list_tasks += maps_physpar histo xyplot

#all: $(list_tasks)
#.PHONY: all


# define rules for reductions
#
$(foreach reduction, $(targets),\
    $(eval $(call DoReduction,$(reduction)))\
)


# define rules for targets
#
$(foreach tgt, $(targets),\
    $(eval $(call PlotMaps,$(tgt)))\
    $(eval $(call PrepareDataset,$(tgt)))\
)


# define rules for findclumps
#
$(foreach tgt, $(targets),\
    $(foreach fc, $(findclump_tags),\
        $(eval $(call Findclumps,$(tgt),$(fc)))\
    ) \
)


# define rules for combined and findclumps_id
#
$(foreach tgt, $(combined),\
    $(eval $(call Align_Dataset,$(tgt)))\
    $(foreach fc, $(findclump_tags),\
        $(foreach phyvr, $(physcalc_tags),\
            $(eval $(call CalcPhysParam,$(tgt),$(fc),$(phyvr)))\
            $(foreach mp, $(comb_maps),\
                $(eval $(call MapPhysParam,$(tgt),$(fc),$(phyvr),$(mp)))\
            )\
            $(foreach ht, $(histo_tags),\
                $(eval $(call HistoPlots,$(tgt),$(fc),$(phyvr),$(ht)))\
            )\
            $(foreach xy, $(xyplots_tags),\
                $(eval $(call XYPlots,$(tgt),$(fc),$(phyvr),$(xy)))\
            )\
        )\
    ) \
)



# other rules
#
clean_list := clean-reduction clean-fits_datasets clean-plotmaps clean-strip
clean_list += clean-findclumps clean-polygonfiles clean-clumps-map
clean_list += clean-align clean-calcs clean-physcatg
clean_list += clean-maps-physpar clean-histo clean-xyplot

clean:	$(clean_list)
.PHONY: clean


#clean_all
#clean_analysis  <=> clean
#clean_reduction (deactivate)
#clean_outputs


# include help rules
#
include $(CFG_DIR)/help.mk



# include list-related rules
#
include $(CFG_DIR)/lists.mk



##
##-- End of rules ------------------------------------------------------
