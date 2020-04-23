##
##  Makefile to run scripts to analyse the JCMT SCUBA2 data of Sh2-61
##
##  O. Morata
##  2020
##

##-- Info --------------------------------------------------------------
PRJ_NAME=Sh2_61-SCUBA2
#HOME_DIR := /lustre/opsw/work/omoratac/Sh2-61/SCUBA2
HOME_DIR := .
SNAME := Sh2_61

targets := j850r0_co j850r0_co_mb j850r0 j850r0_mb j850r1 j850r1_mb
targets += j450r0 j450r0_mb j450r1 j450r1_mb
fcs := fw_01 fw_02
combined := j850r0_co_mb__j450r0_mb

comb_maps := ratio tdust N mass
#
##-- End info ----------------------------------------------------------

# names of directories
#
BIN_DIR := $(HOME_DIR)/src
CFG_DIR := $(HOME_DIR)/config
DATA_DIR := $(HOME_DIR)/results
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


##-- Template definition -----------------------------------------------

define PlotMaps
# Template to plot maps for targets
#
# Parameter: 1- target
#
$(eval map_dir := $(RES_DIR)/analysis_maps)

$(eval tgt_dir := $(RES_DIR)/analysis_maps)

$(eval orig_file := $(tgt_dir)/$(SNAME)-$(1).fits)

$(eval map_file := $(map_dir)/$(SNAME)-$(1)-map.pdf)
$(eval cfg_file := $(CFG_DIR)/figures/$(SNAME)-$(1)-map.yml)

$(map_file): $(orig_file) $$(wildcard $$(cfg_file))
	@if [ -f $(cfg_file) ]; then \
	     $(EXT_DIR)/dbxmap.py \
                 -c $(cfg_file) \
                 -o $(map_dir) \
                 -w $(map_dir) ;\
         else \
             echo -e "\n++ Ignoring rule $(out_fc)" ;\
             echo -e "    No cfg file $(cfg_file)" ;\
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

$(eval tgt_dir := $(DATA_DIR)/$(1))

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


$(eval fits_origfile := $(out_dir)/$(SNAME)-$(1)-reduc.fits)
$(eval fits_origsnrfile := $(out_dir)/$(SNAME)-$(1)-reduc_snr.fits)

$(fits_origfile): $(orig_file)
	@ $(BIN_DIR)/prepare_maps.sh \
                 -f $(orig_file) \
                 -o $(fits_origfile) \
                 -t "tofits" 

$(fits_origsnrfile): $(orig_snrfile)
	@ $(BIN_DIR)/prepare_maps.sh \
                 -f $(orig_snrfile) \
                 -o $(fits_origsnrfile) \
                 -t "tofits"

tofits-$(1): $(fits_origfile) $(fits_origsnrfile)
.PHONY: tofits-$(1)

tofits: tofits-$(1)
.PHONY: tofits




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


clean-fits_datasets-$(1):
	@rm -fv $(fits_origfile)
	@rm -fv $(fits_origsnrfile)
	@rm -fv $(fits_stripfile)
	@rm -fv $(fits_stripsnrfile)
.PHONY: clean-fits_datasets-$(1)

clean-fits_datasets: clean-fits_datasets-$(1)
.PHONY: clean-fits_datasets

endef




define Findclumps
# Template to find clumps in continuum maps
#
#  Parameters: 1- target; 2- findclump tag
#
$(eval analysis_dir := $(DATA_DIR)/analysis_maps)
$(eval findclumps_dir := $(RES_DIR)/findclumps)


$(eval cfg_file := $(CFG_DIR)/analysis/$(SNAME)-$(1)-$(2).cfg)

$(eval in_fc := $(analysis_dir)/$(SNAME)-$(1).sdf)
$(eval insnr_fc := $(analysis_dir)/$(SNAME)-$(1)-snr.sdf)

$(eval out_fc := $(findclumps_dir)/$(SNAME)-$(1)-$(2)-clumps.sdf)
$(eval out_fc_fits := $(findclumps_dir)/$(SNAME)-$(1)-$(2)-clumps.fits)


$(out_fc): $$(wildcard $$(cfg_file)) $(in_fc) $(insnr_fc)
	@if [ -f $(cfg_file) ]; then \
	     sh $(BIN_DIR)/findclumps.sh \
                 -c $(cfg_file) \
                 -o $(findclumps_dir) \
                 -i $(analysis_dir) \
                 -d $(CFG_DIR)/analysis; \
         else \
             echo -e "\n++ Ignoring rule $(out_fc)" ;\
             echo -e "    No cfg file $(cfg_file)" ;\
         fi


$(out_fc_fits): $$(wildcard $$(cfg_file)) $(out_fc)
	@if [ -f $(cfg_file) ]; then \
             $(BIN_DIR)/prepare_maps.sh \
                 -f $(out_fc) \
                 -o $(out_fc_fits) \
                 -t "tofits" ;\
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
$(out_shapes): $$(wildcard $$(out_catalog))
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

$(cl_mapfile): $(map_file) $$(wildcard $$(cfg_mapfile))
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




define Align_Dataset
# Template to align one dataset to the grid of another one
# 
# 1- combined string

$(eval ref_tgt := $(firstword $(subst __, ,$(1))))
$(eval sec_tgt := $(lastword $(subst __, ,$(1))))

$(eval outdir := $(DATA_DIR)/analysis_maps)


$(eval align_ref := $(outdir)/$(SNAME)-$(ref_tgt).sdf)

$(eval align_tgt := $(outdir)/$(SNAME)-$(sec_tgt).sdf)
$(eval align_snrtgt := $(outdir)/$(SNAME)-$(sec_tgt)-snr.sdf)

$(eval aligned_file := $(outdir)/$(SNAME)-$(sec_tgt)-aligned_to-$(ref_tgt).sdf)
$(eval aligned_snrfile :=    \
    $(outdir)/$(SNAME)-$(sec_tgt)-snr-aligned_to-$(ref_tgt).sdf)

$(eval aligned_fits := $(outdir)/$(SNAME)-$(sec_tgt)-aligned_to-$(ref_tgt).fits)
$(eval aligned_snrfits :=   \
    $(outdir)/$(SNAME)-$(sec_tgt)-snr-aligned_to-$(ref_tgt).fits)



$(aligned_file): $(align_tgt) $(align_ref)
	$(BIN_DIR)/prepare_maps.sh \
             -f $(align_tgt) \
             -o $(aligned_file) \
             -r $(align_ref) \
             -t "align" 

$(aligned_snrfile): $(align_snrtgt) $(align_ref)
	@$(BIN_DIR)/prepare_maps.sh \
             -f $(align_snrtgt) \
             -o $(aligned_snrfile) \
             -r $(align_ref) \
             -t "align" 


$(aligned_fits): $(aligned_file) 
	@ $(BIN_DIR)/prepare_maps.sh \
                 -f $(aligned_file) \
                 -o $(aligned_fits) \
                 -t "tofits" 

$(aligned_snrfits): $(aligned_snrfile) 
	@ $(BIN_DIR)/prepare_maps.sh \
                 -f $(aligned_snrfile) \
                 -o $(aligned_snrfits) \
                 -t "tofits" 

.INTERMEDIATE: $(aligned_file) $(aligned_snrfile)

align-$(sec_tgt)-to-$(ref_tgt): $(aligned_fits) $(aligned_snrfits)
.PHONY: align-$(sec_tgt)-to-$(ref_tgt)

.PHONY: clean-align-$(1)
clean-align-$(1):
	@rm -fv $(aligned_fits)
	@rm -fv $(aligned_snrfits)

.PHONY: clean-align
clean-align: clean-align-$(1)

.PHONY: clean-align-$(sec_tgt)
clean-align-$(sec_tgt): clean-align-$(1)


endef




define CalcPhysParam
# Template to calculate the physical parameters from the datasets
#
# 1- combined string, 2- findclump_id

$(eval ref_tgt := $(firstword $(subst __, ,$(1))))
$(eval sec_tgt := $(lastword $(subst __, ,$(1))))

$(eval outdir := $(DATA_DIR)/analysis_maps)


$(eval aligned_fits := $(outdir)/$(SNAME)-$(sec_tgt)-aligned_to-$(ref_tgt).fits)
$(eval aligned_snrfits :=    \
    $(outdir)/$(SNAME)-$(sec_tgt)-snr-aligned_to-$(ref_tgt).fits)


#
# Calculation of physical parameters
#

$(eval ffile := $(outdir)/$(SNAME)-$(ref_tgt).fits)
$(eval ffile_snr := $(outdir)/$(SNAME)-$(ref_tgt)-snr.fits)
$(eval clfile := $(DATA_DIR)/findclumps/$(SNAME)-$(ref_tgt)-$(2)-clumps.fits)

$(eval cfg_file := $(CFG_DIR)/analysis/$(SNAME)-$(1)-$(2)-phys_calc.yaml)

$(eval calc_log := $(outdir)/calcs-$(1)-$(2).log)

$(eval calc_refs :=    \
    $(ffile) $(ffile_snr) $(aligned_fits) $(aligned_snrfits) $(clfile))

$(calc_log):  $$(wildcard $$(cfg_file)) $(calc_refs)
	@if [ -f $(cfg_file) ]; then \
	     $(BIN_DIR)/calc_phys_parm.py \
                 -c $(cfg_file) \
                 --clumps \
                 -l $(calc_log) \
                 -w $(RES_DIR);\
          else \
             echo -e "\n++ Ignoring rule $(out_fc)" ;\
             echo -e "    No cfg file $(cfg_file)" ;\
         fi

.PHONY: calcs-$(1)-$(2)
calcs-$(1)-$(2) : $(calc_log)

.PHONY: calcs-$(1)
calcs-$(1) : calcs-$(1)-$(2)

.PHONY: calcs-$(ref_tgt)
calcs-$(ref_tgt) : calcs-$(1)


.PHONY: clean-calcs-$(1)-$(2)
clean-calcs-$(1)-$(2):
	@rm -fv $(calc_log)
	@rm -fv $(outdir)/$(SNAME)-$(1)-$(2)-ratio.fits
	@rm -fv $(outdir)/$(SNAME)-$(1)-$(2)-tdust.fits
	@rm -fv $(outdir)/$(SNAME)-$(1)-$(2)-mass.fits
	@rm -fv $(outdir)/$(SNAME)-$(1)-$(2)-clump_table.fits

.PHONY: clean-calcs-$(1)
clean-calcs-$(1): clean-calcs-$(1)-$(2)

.PHONY: clean-calcs
clean-calcs: clean-calcs-$(1)

.PHONY: clean-calcs-$(ref_tgt)
clean-calcs-$(ref_tgt): clean-calcs-$(1)

endef




define MapPhysParam
# Template to plot maps of the calculated physical parameter
#
# Arguments: 1- tgt, 2- findclump id, 3- parameter
#
$(eval tgt_dir := $(RES_DIR)/analysis_maps)
$(eval out_dir := $(RES_DIR)/analysis_maps)

$(eval orig_file := $(tgt_dir)/$(SNAME)-$(1)-$(2)-$(3).fits)

$(eval out_file := $(out_dir)/$(SNAME)-$(1)-$(2)-$(3)-map.pdf)
$(eval cfg_file := $(CFG_DIR)/analysis/$(SNAME)-$(1)-$(2)-$(3)-map.yml)

$(out_file): $(orig_file) $$(wildcard $$(cfg_file))
	@if [ -f $(cfg_file) ]; then \
	     $(EXT_DIR)/dbxmap.py \
                 -c $(cfg_file) \
                 -o $(out_dir) \
                 -w $(tgt_dir) ;\
         else \
             echo -e "\n++ Ignoring rule $(out_file)" ;\
             echo -e "    No cfg file $(cfg_file)" ;\
         fi


map-physpar-$(1)-$(2)-$(3): $(out_file)
.PHONY: map-physpar-$(1)-$(2)-$(3)

maps-physpar-$(1)-$(2): map-physpar-$(1)-$(2)-$(3)
.PHONY: maps-physpar-$(1)-$(2)

maps-physpar-$(1): map-physpar-$(1)-$(2)
.PHONY: maps-physpar-$(1)

maps-physpar: map-physpar-$(1)
.PHONY: maps-physpar


clean-map-physpar-$(1)-$(2)-$(3):
	@rm -fv $(out_file)
.PHONY: clean-map-physpar-$(1)-$(2)-$(3)


clean-maps-physpar-$(1)-$(2): clean-map-physpar-$(1)-$(2)-$(3)
.PHONY: clean-maps-physpar-$(1)-$(2)

clean-maps-physpar-$(1): clean-map-physpar-$(1)-$(2)
.PHONY: clean-maps-physpar-$(1)

clean-maps-physpar: clean-map-physpar-$(1)
.PHONY: clean-maps-physpar-$(1)

endef



##-- End of template definition ----------------------------------------


# define rules for targets
#
$(foreach tgt, $(targets),\
    $(eval $(call PlotMaps,$(tgt)))\
    $(eval $(call PrepareDataset,$(tgt)))\
)

# define rules for findclumps
#
$(foreach tgt, $(targets),\
    $(foreach fc, $(fcs),\
        $(eval $(call Findclumps,$(tgt),$(fc)))\
    ) \
)

# define rules for combined and findclumps_id
#
$(foreach tgt, $(combined),\
    $(eval $(call Align_Dataset,$(tgt)))\
    $(foreach fc, $(fcs),\
        $(eval $(call CalcPhysParam,$(tgt),$(fc)))\
        $(foreach mp, $(comb_maps),\
           $(eval $(call MapPhysParam,$(tgt),$(fc),$(mp)))\
        )\
    ) \
)



# other rules
#
clean_list := clean-strip clean-align clean-findclumps clean-calcs
clean_list += clean-plotmaps clean-maps_physpar clean-clumps-map

clean:	$(clean_list)
.PHONY: clean



help:
	@echo;echo  "Makefile to analyse the data of $(PRJ_NAME)"
	@echo "-------------------------------------------"
	@echo "  pre-defined variables:"
	@echo "             Project Name : $(PRJ_NAME)"
	@echo "              File prefix : $(SNAME)"
	@echo "                  targets : $(targets)"
	@echo "           findclumps ids : $(fcs)"
	@echo "         combined targets : $(combined)"
	@echo " physical parameters maps : $(comb_maps)"
	@echo;echo;echo "  More help options:"
	@echo "      make help_dirs  -  information on directories"
	@echo "      make help_rules -  information on defined rules"
	@echo



help_dirs:
	@echo
	@echo " Directory set-up for project $(PRJ_NAME):"
	@echo "---------------------------------------------------------"
	@echo "   project home : $(HOME_DIR)"
	@echo "           bin  : $(BIN_DIR)"
	@echo "  configuration : $(CFG_DIR)"
	@echo "           data : $(DATA_DIR)"
	@echo "        results : $(RES_DIR)"
	@echo "   external bin : $(EXT_DIR)"
	@echo



help_rules:
	@echo;echo "-----------------"
	@echo "  Defined rules"
	@echo "-----------------"
	@echo;echo "  (Warning: not all the following rules may be available."
	@echo "   It will depend on the definition of the corresponding"
	@echo "   configuration files)"
	@echo;echo " The general actions are:"
	@echo;echo "    make plotmaps  --  plot maps of targets"
	@echo "    make strip  --  strip third axis from sdf files of targets"
	@echo "    make tofits  --  transform .sdf files of targets to .fits"
	@echo "    make tofits-strip  --  transform stripped files to .fits "
	@echo "    [ make align ] -- align map to reference"
	@echo "    make findclumps_snr  --  find clumps in map using snr map"
	@echo "    make polygonfiles  --  extract shapes of clumps"
	@echo "    make clumps-map  --  plot overlay of clumps on emission map"
	@echo "    [ make calcs ] -- calculate physical parameters from emission and clumps"
	@echo "    make maps-physpar  -- plot map of physical parameters"
	@echo "   clean options: clean-plotmaps clean-strip clean-fits_datasets clean-findclumps"
	@echo "        clean-clumps-map clean-align clean-maps-physpar clean-calcs"
	@echo;echo " In more detail:"
	@echo;echo " + rules depending on the target"
	@echo "    (targets: $(targets))"
	@echo "     make plotmap-[target]"
	@echo "     make strip-[target]"
	@echo "     make tofits-[target]"
	@echo "     make tofits_strip-[target]"
	@echo "     make findclumps_snr-[target]"
	@echo "     make polygonfiles-[target]"
	@echo "     make clumps-map-[target]"
	@echo "    clean options: clean-plotmap-[target] clean-strip-[target]"
	@echo "         clean-fits_datasets-[target] clean-findclumps-[target]"
	@echo "         clean-clumps-map-[target] clean-calcs"
	@echo;echo " + rules depending on the target and findclumps_id"
	@echo "    (ids: $(fcs))"
	@echo "     make findclumps_snr-[target]-[id]"
	@echo "     make polygonfiles-[target]-[id]"
	@echo "     make clumps-map-[target]-[id]"
	@echo "    clean options: clean-findclumps-[target]-[id] clean-clumps-map-[target]-[id]"
	@echo;echo " + rules depending on the combination map"
	@echo "    (comb. maps: $(combined))"
	@echo "     make align-[secondary_target]-to-[reference_target]"
	@echo "     make calcs-[comb_map]"
	@echo "     make calcs-[reference_target]"
	@echo "     make maps-physpar-[comb.map]"
	@echo "    clean options: clean-align-[comb_map] clean-align-[secondary-target]"
	@echo "         clean-calcs-[comb_map] clean-calcs-[reference_target]"
	@echo "         clean-maps-physpar-[comb.map]"
	@echo;echo " + rules depending on combination map and findclump_id"
	@echo "     make calcs-[combmap]-[id]"
	@echo "     make maps-physpar-[combmap]-[id]"
	@echo "    clean options: clean-calcs-[combmap]-[id] clean-maps-physpar-[combmap]-[id]"
	@echo;echo " + rules depending on the comb. map, findclumps_id, and physical parameter"
	@echo "    (phys. param: $(comb_maps))"
	@echo "     make map-physpar-[comb.map]-[id]-[physpar]"
	@echo "    clean options:  clean-map-physpar[comb.map]-[id]-[physpar]"
	@echo

.PHONY: help help_dirs help_rules




.SHELLFLAGS = -ec
list:
# lists all the rules in the Makefile
#
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2> /dev/null |\
           awk -v RS= -F: '/^# File/,/^# Finished Make data base/  \
             {if ($$1 !~ "^[#.]") {print $$1}}' | sort | \
           egrep -v -e '^[^[:alnum:]]' -e '^$@$$' 
##


list_deps:

	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null |\
           awk -v RS= -F: '/^# File/,/^# Finished Make data base/  \
             {if ($$1 !~ "^[#.]") {print $$1" ==> "$$2}}' | sort | \
           egrep -v -e '^[^[:alnum:]]' -e '^$@$$'
##


list_files:
# lists all the files created in the rules in the Makefile
#
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | \
           awk -v RS= -F: '/^# File/,/^# Finished Make data base/ \
               {if ($$1 !~ "^[#.]") {print $$1}}' | sort | \
           egrep -e '\/'

.PHONY: list list_deps list_files

##
##-- End of rules ------------------------------------------------------
