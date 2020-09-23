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


# functions to prepare daata for further processing and to plot the maps
# of the reduced data
#
include $(CFG_DIR)/prepdata.mk


# include function to include rules to look for clumps
#
include $(CFG_DIR)/findclumps.mk


# include function to calculate the physical parameters of the data
#
include $(CFG_DIR)/physcalcs.mk


# function to include the rules that plot maps and plots of the
# calculated physical parameters
#
include $(CFG_DIR)/physplots.mk



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
