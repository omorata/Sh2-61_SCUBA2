##
##  Makefile to run scripts to reduce the JCMT SCUBA2 data of Sh2-61
##
##  O. Morata
##  2018-20
##

##-- Info --------------------------------------------------------------
HOME_DIR=.
SNAME=Sh2_61

jointreductions = j850r0 j850r1 j450r0 j450r1
jointreductions += j850r0_mb j850r1_mb j450r0_mb j450r1_mb
jointreductions += j850r0_co j850r1_co
jointreductions += j850r0_co_mb

#
##-- End info ----------------------------------------------------------

# names of directories
#
BIN=$(HOME_DIR)/src
CFG_DIR=$(HOME_DIR)/config/reduction
DATA_DIR=$(HOME_DIR)/data
RES_DIR=$(HOME_DIR)/results

# defaults
#
SHELL := bash
.DELETE_ON_ERROR:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-varibles
MAKEFLAFS += --no-builtin_rules

export


##-- Template definition -----------------------------------------------


define DoReduction
# Template to process all days of the same reduction on the same step
#
#  Parameter: 1- target
#

$(eval reduc_file := $(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc.sdf)
$(eval snr_file := $(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc_snr.sdf)


$(reduc_file): $(wildcard $(CFG_DIR)/reduc-$(1).cfg)
	. $(BIN)/reduce_raw.sh $(CFG_DIR)/reduc-$(1).cfg


$(snr_file): $(reduc_file) 
	$(BIN)/post_scuba2.sh \
		-a snr  \
		-d $(RES_DIR)/$(1)  \
		-i $(SNAME)-$(1)-reduc


$(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc_crop.sdf: $(reduc_file)
	$(BIN)/post_scuba2.sh  \
		-a crop  \
		-d $(RES_DIR)/$(1) \
		-p $(CFG_DIR)/recipe-$(1).cfg  \
		-i $(SNAME)-$(1)-reduc

$(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc_snr_crop.sdf: $(snr_file) 
	$(BIN)/post_scuba2.sh  \
		-a crop  \
		-d $(RES_DIR)/$(1) \
		-p $(CFG_DIR)/recipe-$(1).cfg  \
		-i $(SNAME)-$(1)-reduc_snr


.PHONY: reduce-$(1) snr-$(1) crop-$(1) snrcrop-$(1)
reduce-$(1) : $(reduc_file)
snr-$(1): $(snr_file)
crop-$(1): $(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc_crop.sdf
snrcrop-$(1): $(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc_snr_crop.sdf

.PHONY: $(1)
$(1): reduce-$(1) snr-$(1) crop-$(1) snrcrop-$(1)

endef


##-- End of template definition ----------------------------------------


# define rules for reductions
#
$(foreach reduction, $(jointreductions),\
    $(eval $(call DoReduction,$(reduction)))\
)


