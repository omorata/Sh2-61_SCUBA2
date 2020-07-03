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
BIN_DIR=$(HOME_DIR)/src
CFG_DIR=$(HOME_DIR)/config
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
# Template to process all observations in a single step
#
#  Parameter: 1- target
#

$(eval reduc_file := $(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc.sdf)
$(eval snr_file := $(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc_snr.sdf)


$(reduc_file): $(wildcard $(CFG_DIR)/reduction/reduc-$(1).cfg)
	. $(BIN_DIR)/reduce_raw.sh $(CFG_DIR)/reduction/reduc-$(1).cfg


$(snr_file): $(reduc_file) 
	$(BIN_DIR)/post_scuba2.sh \
		-a snr  \
		-d $(RES_DIR)/$(1)  \
		-i $(SNAME)-$(1)-reduc


$(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc_crop.sdf: $(reduc_file)
	$(BIN_DIR)/post_scuba2.sh  \
		-a crop  \
		-d $(RES_DIR)/$(1) \
		-p $(CFG_DIR)/reduction/recipe-$(1).cfg  \
		-i $(SNAME)-$(1)-reduc

$(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc_snr_crop.sdf: $(snr_file) 
	$(BIN_DIR)/post_scuba2.sh  \
		-a crop  \
		-d $(RES_DIR)/$(1) \
		-p $(CFG_DIR)/reduction/recipe-$(1).cfg  \
		-i $(SNAME)-$(1)-reduc_snr


map-$(1) : $(reduc_file)
.PHONY: map-$(1)

snr-$(1): $(snr_file)
.PHONY: snr-$(1)

crop-$(1): $(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc_crop.sdf
.PHONY: crop-S(1)

snrcrop-$(1): $(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc_snr_crop.sdf
.PHONY: snrcrop-$(1)

reduce-$(1): map-$(1) snr-$(1) crop-$(1) snrcrop-$(1)
.PHONY: reduce-$(1)


clean-map-$(1):
	@rm -vf $(reduc_file)
.PHONY: clean-map-$(1)

clean-snr-$(1):
	@rm -vf $(snr_file)
.PHONY: clean-snr-$(1)

clean-crop-$(1):
	@rm -vf $(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc_crop.sdf
.PHONY: clean-crop-$(1)

clean-snrcrop-$(1):
	@rm -vf $(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc_snr_crop.sdf
.PHONY: clean-snrcrop-$(1)

clean-reduce-$(1): clean-map-$(1) clean-snr-$(1) clean-crop-$(1)
clean-reduce-$(1): clean-snrcrop-$(1)

.PHONY: clean-reduce-$(1)

clean: clean-reduce-$(1)
.PHONY: clean

endef


##-- End of template definition ----------------------------------------


# define rules for reductions
#
$(foreach reduction, $(jointreductions),\
    $(eval $(call DoReduction,$(reduction)))\
)


