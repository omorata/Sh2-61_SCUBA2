##
##  Makefile to run scripts to reduce the JCMT SCUBA2 data of Sh2-61
##
##  O. Morata
##  2018-19
##

##-- Info --------------------------------------------------------------
HOME_DIR=/lustre/opsw/work/omoratac/Sh2-61/SCUBA2
SNAME=Sh2_61

#dayreductions = 850_r0
#days = 0424 0425

jointreductions = j850_r0 j850_r1 j450_r0 j450_r1
jointreductions += j850_r0mf j850_r1mf j450_r0mf j450_r1mf
jointreductions += j850_r0_contamination j850_r1_contamination
jointreductions += j850_r0_contamination_mf

#
##-- End info ----------------------------------------------------------

# names of directories
#
BIN=$(HOME_DIR)/scripts
CFG_DIR=$(HOME_DIR)/config/reduction
DATA_DIR=$(HOME_DIR)/data
RES_DIR=$(HOME_DIR)/results


export


##-- Template definition -----------------------------------------------

# Template to define the rules to process steps of each different day for
# the given reduction
#
define Days_Template

.PHONY: reduce-$(1)-$(2) cal-$(1)-$(2) checkcal-$(1)-$(2)


reduce-$(1)-$(2) : $(RES_DIR)/$(1)/$(SNAME)-$(1)-$(2)-reduc.sdf

cal-$(1)-$(2): $(RES_DIR)/$(1)/$(SNAME)-$(1)-$(2)-reduc_cal.sdf


$(RES_DIR)/$(1)/$(SNAME)-$(1)-$(2)-reduc.sdf: $(CFG_DIR)/reduction-$(1)-$(2).cfg
	. $(BIN)/reduce_raw.sh $(CFG_DIR)/reduction-$(1)-$(2).cfg


$(RES_DIR)/$(1)/$(SNAME)-$(1)-$(2)-reduc_cal.sdf: $(RES_DIR)/$(1)/$(SNAME)-$(1)-$(2)-reduc.sdf
	$(BIN)/post_scuba2.sh  \
		-a cal  \
		-d $(RES_DIR)/$(1)/ \
		-p $(CFG_DIR)/par$(1)-$(2).cfg  \
		-i $(SNAME)-$(1)-$(2)-reduc


checkcal-$(1)-$(2):
	. $(BIN)/check_fcf.sh \
		-d $(RES_DIR)/$(1)  \
		-c $(CFG_DIR)/reduction-$(1)-$(2).cfg

endef



# Template to run the rest of the steps of the defined reduction
#
define Tag_Template

.PHONY: reduce-$(1)-days cal-$(1)-days checkcal-$(1)
.PHONY: coadd-$(1) snr-$(1) crop-$(1)


reduce-$(1)-days: $(reduce_$(1)_list)

cal-$(1)-days: $(cal_$(1)_list)

checkcal-$(1): $(checkcal_$(1)_list)

coadd-$(1): $(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd.sdf

snr-$(1): $(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd_snr.sdf

crop-$(1): $(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd_crop.sdf

$(eval mosaic_dir := $(CFG_DIR)/mosaic_lists)


$(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd.sdf:  $(jfiles_$(1)_list) $(mosaic_dir)/$(1)-mosaic.list
	$(BIN)/post_scuba2.sh  \
		-a add  \
		-d $(RES_DIR)/$(1) \
		-p $(CFG_DIR)/par$(1).cfg \
		-l $(mosaic_dir)/$(1)-mosaic.list  \
		-o $(SNAME)-$(1)-coadd


$(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd_snr.sdf: $(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd.sdf
	$(BIN)/post_scuba2.sh \
		-a snr  \
		-d $(RES_DIR)/$(1)  \
		-i $(SNAME)-$(1)-coadd


$(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd_crop.sdf: $(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd.sdf
	$(BIN)/post_scuba2.sh  \
		-a crop  \
		-d $(RES_DIR)/$(1) \
		-p $(CFG_DIR)/par$(1).cfg  \
		-i $(SNAME)-$(1)-coadd

endef



# Template to process all days of the same reduction on the same step
#
define Joint_Template

.PHONY: $(1) reduce-$(1) snr-$(1) crop-$(1)

$(eval reduc_file := $(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc.sdf)


$(1): reduce-$(1) snr-$(1) crop-$(1)


reduce-$(1) : $(reduc_file)

snr-$(1): $(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc_snr.sdf

crop-$(1): $(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc_crop.sdf


$(reduc_file): $(CFG_DIR)/reduction-$(1).cfg
	. $(BIN)/reduce_raw.sh $(CFG_DIR)/reduction-$(1).cfg


$(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc_snr.sdf: $(reduc_file) 
	$(BIN)/post_scuba2.sh \
		-a snr  \
		-d $(RES_DIR)/$(1)  \
		-i $(SNAME)-$(1)-reduc


$(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc_crop.sdf: $(reduc_file)
	$(BIN)/post_scuba2.sh  \
		-a crop  \
		-d $(RES_DIR)/$(1) \
		-p $(CFG_DIR)/par$(1).cfg  \
		-i $(SNAME)-$(1)-reduc

endef


##-- End of template definition ----------------------------------------

##-- Generate lists ----------------------------------------------------

# loop to generate the lists needed for the rules that join data from more
# than one day
#
$(foreach reduction, $(dayreductions),\
    $(eval reduce_$(reduction)_list = ) \
    $(eval cal_$(reduction)_list = ) \
    $(eval checkcal_$(reduction)_list = ) \
    $(eval jfiles_$(reduction)_list = ) \
    $(foreach day, $(days),\
        $(eval reduce_$(reduction)_list += reduce-$(reduction)-$(day)) \
        $(eval cal_$(reduction)_list += cal-$(reduction)-$(day)) \
        $(eval checkcal_$(reduction)_list += checkcal-$(reduction)-$(day)) \
        $(eval jfiles_$(reduction)_list += $(RES_DIR)/$(reduction)/$(SNAME)-$(reduction)-$(day)-reduc_cal.sdf) \
))

##-- Enf of lists ------------------------------------------------------


# define rules for joint reductions
#
$(foreach reduction, $(jointreductions),\
    $(eval $(call Joint_Template,$(reduction)))\
)


#850_r0: reduce-850_r0-days cal-850_r0-days coadd-850_r0 snr-850_r0 crop-850_r0
#850_r1: reduce-850_r1-days cal-850_r1-days coadd-850_r1 snr-850_r1 crop-850_r1

#450_r0: reduce-450_r0-days cal-450_r0-days coadd-450_r0 snr-450_r0 crop-450_r0
#450_r1: reduce-450_r1-days cal-450_r1-days coadd-450_r1 snr-450_r1 crop-450_r1

#850_r1mf: reduce-850_r1mf-days cal-850_r1mf-days coadd-850_r1mf snr-850_r1mf crop-850_r1mf mf-850_r1mf-days


# define rules for day reductions
#
$(foreach reduction, $(dayreductions),\
    $(foreach day, $(days), \
        $(eval $(call Days_Template,$(reduction),$(day)))\
))

$(foreach reduction, $(dayreductions),\
    $(eval $(call Tag_Template,$(reduction)))\
)

.PHONY: ratios
.PHONY: ratio-j450_r0mf_j850_r0mf ratio-j450_r0mf_j850_r1mf


ratios: ratio-j450_r0mf_j850_r0mf ratio-j450_r0mf_j850_r1mf ratio-j450_r0mf_j850_r0_contamination_mf

$(eval rj450r0mfj850r0mf := $(RES_DIR)/ratios/ratio-j450_r0mf_j850_r0mf.sdf)

ratio-j450_r0mf_j850_r0mf: $(rj450r0mfj850r0mf)

$(rj450r0mfj850r0mf): $(CFG_DIR)/ratio-j450_r0mf_j850_r0mf.cfg $(RES_DIR)/j450_r0mf/$(SNAME)-j450_r0mf-reduc.sdf $(RES_DIR)/j850_r0mf/$(SNAME)-j850_r0mf-reduc.sdf
	sh $(BIN)/mkratios.sh \
		-c $(CFG_DIR)/ratio-j450_r0mf_j850_r0mf.cfg \
		-d $(RES_DIR)/ratios


$(eval rj450r0mfj850r1mf := $(RES_DIR)/ratios/ratio-j450_r0mf_j850_r1mf.sdf)

ratio-j450_r0mf_j850_r1mf: $(rj450r0mfj850r1mf)

$(rj450r0mfj850r1mf): $(CFG_DIR)/ratio-j450_r0mf_j850_r1mf.cfg $(RES_DIR)/j450_r0mf/$(SNAME)-j450_r0mf-reduc.sdf $(RES_DIR)/j850_r1mf/$(SNAME)-j850_r1mf-reduc.sdf
	sh $(BIN)/mkratios.sh \
		-c $(CFG_DIR)/ratio-j450_r0mf_j850_r1mf.cfg \
		-d $(RES_DIR)/ratios

$(eval rj450r0mfj850r0contmf := $(RES_DIR)/ratios/ratio-j450_r0mf_j850_r0_contamination_mf.sdf)

ratio-j450_r0mf_j850_r0_contamination_mf: $(rj450r0mfj850r0contmf)

$(rj450r0mfj850r0contmf): $(CFG_DIR)/ratio-j450_r0mf_j850_r0_contamination_mf.cfg $(RES_DIR)/j450_r0mf/$(SNAME)-j450_r0mf-reduc.sdf $(RES_DIR)/j850_r0_contamination_mf/$(SNAME)-j850_r0_contamination_mf-reduc.sdf
	sh $(BIN)/mkratios.sh \
		-c $(CFG_DIR)/ratio-j450_r0mf_j850_r0_contamination_mf.cfg \
		-d $(RES_DIR)/ratios


