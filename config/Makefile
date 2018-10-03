##
##  Makefile to run scripts to reduce the JCMT SCUBA2 data of Sh2-61
##
##  O. Morata
##  2018
##

##-- Info --------------------------------------------------------------
HOME_DIR=/data/omorata/JCMT/SCUBA2/Sh2-61
SNAME=Sh2_61

dayruns = 850_r0 850_r1 450_r0 450_r1 850_r1mf 450_r1mf
mfdayruns = 850_r1mf 450_r1mf
days = 0424 0425

jointruns = j850_r0 j850_r1 j450_r0 j450_r1

#
##-- End info ----------------------------------------------------------

# names of directories
#
BIN=$(HOME_DIR)/scripts
CFG_DIR=$(HOME_DIR)/config
DATA_DIR=$(HOME_DIR)/data
RES_DIR=$(HOME_DIR)/results


export


##-- Template definition -----------------------------------------------

# Template to define the rules to process steps of each different day for
# the given run
#
define Days_Template

.PHONY: reduce-$(1)-$(2) cal-$(1)-$(2) checkcal-$(1)-$(2)


reduce-$(1)-$(2) : $(RES_DIR)/$(1)/$(SNAME)-$(1)-$(2)-reduc.sdf

cal-$(1)-$(2): $(RES_DIR)/$(1)/$(SNAME)-$(1)-$(2)-reduc_cal.sdf


$(RES_DIR)/$(1)/$(SNAME)-$(1)-$(2)-reduc.sdf: $(CFG_DIR)/run-$(1)-$(2).cfg
	. $(BIN)/reduce_raw.sh $(CFG_DIR)/run-$(1)-$(2).cfg


$(RES_DIR)/$(1)/$(SNAME)-$(1)-$(2)-reduc_cal.sdf: $(RES_DIR)/$(1)/$(SNAME)-$(1)-$(2)-reduc.sdf
	$(BIN)/post_scuba2.sh  -a cal  -d $(RES_DIR)/$(1)/ \
		-p $(CFG_DIR)/par$(1)-$(2).cfg  -i $(SNAME)-$(1)-$(2)-reduc


checkcal-$(1)-$(2):
	. $(BIN)/check_fcf.sh -d $(RES_DIR)/$(1) -c $(CFG_DIR)/run-$(1)-$(2).cfg

endef



# Template to run the rest of the steps of the defined run
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


$(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd.sdf:  $(jfiles_$(1)_list)  $(CFG_DIR)/$(1)-mosaic.list
	$(BIN)/post_scuba2.sh  -a add  -d $(RES_DIR)/$(1) \
		-p $(CFG_DIR)/par$(1).cfg \
		-l $(CFG_DIR)/$(1)-mosaic.list  -o $(SNAME)-$(1)-coadd


$(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd_snr.sdf: $(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd.sdf
	$(BIN)/post_scuba2.sh -a snr  -d $(RES_DIR)/$(1)  \
		-i $(SNAME)-$(1)-coadd


$(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd_crop.sdf: $(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd.sdf
	$(BIN)/post_scuba2.sh  -a crop  -d $(RES_DIR)/$(1) \
		-p $(CFG_DIR)/par$(1).cfg  -i $(SNAME)-$(1)-coadd

endef



# Template to process all days of the same run on the same step
#
define Joint_Template

.PHONY: $(1) reduce-$(1) snr-$(1) crop-$(1)


$(1): reduce-$(1) snr-$(1) crop-$(1)

reduce-$(1) : $(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc.sdf

snr-$(1): $(RES_DIR)/$(1)/$(SNAME)-$(1)-snr.sdf

crop-$(1): $(RES_DIR)/$(1)/$(SNAME)-$(1)-crop.sdf


$(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc.sdf: $(CFG_DIR)/run-$(1).cfg
	. $(BIN)/reduce_raw.sh $(CFG_DIR)/run-$(1).cfg


$(RES_DIR)/$(1)/$(SNAME)-$(1)-snr.sdf: $(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc.sdf
	$(BIN)/post_scuba2.sh -a snr  -d $(RES_DIR)/$(1)  \
		-i $(SNAME)-$(1)-reduc


$(RES_DIR)/$(1)/$(SNAME)-$(1)-crop.sdf: $(RES_DIR)/$(1)/$(SNAME)-$(1)-reduc.sdf
	$(BIN)/post_scuba2.sh  -a crop  -d $(RES_DIR)/$(1) \
		-p $(CFG_DIR)/par$(1).cfg  -i $(SNAME)-$(1)-reduc

endef


##-- End of template definition ----------------------------------------

##-- Generate lists ----------------------------------------------------

# loop to generate the lists needed for the rules that join data from more
# than one day
#
$(foreach run, $(dayruns),\
    $(eval reduce_$(run)_list = ) \
    $(eval cal_$(run)_list = ) \
    $(eval checkcal_$(run)_list = ) \
    $(eval jfiles_$(run)_list = ) \
    $(foreach day, $(days),\
        $(eval reduce_$(run)_list += reduce-$(run)-$(day)) \
        $(eval cal_$(run)_list += cal-$(run)-$(day)) \
        $(eval checkcal_$(run)_list += checkcal-$(run)-$(day)) \
        $(eval jfiles_$(run)_list += $(RES_DIR)/$(run)/$(SNAME)-$(run)-$(day)-reduc_cal.sdf) \
))

##-- Enf of lists ------------------------------------------------------


850_r0: reduce-850_r0-days cal-850_r0-days coadd-850_r0 snr-850_r0 crop-850_r0
850_r1: reduce-850_r1-days cal-850_r1-days coadd-850_r1 snr-850_r1 crop-850_r1

450_r0: reduce-450_r0-days cal-450_r0-days coadd-450_r0 snr-450_r0 crop-450_r0
450_r1: reduce-450_r1-days cal-450_r1-days coadd-450_r1 snr-450_r1 crop-450_r1

850_r1mf: reduce-850_r1mf-days cal-850_r1mf-days coadd-850_r1mf snr-850_r1mf crop-850_r1mf mf-850_r1mf-days


# define rules for day runs
#
$(foreach run, $(dayruns),\
    $(foreach day, $(days), \
        $(eval $(call Days_Template,$(run),$(day)))\
))

$(foreach run, $(dayruns),\
    $(eval $(call Tag_Template,$(run)))\
)


# define rules for joint runs
#
$(foreach run, $(jointruns),\
    $(eval $(call Joint_Template,$(run)))\
)

