##
##  Makefile to run scripts to reduce the JCMT SCUBA2 data of Sh2-61
##
##  O. Morata
##  2018
##

##-- Info --------------------------------------------------------------
HOME_DIR=/data/omorata/JCMT/SCUBA2/Sh2-61
SNAME=Sh2_61

#
##-- End info ----------------------------------------------------------

# names of directories
#
BIN=$(HOME_DIR)/scripts
CFG_DIR=$(HOME_DIR)/config
DATA_DIR=$(HOME_DIR)/data
RES_DIR=$(HOME_DIR)/results
TMP_DIR=$(HOME_DIR)/tmpdir


export

reducedmaps = 850_r0
reducedmaps = 850_r0 850_r1 850_r1mf 450_r0 450_r1 450_r1mf
days = 0424 0425


##-- Template definition -----------------------------------------------


# Template to define the rules to process steps of each different day for
# the given run
#
define Days_Template

reduce-$(1)-$(2) : $(RES_DIR)/$(1)/$(SNAME)-$(1)-$(2)-reduc.sdf

$(RES_DIR)/$(1)/$(SNAME)-$(1)-$(2)-reduc.sdf: $(CFG_DIR)/run-$(1)-$(2).cfg
	. $(BIN)/reduce_raw.sh $(CFG_DIR)/run-$(1)-$(2).cfg


cal-$(1)-$(2): $(RES_DIR)/$(1)/$(SNAME)-$(1)-$(2)-reduc_cal.sdf

$(RES_DIR)/$(1)/$(SNAME)-$(1)-$(2)-reduc_cal.sdf: $(RES_DIR)/$(1)/$(SNAME)-$(1)-$(2)-reduc.sdf
	$(BIN)/post_scuba2.sh  -a cal  -d $(RES_DIR)/$(1)/ \
		-p $(CFG_DIR)/par$(1)-$(2).cfg  -i $(SNAME)-$(1)-$(2)-reduc


checkcal-$(1)-$(2):
	. $(BIN)/check_fcf.sh -d $(RES_DIR)/$(1) -c $(CFG_DIR)/run-$(1)-$(2).cfg

endef



#define Join_Template
#reduce-$(1)-days: reduce-$(1)-0424 reduce-$(1)-0425

#cal-$(1)-days: cal-$(1)-0424 cal-$(1)-0425

#$(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd.sdf:  $(RES_DIR)/$(1)/$(SNAME)-$(1)-0424-reduc_cal.sdf $(RES_DIR)/$(1)/$(SNAME)-$(1)-0425-reduc_cal.sdf   $(CFG_DIR)/$(1)-mosaic.list
#	$(BIN)/post_scuba2.sh  -a add  -d $(RES_DIR)/$(1) \
#		-p $(CFG_DIR)/par850_process.lis \
#		-l $(CFG_DIR)/$(1)-mosaic.list  -o $(SNAME)-$(1)-coadd


#checkcal-$(1): checkcal-$(1)-0424 checkcal-$(1)-0425

#endef



# Template to run the rest of the steps of the defined run
#
define Tag_Template
reduce-$(1)-days: $(reduce_$(1)_list)

cal-$(1)-days: $(cal_$(1)_list)

checkcal-$(1): $(checkcal_$(1)_list)

coadd-$(1): $(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd.sdf

$(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd.sdf:  $(jfiles_$(1)_list)  $(CFG_DIR)/$(1)-mosaic.list
	$(BIN)/post_scuba2.sh  -a add  -d $(RES_DIR)/$(1) \
		-p $(CFG_DIR)/par$(1).cfg \
		-l $(CFG_DIR)/$(1)-mosaic.list  -o $(SNAME)-$(1)-coadd


snr-$(1): $(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd_snr.sdf

$(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd_snr.sdf: $(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd.sdf
	$(BIN)/post_scuba2.sh -a snr  -d $(RES_DIR)/$(1)  \
		-i $(SNAME)-$(1)-coadd


crop-$(1): $(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd_crop.sdf

$(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd_crop.sdf: $(RES_DIR)/$(1)/$(SNAME)-$(1)-coadd.sdf
	$(BIN)/post_scuba2.sh  -a crop  -d $(RES_DIR)/$(1) \
		-p $(CFG_DIR)/par$(1).cfg  -i $(SNAME)-$(1)-coadd



endef


##-- End of template definition ----------------------------------------

##-- Generate lists ----------------------------------------------------

# loop to generate the lists needed for the rules that join data from more
# than one day
#
$(foreach run, $(reducedmaps),\
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


$(foreach run, $(reducedmaps),\
    $(foreach day, $(days), \
        $(eval $(call Days_Template,$(run),$(day)))\
))

$(foreach run, $(reducedmaps),\
    $(eval $(call Tag_Template,$(run)))\
)

