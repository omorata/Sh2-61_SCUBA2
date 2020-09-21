#  reduction.mk
#  O. Morata 2020
#
#  rules to reduce the JCMT SCUBA2 data of Sh2-61
#
#  To be included by the main Makefile
#

define DoReduction
# Template to reduce all observations in a single step
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
.PHONY: crop-$(1)

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

clean-reduction-$(1): clean-map-$(1) clean-snr-$(1) clean-crop-$(1)
clean-reduction-$(1): clean-snrcrop-$(1)

.PHONY: clean-reduction-$(1)

clean-reduction: clean-reduction-$(1)
.PHONY: clean-reduction

endef
