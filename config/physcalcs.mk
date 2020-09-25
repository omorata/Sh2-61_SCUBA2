#  physcalcs.mk
#  O. Morata 2020
#
#  rules to calculate the physical parameters dervied from the JCMT
#  SCUBA2 850 and 450 micron data of Sh2-61
#
#  To be included by the main Makefile



define Align_Dataset
# Template to align one dataset to the grid of another one
# 
# 1- combined string

$(eval ref_tgt := $(firstword $(subst __, ,$(1))))
$(eval sec_tgt := $(lastword $(subst __, ,$(1))))

$(eval outdir := $(RES_DIR)/analysis_maps)


$(eval align_ref := $(outdir)/$(SNAME)-$(ref_tgt).sdf)

$(eval align_tgt := $(outdir)/$(SNAME)-$(sec_tgt).sdf)
$(eval align_snrtgt := $(outdir)/$(SNAME)-$(sec_tgt)-snr.sdf)

$(eval aligned_file :=  \
    $(outdir)/$(SNAME)-$(sec_tgt)-aligned_to-$(ref_tgt).sdf)
$(eval aligned_snrfile :=    \
    $(outdir)/$(SNAME)-$(sec_tgt)-snr-aligned_to-$(ref_tgt).sdf)

$(eval aligned_fits :=  \
    $(outdir)/$(SNAME)-$(sec_tgt)-aligned_to-$(ref_tgt).fits)
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
# 1- combined string, 2- findclump_id, 3- physical parameters
# calculation variant
#

$(eval ref_tgt := $(firstword $(subst __, ,$(1))))
$(eval sec_tgt := $(lastword $(subst __, ,$(1))))

$(eval outdir := $(RES_DIR)/analysis_maps)


$(eval aligned_fits :=  \
    $(outdir)/$(SNAME)-$(sec_tgt)-aligned_to-$(ref_tgt).fits)
$(eval aligned_snrfits :=    \
    $(outdir)/$(SNAME)-$(sec_tgt)-snr-aligned_to-$(ref_tgt).fits)


#
# Calculation of physical parameters
#

$(eval ffile := $(outdir)/$(SNAME)-$(ref_tgt).fits)
$(eval ffile_snr := $(outdir)/$(SNAME)-$(ref_tgt)-snr.fits)
$(eval clfile := $(RES_DIR)/findclumps/$(SNAME)-$(ref_tgt)-$(2)-clumps.fits)

$(eval cfg_file := $(CFG_DIR)/analysis/$(SNAME)-$(1)-$(2)_$(3)-phys_calc.yaml)

$(eval calc_log := $(outdir)/calcs-$(1)-$(2)_$(3).log)

$(eval calc_refs :=    \
    $(ffile) $(ffile_snr) $(aligned_fits) $(aligned_snrfits) $(clfile))

$(calc_log):  $$(wildcard $$(cfg_file)) $(calc_refs)
	@if [ -f $(cfg_file) ]; then \
	    $(BIN_DIR)/calc_phys_parm.py \
                -c $(cfg_file) \
                --clumps \
                -l $(calc_log) \
                -w $(RES_DIR);\
        fi
          #else \
          #   echo -e "\n++ Ignoring rule $(out_fc)" ;\
          #   echo -e "    No cfg file $(cfg_file)" ;\
          #fi


.PHONY: calcs-$(1)-$(2)_$(3)
calcs-$(1)-$(2)_$(3) : $(calc_log)

.PHONY: calcs-$(1)-$(2)
calcs-$(1)-$(2) : calcs-$(1)-$(2)_$(3)

.PHONY: calcs-$(1)
calcs-$(1): calcs-$(1)-$(2)

.PHONY: calcs-$(ref_tgt)
calcs-$(ref_tgt) : calcs-$(1)


.PHONY: clean-calcs-$(1)-$(2)_$(3)
clean-calcs-$(1)-$(2)_$(3):
	@rm -fv $(calc_log)
	@rm -fv $(outdir)/$(SNAME)-$(1)-$(2)_$(3)-ratio.fits
	@rm -fv $(outdir)/$(SNAME)-$(1)-$(2)_$(3)-tdust.fits
	@rm -fv $(outdir)/$(SNAME)-$(1)-$(2)_$(3)-mass.fits
	@rm -fv $(outdir)/$(SNAME)-$(1)-$(2)_$(3)-N.fits
	@rm -fv $(outdir)/$(SNAME)-$(1)-$(2)_$(3)-clump_table.fits

.PHONY: clean-calcs-$(1)-$(2)
clean-calcs-$(1)-$(2): clean-calcs-$(1)-$(2)_$(3)

.PHONY: clean-calcs-$(1)
clean-calcs-$(1): clean-calcs-$(1)-$(2)

.PHONY: clean-calcs
clean-calcs: clean-calcs-$(1)

.PHONY: clean-calcs-$(ref_tgt)
clean-calcs-$(ref_tgt): clean-calcs-$(1)


# print catalog of physical parameters of clumps
#
$(eval table_name := $(outdir)/$(SNAME)-$(1)-$(2)_$(3)-clump_table)

$(table_name).txt: $$(wilcard $$(table_name).fits) $(calc_log)
	@if [ -f $(table_name).fits ]; then \
            $(BIN_DIR)/print_catalog.py \
                -t 'phys' \
                -i $(table_name).fits \
                -o $(table_name).txt; \
        fi

print_physcatg-$(1)-$(2)_$(3): $(table_name).txt

.PHONY: print_physcatg-$(1)-$(2)_$(3)


print_physcatg-$(1)-$(2): print_physcatg-$(1)-$(2)_$(3)

.PHONY: print_physcatg-$(1)-$(2)


print_physcatg-$(1): print_physcatg-$(1)-$(2)

.PHONY: print_physcatg-$(1)


print_physcatg: print_physcatg-$(1)

.PHONY: print_physcatg


print_physcatg-$(ref_tgt): print_physcatg-$(1)

.PHONY: print_physcatg-$(ref_tgt)



clean-physcatg-$(1)-$(2)_$(3):
	@rm -fv $(table_name).txt

.PHONY: clean-physcatg-$(1)-$(2)_$(3)


clean-physcatg-$(1)-$(2): clean-physcatg-$(1)-$(2)_$(3)

.PHONY: clean-physcatg-$(1)-$(2)


clean-physcatg-$(1): clean-physcatg-$(1)-$(2)

.PHONY: clean-physcatg-$(1)


clean-physcatg: clean-physcatg-$(1)

.PHONY: clean-physcatg


clean-physcatg-$(ref_tgt): clean-physcatg-$(1)

.PHONY: clean-physcatg-$(ref_tgt)


endef
