#  findclumps.mk
#  O. Morata 2020
#
#  rules to look for clumps with Starlink's findclump command in the
#  JCMT SCUBA2 data of Sh2-61
#
#  To be included by the main Makefile



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
        fi
        #else \
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
        fi
        #else \
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





# print catalog of findclumps parameters 
#
$(eval table_name := $(findclumps_dir)/$(SNAME)-$(1)-$(2)-catalog)

$(table_name).txt: $$(wilcard $$(table_name).fits)
	@if [ -f $(table_name).fits ]; then \
            $(BIN_DIR)/print_catalog.py \
                -t 'findclumps' \
                -i $(table_name).fits \
                -o $(table_name).txt \
                -f "['PIDENT', 'Peak1', 'Peak2', 'Cen1', 'Cen2', 'Size1', \
                    'Size2', 'Sum', 'Peak', 'Volume']"; \
        fi


print-clumpcatg-$(1)-$(2): $(table_name).txt

.PHONY: print-clumpcatg-$(1)-$(2)


print-clumpcatg-$(1): print-clumpcatg-$(1)-$(2)

.PHONY: print-clumpcatg-$(1)


print-clumpcatg: print-clumpcatg-$(1)

.PHONY: print-clumpcatg


print-clumpcatg-$(ref_tgt): print-clumpcatg-$(1)

.PHONY: print-clumpcatg-$(ref_tgt)



clean-clumpcatg-$(1)-$(2):
	@rm -fv $(table_name).txt

.PHONY: clean-clumpcatg-$(1)-$(2)


clean-clumpcatg-$(1): clean-clumpcatg-$(1)-$(2)

.PHONY: clean-clumpcatg-$(1)


clean-clumpcatg: clean-clumpcatg-$(1)

.PHONY: clean-clumpcatg


clean-clumpcatg-$(ref_tgt): clean-clumpcatg-$(1)

.PHONY: clean-clumpcatg-$(ref_tgt)



endef
