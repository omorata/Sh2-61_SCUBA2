#  help.mk
#  O. Morata 2020
#
#  rules for the help information of the Makefile used to reduce and
#  analyse the JCMT SCUBA2 data of Sh2-61
#
#  To be included by the main Makefile
#

help:
	@echo
	@echo " ----------------------------------------------------------"
	@echo "  Makefile to reduce and analyse the data of $(PRJ_NAME)"
	@echo " ----------------------------------------------------------"
	@echo "  Pre-defined variables:"
	@echo "      Project Name : $(PRJ_NAME)"
	@echo "      File prefix : $(SNAME)"
	@echo "      targets : $(targets)"
	@echo "      findclumps ids : $(fcs)"
	@echo "      combined targets : $(combined)"
	@echo "      physical parameter calculation tags : $(physcalc_tags)"
	@echo "      physical parameters maps : $(comb_maps)"
	@echo "      histogram plot tags : $(histo_tags)"
	@echo "      xyplot tags : $(xyplots_tags)"
	@echo;echo;echo "  Help options:"
	@echo "            make help - this help"
	@echo "       make help_dirs - information on directories"
	@echo "      make help_rules - information on defined rules"
	@echo;echo;echo "  List options:"
	@echo "            make list - list all the rules"
	@echo "      make list_files - list all the files (possibly) created"
	@echo "       make list_deps - list the rules and their dependencies"
	@echo;echo



help_dirs:
	@echo
	@echo " ---------------------------------------------------------"
	@echo "  Directory set-up for project $(PRJ_NAME):"
	@echo " ---------------------------------------------------------"
	@echo "    project home : $(HOME_DIR)"
	@echo "            bin  : $(BIN_DIR)"
	@echo "   configuration : $(CFG_DIR)"
	@echo "            data : $(DATA_DIR)"
	@echo "         results : $(RES_DIR)"
	@echo "    external bin : $(EXT_DIR)"
	@echo;echo



help_rules:
	@echo;echo " -----------------"
	@echo "   Defined rules"
	@echo " -----------------"
	@echo
	@echo "  (Warning: not all the following rules may be available."
	@echo "   In many cases, it will depend on the definition of the"
	@echo "   corresponding configuration files)"
	@echo;echo " The general actions are:"
	@echo;echo "    make reduce -- do the reduction of the targets"
	@echo
	@echo "    make plotmaps  --  plot maps of targets"
	@echo "    make strip  --  strip third axis from sdf files of targets"
	@echo "    make tofits  --  transform .sdf files of targets to .fits"
	@echo "    make tofits-strip  --  transform stripped files to .fits "
	@echo "    [ make align ] -- align map to reference"
	@echo "    make findclumps_snr  --  find clumps in map using snr map"
	@echo "    make polygonfiles  --  extract shapes of clumps"
	@echo "    make clumps-map  --  plot overlay of clumps on emission map"
	@echo "    [ make calcs ] -- calculate physical parameters from" \
            "emission and clumps"
	@echo "    make maps-physpar  -- plot map of physical parameters"
	@echo;echo "   clean options: clean-reduce clean-plotmaps clean-strip" \
            "clean-fits_datasets"
	@echo "        clean-findclumps clean-clumps-map clean-align" \
            "clean-maps-physpar clean-calcs"
	@echo;echo " In more detail:"
	@echo;echo " + rules depending on the target"
	@echo "    (targets: $(targets))"
	@echo;echo "     make reduce-[target] -- includes the following:"
	@echo "         make map-[target]  --  get the reduced map"
	@echo "         make snr-[target]  --  get the snr map"
	@echo "         make crop-[target]  -- get the cropped map"
	@echo "         make snrcrop-[target]  -- get the cropped snr map"
	@echo;echo "     make plotmap-[target]"
	@echo "     make strip-[target]"
	@echo "     make tofits-[target]"
	@echo "     make tofits_strip-[target]"
	@echo "     make findclumps_snr-[target]"
	@echo "     make polygonfiles-[target]"
	@echo "     make clumps-map-[target]"
	@echo;echo "    clean options: clean-plotmap-[target]" \
            "clean-strip-[target]"
	@echo "         clean-fits_datasets-[target] clean-findclumps-[target]"
	@echo "         clean-clumps-map-[target] clean-calcs"
	@echo
	@echo;echo " + rules depending on the target and findclumps_id"
	@echo "    (ids: $(fcs))"
	@echo;echo "     make findclumps_snr-[target]-[id]"
	@echo "     make polygonfiles-[target]-[id]"
	@echo "     make clumps-map-[target]-[id]"
	@echo;echo "    clean options: clean-findclumps-[target]-[id]" \
            "clean-clumps-map-[target]-[id]"
	@echo
	@echo;echo " + rules depending on the combination map"
	@echo "    (comb. maps: $(combined))"
	@echo;echo "     make align-[secondary_target]-to-[reference_target]"
	@echo "     make calcs-[comb_map]"
	@echo "     make calcs-[reference_target]"
	@echo "     make maps-physpar-[comb.map]"
	@echo;echo "    clean options: clean-align-[comb_map]" \
            "clean-align-[secondary-target]"
	@echo "         clean-calcs-[comb_map] clean-calcs-[reference_target]"
	@echo "         clean-maps-physpar-[comb.map]"
	@echo
	@echo;echo " + rules depending on combination map and findclump_id"
	@echo;echo "     make calcs-[combmap]-[id]"
	@echo "     make maps-physpar-[combmap]-[id]"
	@echo;echo "    clean options: clean-calcs-[combmap]-[id]" \
            "clean-maps-physpar-[combmap]-[id]"
	@echo
	@echo;echo " + rules depending on the comb. map, findclumps_id," \
            "and physical parameter"
	@echo "    (phys. param: $(comb_maps))"
	@echo;echo "     make map-physpar-[comb.map]-[id]-[physpar]"
	@echo
	@echo "    clean options:  clean-map-physpar[comb.map]-[id]-[physpar]"
	@echo;echo

.PHONY: help help_dirs help_rules
