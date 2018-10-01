##
##  Makefile to run scripts to reduce the JCMT SCUBA2 data of Sh2-61
##
##  O. Morata
##  2018
##

##-- Info --------------------------------------------------------------
HOME_DIR=/data/omorata/JCMT/SCUBA2/Sh2-61

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

reducedmaps = 850-r0


##-----
map850_r0: $(RES_DIR)/850-r0/Sh2_61-850-r0.sdf

$(RES_DIR)/850-r0/Sh2_61-850-r0.sdf: $(RES_DIR)/850-ro/gs20150424_32_850_reduced.sdf
	. $(BIN)/process_scuba2_maps.sh -c no -d $(RES_DIR)/850-r0 \
		-i gs20150424_32_850_reduced  \
		-o Sh2_61 -t -r0 -w 850 -x \
		-p $(CFG_DIR)/par850_process.lis		

$(RES_DIR)/850-r0/gs20150424_32_850_reduced.sdf: $(CFG_DIR)/run-850-r0
	. $(BIN)/reduce_raw.sh $(CFG_DIR)/run-850-r0


##-----

reduce-days-850: reduce-0424-850 reduce-0425-850

#.........

reduce-0424-850: $(RES_DIR)/850-r0/gs20150424_32_850_reduced.sdf

$(RES_DIR)/850-r0/gs20150424_32_850_reduced.sdf: $(CFG_DIR)/run-0424-850-r0
	. $(BIN)/reduce_raw.sh $(CFG_DIR)/run-0424-850-r0

#.........
reduce-0425-850: $(RES_DIR)/850-r0/gs20150425_59_850_reduced.sdf

$(RES_DIR)/850-r0/gs20150425_59_850_reduced.sdf: $(CFG_DIR)/run-0425-850-r0
	. $(BIN)/reduce_raw.sh $(CFG_DIR)/run-0425-850-r0

#.........

checkcal-0424:
	. $(BIN)/check_fcf.sh -d $(RES_DIR)/850-r0/ \
	    -i gs20150424_16_850_reduced.sdf

checkcal-0425:
	. $(BIN)/check_fcf.sh -d $(RES_DIR)/850-r0/ \
	    -i gs20150425_63_850_reduced.sdf

$(RES_DIR)/850-r0/Sh2_61-0424-850-r0.sdf: \
    $(RES_DIR)/850-r0/gs20150424_32_850_reduced.sdf
	. $(BIN)/process_scuba2_maps.sh -c y -d $(RES_DIR)/850-r0 \
		-i gs20150424_32_850_reduced  \
		-o Sh2_61 -t 850-r0  -x \
		-p $(CFG_DIR)/par850_process.lis		
