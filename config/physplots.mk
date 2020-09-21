#  physplots.mk
#  O. Morata 2020
#
#  rules to produce plots out of the physical parameters calculated on
#  the JCMT SCUBA2 data of Sh2-61
#
#  To be included by the main Makefile


define HistoPlots
# Template to plot histograms from the pixel maps of calculated physical
# parameters
#
# Arguments: 1- tgt, 2- findclump id, 3- physcalc tag, 4-histoplot tag
#
# The plots depend on calcs.log as a generic way  and on the wildcard of
# a .yml file with the same stem as the pdf
#
# The output name must have histo between arg 3 and arg 4
#
# models
#  Sh2_61-j850r0_co_mb__j450ro_mb-fw_01_ref-histo-N_all.pdf
#
endef



define XYPlots
# Template to plot xyplots from the pixel maps of calculated physical
# parameters
#
# Arguments: 1- tgt, 2- findclump id, 3- physcalc tag, 4-xyplot tag
#
# the plots depend on calcs.log as a generic way  and on the wildcard of
# a .yml file with the same stem as the pdf
#
# The output name must have xyplot between arg 3 and arg 4
#
# models
#  Sh2_61-j850r0_co_mb__j450ro_mb-fw_01_ref-xyplot-N_vs_tdust_all.pdf
#
endef
