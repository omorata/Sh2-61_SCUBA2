#!/bin/bash
fclumps_infile=analysis_maps/Sh2_61-j850_r0_contamination_mf-snr.sdf
fclumps_method=FellWalker 
fclumps_config=fw-Sh2-61-j850_r0_contamination_mf-snr6.cfg
fclumps_deconv=FALSE
fclumps_logfile=Sh2-61-j850_r0_contamination_mf-fw-snr6.log
fclumps_wcspar=TRUE 
fclumps_shape=Polygon
fclumps_outfile=Sh2-61-j850_r0_contamination_mf-fw-snr6.sdf
fclumps_outcat=Sh2-61-j850_r0_contamination_mf-fw-snr6.fits
##
fclumps_extract=y
extract_data=analysis_maps/Sh2_61-j850_r0_contamination_mf.sdf
extract_out=Sh2-61-j850_r0_contamination_mf-fw-snr6-extr.sdf
extract_outcat=Sh2-61-j850_r0_contamination_mf-fw-snr6-extr.fits
extract_deconv=TRUE
extract_logfile=Sh2-61-j850_r0_contamination_mf-fw-snr6-extr.log
extract_shape=Polygon
extract_wcspar=TRUE
