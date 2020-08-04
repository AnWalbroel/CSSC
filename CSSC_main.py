import numpy as np
import os
import datetime
import glob
import netCDF4 as nc
import pdb
import pandas as pd
from concat_mwr import run_concat_mwr
from dropsonde_raw_gap_filler import run_dropsonde_raw_gap_filler
from get_sst_data import download_sst_data
from HALO_raw_dropsonde_to_TB import run_HALO_raw_dropsonde_to_TB
from TB_statistics_raw import run_TB_statistics_raw



'''
	This is the control room for the Clear Sky Sonde Comparison (CSSC) code package where
	you assign the required variables (mostly paths of data).
	More information about each to be manually assigned variable can be found in 
	'clear_sky_sonde_comparison_README.txt'.
'''


# # # # # # # # 	1. "concat_mwr.py": Due to changing conventions how the MWR data is saved you have to 
# # # # # # # # 	manually edit 3 lines in "concat_mwr.py" (marked with "##################" and 
# # # # # # # # 	explained in the README file).
# # # # # # # mwr_base_path = "/data/hamp/flights/EUREC4A/"
# # # # # # # mwr_out_path = "/work/walbroel/CSSC_test/halo/mwr/"
# # # # # # # print("Running concat_mwr.py ..........\n")
# # # # # # # run_concat_mwr(mwr_base_path, mwr_out_path)


# # # # # # # # 	2. "dropsonde_raw_gap_filler.py":
# # # # # # # path_raw_sondes = "/work/walbroel/CSSC_test/sonde_half_raw/"
# # # # # # # data_out_path_halo = "/work/walbroel/CSSC_test/sonde_half_raw_v01/"
# # # # # # # path_BAH_data = "/data/hamp/flights/EUREC4A/unified/"	# BAHAMAS data path
# # # # # # # print("Running dropsonde_raw_gap_filler.py ..........\n")
# # # # # # # run_dropsonde_raw_gap_filler(path_raw_sondes, data_out_path_halo, path_BAH_data)


# # # # # # # #	3. Either manually select and download the SST data from the website given in
# # # # # # # #	the README file or assign latitude and longitude boundaries and required dates 
# # # # # # # #	for an automated download:
# # # # # # # path_sst = "/work/walbroel/CSSC_test/sst_slice/"
# # # # # # # lat_bound = [0, 40]
# # # # # # # lon_bound = [-70, 0]
# # # # # # # start_date = "2020-01-19"
# # # # # # # end_date = "2020-02-18"
# # # # # # # print("Running get_sst_data.py ..........\n")
# # # # # # # download_sst_data(path_sst, lat_bound, lon_bound, start_date, end_date)


# # # # # # # #	4. "HALO_raw_dropsonde_to_TB.py":
# # # # # # # path_halo_dropsonde = "/work/walbroel/CSSC_test/sonde_half_raw_v01/"	# interpolated dropsonde output path
# # # # # # # path_sst_data = "/work/walbroel/CSSC_test/sst_slice/"		# SST data path
# # # # # # # path_BAH_data = "/data/hamp/flights/EUREC4A/unified/"	# BAHAMAS data path
# # # # # # # pam_out_path = "/work/walbroel/CSSC_test/pamtra_half_raw_PQC/"
# # # # # # # pamtra_datadir = '/work/walbroel/pamtra/'
# # # # # # # print("Running HALO_raw_dropsonde_to_TB.py ..........\n")
# # # # # # # run_HALO_raw_dropsonde_to_TB(path_halo_dropsonde, path_sst_data, path_BAH_data, pam_out_path, pamtra_datadir)


#	5. "TB_statistics_raw.py":
path_mwr = "/work/walbroel/CSSC_test/halo/mwr/"		# path of the concatenated mwr files
path_pam_ds = "/work/walbroel/CSSC_test/pamtra_half_raw_PQC/"	# output path of pamtra
path_BAH_data = "/data/hamp/flights/EUREC4A/unified/"	# BAHAMAS data path
out_path = "/work/walbroel/CSSC_test/sonde_comparison_half_raw/"
plot_path = "/work/walbroel/CSSC_test/plots/TB_stat/"
scatterplot_name = "TB_scatterplot"
bias_ev_plotname = "TB_abs_biasevolution"
output_filename = "clear_sky_sonde_comparison_ALL"
print("Running TB_statistics_raw.py ..........\n")
run_TB_statistics_raw(path_mwr, path_pam_ds, path_BAH_data, out_path, plot_path, scatterplot_name, 
		bias_ev_plotname, output_filename)