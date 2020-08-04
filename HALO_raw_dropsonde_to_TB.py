import numpy as np
import os
import datetime
import glob
from general_importer import readNCraw_V01, import_GHRSST, import_BAHAMAS_unified		# "General Importer reporting for duty"
import pdb
import matplotlib.pyplot as plt



def run_HALO_raw_dropsonde_to_TB(path_halo_dropsonde, path_sst_data, path_BAH_data, pam_out_path, pamtra_datadir):

	os.environ['PAMTRA_DATADIR'] = pamtra_datadir
	import pyPamtra

	HALO_sondes_NC = sorted(glob.glob(path_halo_dropsonde + "*v01.nc"))
	SST_files_NC = sorted(glob.glob(path_sst_data + "*.nc.nc4"))
	BAH_files_NC = sorted(glob.glob(path_BAH_data + "bahamas*.nc"))

	for filename_in in HALO_sondes_NC:
		# Import sonde data:
		# # filename_in = HALO_sondes_NC[0]
		sonde_dict = readNCraw_V01(filename_in)

		# Import GHRSST CMC data:
		convention_timedelta_sec = (datetime.datetime(2020,1,1) - datetime.datetime(1970,1,1)).total_seconds()
		dropsonde_date = datetime.datetime.utcfromtimestamp(sonde_dict['launch_time'] + convention_timedelta_sec).strftime("%Y%m%d")
		sst_filename = [sst_file for idx, sst_file in enumerate(SST_files_NC) if dropsonde_date in sst_file]
		sst_keys = ['time', 'lat', 'lon', 'analysed_sst', 'analysis_error']
		sst_dict = import_GHRSST(sst_filename[0], sst_keys)

		# Import altitude and time data from BAHAMAS:
		bah_filename = [bah_file for idx, bah_file in enumerate(BAH_files_NC) if dropsonde_date in bah_file]
		bah_keys = ['time', 'altitude']
		bah_dict = import_BAHAMAS_unified(bah_filename[0], bah_keys)
		bah_dict['time'] = np.rint(bah_dict['time']).astype(float)		# must be done to avoid small fractions of seconds
		bah_dict['time'] = (datetime.datetime(1970,1,1) - datetime.datetime(2020,1,1)).total_seconds() + bah_dict['time']	# convert to seconds since 2020-01-01 00:00:00
		

		n_launches = 1		# number of sonde launches
		n_alt = len(sonde_dict['Z'])		# number of height levels

		# find the sonde launches that produced too many nan values so that cannot run: use the RH, T, P for that:
		if not (np.all([~np.isnan(sonde_dict['T']), ~np.isnan(sonde_dict['P']), ~np.isnan(sonde_dict['RH'])])):
			continue

		# assert np.all(~np.isnan(sonde_dict['RH']))


		# HAMP FREQUENCIES:
		frq = [22.2400,23.0400,23.8400,25.4400,26.2400,27.8400,31.4000,50.3000,51.7600,52.8000,53.7500,54.9400,56.6600,58.0000,90.0000,110.250,114.550,116.450,117.350,120.150,121.050,122.950,127.250,170.810,175.810,178.310,179.810,180.810,181.810,182.710,183.910,184.810,185.810,186.810,188.310,190.810,195.810]

		# create pamtra object; change settings:
		pam = pyPamtra.pyPamtra()

		pam.nmlSet['hydro_adaptive_grid'] = True
		pam.nmlSet['add_obs_height_to_layer'] = False		# adds observation layer height to simulation height vector
		pam.nmlSet['passive'] = True						# passive simulation
		pam.nmlSet['active'] = False						# False: no radar simulation
		# pam.nmlSet['gas_mod'] = 'L93'						# default: 'R98'

		pamData = dict()
		shape2d = [1, 1]

		# use highest non nan values of sonde for location information:
		if ~np.isnan(sonde_dict['reference_lon']):
			reflon = sonde_dict['reference_lon']
		else:
			reflon = np.asarray([sonde_dict['lon'][np.argwhere(~np.isnan(sonde_dict['lon']))[-1]] ]).flatten()

		if ~np.isnan(sonde_dict['reference_lat']):
			reflat = sonde_dict['reference_lat']
		else:
			reflat = np.asarray([sonde_dict['lat'][np.argwhere(~np.isnan(sonde_dict['lat']))[-1]] ]).flatten()

		pamData['lon'] = np.broadcast_to(reflon, shape2d)
		pamData['lat'] = np.broadcast_to(reflat, shape2d)
		pamData['timestamp'] = np.broadcast_to(sonde_dict['launch_time'], shape2d)

		# to get the obs_height: average BAHAMAS altitude over +/- 10 seconds around launch_time:
		# find time index of the sonde launches:
		bah_launch_idx = np.asarray([np.argwhere(bah_dict['time'] == pamData['timestamp'][0])]).flatten()		# had some dimensions too many -> flattened
		drop_alt = np.floor(np.asarray([np.mean(bah_dict['altitude'][i-10:i+10]) for i in bah_launch_idx])/100)*100
		obs_height = drop_alt

		print(dropsonde_date + "\n")

		# surface type & reflectivity:
		pamData['sfc_type'] = np.zeros(shape2d)			# 0: ocean, 1: land
		pamData['sfc_refl'] = np.chararray(shape2d)
		pamData['sfc_refl'][:] = 'F'
		pamData['sfc_refl'][pamData['sfc_type'] == 1] = 'S'
		
		pamData['obs_height'] = np.broadcast_to(obs_height, shape2d + [len(obs_height), ])			# must be changed to the actual top of soundings (or mwr altitude / bahamas altitude)

		# meteorolog. surface information:
		# to find the SST: use the designated lat,lon in pamData to find the closest entry in the GHRSST dataset:
		dlat = np.asarray([sst_dict['lat'] - pamData['lat'][0,0]])
		dlon = np.asarray([sst_dict['lon'] - pamData['lon'][0,0]])	# for each sonde, for each sst entry
		distance_lat_squared = (2*np.pi*6371000*dlat/360)**2
		distance_lon_squared = (2*np.pi*6371000*np.cos(2*np.pi*pamData['lat'][0,0]/360)*dlon/360)**2
		i_lat = np.argmin(distance_lat_squared)	# contains index of sst_dict['lat'] which had a min distance to pamData[lat] for each sonde
		i_lon = np.argmin(distance_lon_squared)	# contains index of sst_dict['lat'] which had a min distance to pamData[lat] for each sonde

		sst = sst_dict['SST'][0,i_lat,i_lon]		# [time, lat, lon]


		pamData['groundtemp'] = np.broadcast_to(sst, shape2d)
		pamData['wind10u'] = np.broadcast_to(sonde_dict['u_wind'][1], shape2d)		# = 1 because sfc would be 0 m; index 1 is 10 m
		pamData['wind10v'] = np.broadcast_to(sonde_dict['v_wind'][1], shape2d)		# = 1 because sfc would be 0 m; index 1 is 10 m

		# 3d variables:
		shape3d = shape2d + [n_alt]
		pamData['hgt_lev'] = np.broadcast_to(sonde_dict['Z'], shape3d)
		pamData['temp_lev'] = np.broadcast_to(sonde_dict['T'][:], shape3d)
		pamData['press_lev'] = np.broadcast_to(sonde_dict['P'][:], shape3d)
		pamData['relhum_lev'] = np.broadcast_to(sonde_dict['RH'][:], shape3d)

		# 4d variables: hydrometeors:
		shape4d = [1, 1, n_alt-1, 5]			# potentially 5 hydrometeor classes with this setting
		pamData['hydro_q'] = np.zeros(shape4d)
		pamData['hydro_q'][...,0] = 0# CLOUD
		pamData['hydro_q'][...,1] = 0# ICE
		pamData['hydro_q'][...,2] = 0# RAIN
		pamData['hydro_q'][...,3] = 0# SNOW
		pamData['hydro_q'][...,4] = 0# GRAUPEL


		# descriptorfile must be included. otherwise, pam.p.nhydro would be 0 which is not permitted. (OLD DESCRIPTOR FILE)
		descriptorFile = np.array([
			  #['hydro_name' 'as_ratio' 'liq_ice' 'rho_ms' 'a_ms' 'b_ms' 'alpha_as' 'beta_as' 'moment_in' 'nbin' 'dist_name' 'p_1' 'p_2' 'p_3' 'p_4' 'd_1' 'd_2' 'scat_name' 'vel_size_mod' 'canting']
			   ('cwc_q', -99.0, 1, -99.0, -99.0, -99.0, -99.0, -99.0, 3, 1, 'mono', -99.0, -99.0, -99.0, -99.0, 2e-05, -99.0, 'mie-sphere', 'khvorostyanov01_drops', -99.0),
			   ('iwc_q', -99.0, -1, -99.0, 130.0, 3.0, 0.684, 2.0, 3, 1, 'mono_cosmo_ice', -99.0, -99.0, -99.0, -99.0, -99.0, -99.0, 'mie-sphere', 'heymsfield10_particles', -99.0),
			   ('rwc_q', -99.0, 1, -99.0, -99.0, -99.0, -99.0, -99.0, 3, 50, 'exp', -99.0, -99.0, 8000000.0, -99.0, 0.00012, 0.006, 'mie-sphere', 'khvorostyanov01_drops', -99.0),
			   ('swc_q', -99.0, -1, -99.0, 0.038, 2.0, 0.3971, 1.88, 3, 50, 'exp_cosmo_snow', -99.0, -99.0, -99.0, -99.0, 5.1e-11, 0.02, 'mie-sphere', 'heymsfield10_particles', -99.0),
			   ('gwc_q', -99.0, -1, -99.0, 169.6, 3.1, -99.0, -99.0, 3, 50, 'exp', -99.0, -99.0, 4000000.0, -99.0, 1e-10, 0.01, 'mie-sphere', 'khvorostyanov01_spheres', -99.0)], 
			  dtype=[('hydro_name', 'S15'), ('as_ratio', '<f8'), ('liq_ice', '<i8'), ('rho_ms', '<f8'), ('a_ms', '<f8'), ('b_ms', '<f8'), ('alpha_as', '<f8'), ('beta_as', '<f8'), ('moment_in', '<i8'), ('nbin', '<i8'), ('dist_name', 'S15'), ('p_1', '<f8'), ('p_2', '<f8'), ('p_3', '<f8'), ('p_4', '<f8'), ('d_1', '<f8'), ('d_2', '<f8'), ('scat_name', 'S15'), ('vel_size_mod', 'S30'), ('canting', '<f8')]
			  )
		for hyd in descriptorFile: pam.df.addHydrometeor(hyd)


		# Create pamtra profile and go:
		pam.createProfile(**pamData)
		print("Starting PAMTRA on '" + filename_in + "':")
		# pam.runParallelPamtra(frq, pp_deltaX=0, pp_deltaY=0, pp_deltaF=1, pp_local_workers='auto')
		pam.runPamtra(frq)

		# save output:
		filename_out = pam_out_path + "pamtra_" + filename_in.replace(path_halo_dropsonde, '')
		pam.writeResultsToNetCDF(filename_out, xarrayCompatibleOutput=True, ncCompression=True)

		
		# # # # # Save the dropsonde launch number:
		# # # # sonde_number_filename = "/work/walbroel/data/" + "sonde_number_" + filename_in[-15:-3]
		# # # # np.save(sonde_number_filename, whichsonde)