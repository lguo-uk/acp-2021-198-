#!/usr/bin/env python
#SBATCH --job-name=var_trunk
#SBATCH --partition=short-serial-4hr
#SBATCH --time=02:00:00
#SBATCH --time-min=01:00:00

import cf
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("white")
import cfplot as cfp
import numpy as np
import os
import glob
import pickle
import cartopy.crs as ccrs
from scipy import stats

import sys
sys.path.insert(1,'/home/users/lguo/cssp_aero/analysis_script')
import local_functions3 as loc_func

'''
The IPCC colour scheme is:
	SSP1-1.9 = Hex #1E9684, RGB( 30, 150, 132); Green
	SSP1-2.6 = Hex #1D3354, RGB( 29,  51,  84); Blue
	SSP2-4.5 = Hex #EADD3D, RGB(234, 221,  61); Yellow
	SSP3-7.0 = Hex #F21111, RGB(242,  17,  17); Red
	SSP5-8.5 = Hex #840B22, RGB(132,  11,  34); Dark red
'''

'''
Sort lists of source_id and variant_label into a dictionary;
Ensemble members from the same source_id are sorted into a list under the same source_id;
'''
def sort_key_value(key_lst,value_lst):
	key_value_tuples=zip(key_lst,value_lst)
	key_value_dir={}
	for key,value in key_value_tuples:
		if key in key_value_dir:
			key_value_dir[key].append(value)
		else:
			key_value_dir[key]=[]
			key_value_dir[key].append(value)
	return key_value_dir

'''
Global variables
'''
scenarios=['ssp119','ssp126','ssp245','ssp370','ssp585']
time_slice=["2025-2034","2035-2044","2045-2054"]
obsdir='/home/users/lguo/cssp_aero/data/obs'
cmipdir='/home/users/lguo/cssp_aero/data/cmip6'
'''
Local Functions
'''
def filter_with_hwi_gt1(od550aer_mm_fld,experiment,period):
	'''
	read in hwi mm
	'''
	fi_hwi_mm="/home/users/lguo/cssp_aero/analysis_script/data/mm_hwi_{}.nc".format(experiment)
	print(fi_hwi_mm)
	hwi_cmip=cf.read(fi_hwi_mm).select_field("hwi")
	'''
	create a copy
	'''
	od550aer_mm_fld_gt1=od550aer_mm_fld.copy()
	'''
	loop through available model of each experiment
	'''
	for im,var in enumerate(od550aer_mm_fld):
		model_nm=var.aux('source_info').array.tolist()[0].decode("utf-8")
		print("Dealing with {}".format(model_nm))
		'''
		select corresponding model HWI
		'''
		try:
			hwi_mm=hwi_cmip.subspace(source_id=str.encode(model_nm)).squeeze()
			'''
			loop through HWI for filtering
			'''
			for it,ihwi in enumerate(hwi_mm.array):
				print(model_nm,it,ihwi)
				if ihwi<=1.:
					od550aer_mm_fld_gt1[im,it,:,:]=cf.masked
		except ValueError:
			'''
			mask out the mm, if there is not HWI index available
			'''
			print("{} does not have HWI data; mask it completely!".format(model_nm))
			od550aer_mm_fld_gt1[im,:,:,:]=cf.masked
	'''
	write out for future use
	'''
	cf.write(od550aer_mm_fld_gt1,"./data/od550aer_{}_{}_gt1.nc".format(experiment,period))

	return od550aer_mm_fld_gt1
'''
Calculate HWI from obs for 1979-2014
'''
if False:
	t850 = cf.read('/home/users/lguo/cssp_aero/analysis_script/data/ta_hwit850_obs.nc')[0]
	t250 = cf.read('/home/users/lguo/cssp_aero/analysis_script/data/ta_hwit250_obs.nc')[0]
	u500_1 = cf.read('/home/users/lguo/cssp_aero/analysis_script/data/ua_hwiu500_1_obs.nc')[0]
	u500_2 = cf.read('/home/users/lguo/cssp_aero/analysis_script/data/ua_hwiu500_2_obs.nc')[0]
	v850 = cf.read('/home/users/lguo/cssp_aero/analysis_script/data/va_hwiv850_obs.nc')[0]

	hwi,delta_t,delta_u,delta_v = loc_func.cal_hwi(t850,t850.collapse('T: mean'),t250,t250.collapse('T: mean'),u500_1,u500_1.collapse('T: mean'),u500_2,u500_2.collapse('T: mean'),v850,v850.collapse('T: mean'))

	cf.write(hwi,'/home/users/lguo/cssp_aero/analysis_script/data/hwi_obs_plus1.nc')
#	cf.write(delta_t,'/home/users/lguo/cssp_aero/analysis_script/data/hwi_delta_t_obs.nc')
#	cf.write(delta_u,'/home/users/lguo/cssp_aero/analysis_script/data/hwi_delta_u_obs.nc')
#	cf.write(delta_v,'/home/users/lguo/cssp_aero/analysis_script/data/hwi_delta_v_obs.nc')
'''
Calculate HWI for each matched model and each ensemble member in CMIP historical
'''
if False:
	mip='cmip'
	mip_var=['ta','ta','ua','ua','va']
	hwi_component=['hwit850','hwit250','hwiu500_1','hwiu500_2','hwiv850']
	'''
	Read in matched model name
	'''
	with open("./data/hwi_match_cmip-scenario.pickle","rb") as fp:
		model_lst=pickle.load(fp)

	hwi_component_dir_model_dir={}
	for i,component in enumerate(hwi_component):
		print("===== "+component+" =====")
		var_fldlst=cf.read('/home/users/lguo/cssp_aero/data/cmip6/{0}/historical/{1}/{2}/{1}*_185001-201412.nc'.format(mip,mip_var[i],component))
		hwi_component_dir_model_dir[component]={}
		for model in model_lst:
			print(model)
			hwi_component_dir_model_dir[component][model]=var_fldlst.select('source_id={}'.format(model))
	'''
	Calculate HWI and concatenate all ensemble members (according to v850, the shortest list)
	'''
	hwi_lst=[];hwi_delta_t_lst=[];hwi_delta_u_lst=[];hwi_delta_v_lst=[]
	for i,model in enumerate(model_lst):
		for j,sample in enumerate(hwi_component_dir_model_dir["hwiv850"][model]):
			print(model,sample.get_property("variant_label"))
			t850=hwi_component_dir_model_dir["hwit850"][model].select_field("variant_label="+sample.get_property("variant_label"))[1548:].subspace(T=cf.djf())
			t250=hwi_component_dir_model_dir["hwit250"][model].select_field("variant_label="+sample.get_property("variant_label"))[1548:].subspace(T=cf.djf())
			u500_1=hwi_component_dir_model_dir["hwiu500_1"][model].select_field("variant_label="+sample.get_property("variant_label"))[1548:].subspace(T=cf.djf())
			u500_2=hwi_component_dir_model_dir["hwiu500_2"][model].select_field("variant_label="+sample.get_property("variant_label"))[1548:].subspace(T=cf.djf())
			v850=hwi_component_dir_model_dir["hwiv850"][model].select_field("variant_label="+sample.get_property("variant_label"))[1548:].subspace(T=cf.djf())
			hwi,delta_t,delta_u,delta_v=loc_func.cal_hwi(t850,t850.collapse('T: mean'),t250,t250.collapse('T: mean'),u500_1,u500_1.collapse('T: mean'),u500_2,u500_2.collapse('T: mean'),v850,v850.collapse('T: mean'))
			'''
			Save for AOD comparison
			'''
			dirout="/home/users/lguo/cssp_aero/data/cmip6/cmip/historical/hwi"
			if not os.path.isdir(dirout): os.makedirs(dirout)
			fonm='{}/hwi_{}_historical_{}_gn_185001-201412.nc'.format(dirout,model,hwi.get_property("variant_label"))
			print(fonm)
			cf.write(hwi,fonm)
			#
			hwi_lst.append(hwi);hwi_delta_t_lst.append(delta_t);hwi_delta_u_lst.append(delta_u);hwi_delta_v_lst.append(delta_v)
	'''
	Calculate model mean
	'''
	hwi_lst_mm=loc_func.fldlst2lstmm(hwi_lst);hwi_delta_t_lst_mm=loc_func.fldlst2lstmm(hwi_delta_t_lst);hwi_delta_u_lst_mm=loc_func.fldlst2lstmm(hwi_delta_u_lst);hwi_delta_v_lst_mm=loc_func.fldlst2lstmm(hwi_delta_v_lst)
	'''
	Write out hwi_lst_mm for fig:idx_scenario_anomaly
	Note, cf.aggregate does not work in cf-python3.
	And, writing a FieldList to a single .nc file stops wporking.
	Therefore, mm are wrote to individual files and are concatenated using necat in the data directory.
	'''
	for i,mdlnm in enumerate(model_lst):
		hwi_lst_mm_normalised=(hwi_lst_mm[i]-hwi_lst_mm[i].collapse('T: mean'))/hwi_lst_mm[i].collapse('T: sd')
		cf.write(hwi_lst_mm_normalised,'/home/users/lguo/cssp_aero/analysis_script/data/hwi_historical_{:02d}.nc'.format(i))
		hwi_delta_t_lst_mm_normalised=(hwi_delta_t_lst_mm[i]-hwi_delta_t_lst_mm[i].collapse('T: mean'))/hwi_delta_t_lst_mm[i].collapse('T: sd')
		cf.write(hwi_delta_t_lst_mm_normalised,'/home/users/lguo/cssp_aero/analysis_script/data/hwi_delta_t_historical_{:02d}.nc'.format(i))
		hwi_delta_u_lst_mm_normalised=(hwi_delta_u_lst_mm[i]-hwi_delta_u_lst_mm[i].collapse('T: mean'))/hwi_delta_u_lst_mm[i].collapse('T: sd')
		cf.write(hwi_delta_u_lst_mm_normalised,'/home/users/lguo/cssp_aero/analysis_script/data/hwi_delta_u_historical_{:02d}.nc'.format(i))
		hwi_delta_v_lst_mm_normalised=(hwi_delta_v_lst_mm[i]-hwi_delta_v_lst_mm[i].collapse('T: mean'))/hwi_delta_v_lst_mm[i].collapse('T: sd')
		cf.write(hwi_delta_v_lst_mm_normalised,'/home/users/lguo/cssp_aero/analysis_script/data/hwi_delta_v_historical_{:02d}.nc'.format(i))
		if i==0:
			hwi_conc=hwi_lst_mm_normalised.array
			hwi_delta_t_conc=hwi_delta_t_lst_mm_normalised.array
			hwi_delta_u_conc=hwi_delta_u_lst_mm_normalised.array
			hwi_delta_v_conc=hwi_delta_v_lst_mm_normalised.array
		else:
			hwi_conc=np.concatenate((hwi_conc,hwi_lst_mm_normalised.array))
			hwi_delta_t_conc=np.concatenate((hwi_delta_t_conc,hwi_delta_t_lst_mm_normalised.array))
			hwi_delta_u_conc=np.concatenate((hwi_delta_u_conc,hwi_delta_u_lst_mm_normalised.array))
			hwi_delta_v_conc=np.concatenate((hwi_delta_v_conc,hwi_delta_v_lst_mm_normalised.array))
	'''
	Write out hwi_conc for ploting
	'''
	with open('/home/users/lguo/cssp_aero/analysis_script/data/hwi_historical_concatenated.npy','wb') as f:
		np.save(f,hwi_conc)
	with open('/home/users/lguo/cssp_aero/analysis_script/data/hwi_delta_t_historical_concatenated.npy','wb') as f:
		np.save(f,hwi_delta_t_conc)
	with open('/home/users/lguo/cssp_aero/analysis_script/data/hwi_delta_u_historical_concatenated.npy','wb') as f:
		np.save(f,hwi_delta_u_conc)
	with open('/home/users/lguo/cssp_aero/analysis_script/data/hwi_delta_v_historical_concatenated.npy','wb') as f:
		np.save(f,hwi_delta_v_conc)
'''
Calculate HWI for local data historical
'''
if False:
	mip='cmip'
	mip_var=['ta','ta','ua','ua','va']
	hwi_component=['hwit850','hwit250','hwiu500_1','hwiu500_2','hwiv850']
	'''
	Read in matched model name
	'''
	model_lst=["MIROC-ES2L"]

	hwi_component_dir_model_dir={}
	for i,component in enumerate(hwi_component):
		print("===== "+component+" =====")
		print('/home/users/lguo/cssp_aero/data/cmip6/{0}/historical/{1}/{2}/{1}*_185001-201412.nc'.format(mip,mip_var[i],component))
		var_fldlst=cf.read('/home/users/lguo/cssp_aero/data/cmip6/{0}/historical/{1}/{2}/{1}*_185001-201412.nc'.format(mip,mip_var[i],component))
		hwi_component_dir_model_dir[component]={}
		for model in model_lst:
			print(model)
			hwi_component_dir_model_dir[component][model]=var_fldlst.select('source_id={}'.format(model))
	'''
	Calculate HWI and concatenate all ensemble members (according to v850, the shortest list)
	'''
	hwi_lst=[];hwi_delta_t_lst=[];hwi_delta_u_lst=[];hwi_delta_v_lst=[]
	for i,model in enumerate(model_lst):
		for j,sample in enumerate(hwi_component_dir_model_dir["hwiv850"][model]):
			print(model,sample.get_property("variant_label"))
			t850=hwi_component_dir_model_dir["hwit850"][model].select_field("variant_label="+sample.get_property("variant_label"))[1548:].subspace(T=cf.djf())
			t250=hwi_component_dir_model_dir["hwit250"][model].select_field("variant_label="+sample.get_property("variant_label"))[1548:].subspace(T=cf.djf())
			u500_1=hwi_component_dir_model_dir["hwiu500_1"][model].select_field("variant_label="+sample.get_property("variant_label"))[1548:].subspace(T=cf.djf())
			u500_2=hwi_component_dir_model_dir["hwiu500_2"][model].select_field("variant_label="+sample.get_property("variant_label"))[1548:].subspace(T=cf.djf())
			v850=hwi_component_dir_model_dir["hwiv850"][model].select_field("variant_label="+sample.get_property("variant_label"))[1548:].subspace(T=cf.djf())
			hwi,delta_t,delta_u,delta_v=loc_func.cal_hwi(t850,t850.collapse('T: mean'),t250,t250.collapse('T: mean'),u500_1,u500_1.collapse('T: mean'),u500_2,u500_2.collapse('T: mean'),v850,v850.collapse('T: mean'))
			'''
			Save for AOD comparison
			'''
			dirout="/home/users/lguo/cssp_aero/data/cmip6/cmip/historical/hwi"
			if not os.path.isdir(dirout): os.makedirs(dirout)
			fonm='{}/hwi_{}_historical_{}_gn_185001-201412.nc'.format(dirout,model,hwi.get_property("variant_label"))
			print(fonm)
			cf.write(hwi,fonm)
'''
Calculate HWI for each matched model in Historical (without normalisation)
'''
if True:
	mip='CMIP'
	mip_var=['ta','ta','ua','ua','va']
	hwi_component=['hwit850','hwit250','hwiu500_1','hwiu500_2','hwiv850']
	'''
	Read in matched model name
	'''
	with open("./data/hwi_match_cmip-scenario.pickle","rb") as fp:
		model_lst=pickle.load(fp)

	hwi_component_dir_model_dir={}
	for i,component in enumerate(hwi_component):
		print("===== "+component+" =====")
		var_fldlst=cf.read('/home/users/lguo/cssp_aero/data/cmip6/{0}/historical/{1}/{2}/{1}*_185001-201412.nc'.format(mip,mip_var[i],component))
		hwi_component_dir_model_dir[component]={}
		for model in model_lst:
			print(model)
			hwi_component_dir_model_dir[component][model]=var_fldlst.select('source_id={}'.format(model))
	'''
	Calculate HWI and concatenate all ensemble members (according to v850, the shortest list)
	'''
	hwi_lst=[];hwi_delta_t_lst=[];hwi_delta_u_lst=[];hwi_delta_v_lst=[]
	for i,model in enumerate(model_lst):
		for j,sample in enumerate(hwi_component_dir_model_dir["hwiv850"][model]):
			print(model,sample.get_property("variant_label"))
			t850=hwi_component_dir_model_dir["hwit850"][model].select_field("variant_label="+sample.get_property("variant_label"))[1548:].subspace(T=cf.djf())
			t250=hwi_component_dir_model_dir["hwit250"][model].select_field("variant_label="+sample.get_property("variant_label"))[1548:].subspace(T=cf.djf())
			u500_1=hwi_component_dir_model_dir["hwiu500_1"][model].select_field("variant_label="+sample.get_property("variant_label"))[1548:].subspace(T=cf.djf())
			u500_2=hwi_component_dir_model_dir["hwiu500_2"][model].select_field("variant_label="+sample.get_property("variant_label"))[1548:].subspace(T=cf.djf())
			v850=hwi_component_dir_model_dir["hwiv850"][model].select_field("variant_label="+sample.get_property("variant_label"))[1548:].subspace(T=cf.djf())
			hwi,delta_t,delta_u,delta_v=loc_func.cal_hwi(t850,t850.collapse('T: mean'),t250,t250.collapse('T: mean'),u500_1,u500_1.collapse('T: mean'),u500_2,u500_2.collapse('T: mean'),v850,v850.collapse('T: mean'))
			'''
			Save for AOD comparison
			'''
			dirout="/home/users/lguo/cssp_aero/data/cmip6/CMIP/historical/hwi_raw"
			if not os.path.isdir(dirout): os.makedirs(dirout)
			fonm='{}/hwi_{}_historical_{}_gn_185001-201412.nc'.format(dirout,model,hwi.get_property("variant_label"))
			print(fonm)
			cf.write(hwi,fonm)
			#
			hwi_lst.append(hwi);hwi_delta_t_lst.append(delta_t);hwi_delta_u_lst.append(delta_u);hwi_delta_v_lst.append(delta_v)

			breakpoint()
'''
Plot HWI histogram for obs and CMIP historical concatenation
'''
if False:
	title=["HWI",r"HWI $\Delta T$",r"HWI $\Delta U$",r"HWI $\Delta V$"]
#	for i,component in enumerate(["hwi","hwi_delta_t","hwi_delta_u","hwi_delta_v"]):
	for i,component in enumerate(["hwi"]):
		hwi_obs=cf.read('/home/users/lguo/cssp_aero/analysis_script/data/{}_obs.nc'.format(component))[0]
		with open('/home/users/lguo/cssp_aero/analysis_script/data/{}_historical_concatenated.npy'.format(component),'rb') as f:
			hwi_cmip=np.load(f)
		'''
		Kolmogorov-Smirnov test: H0: Two samples are drawn from the same distribution.
		'''
		d,pvalue=stats.ks_2samp(hwi_obs.array,hwi_cmip)
		print("The P-value of Kolmogorov-Smirnov test is:",pvalue)

		kwargs=dict(bins=50,hist_kws={'alpha':.6},kde_kws={'linewidth':2})
		sns.distplot(hwi_obs,color="dodgerblue",label="ERA5",**kwargs)
		sns.distplot(hwi_cmip,color="orange",label="CMIP6",**kwargs)
		plt.ylim(0,1.)
		plt.xlim(-5,5)
		plt.xlabel(title[i],fontsize=15)
		plt.ylabel("Probability Density Function",fontsize=15)
		plt.gca().tick_params(labelsize="15",direction="out",bottom=True,left=True)
		plt.legend(fontsize=15)
		plt.savefig('./plot/hwi_scenariomip/fig_{}_historical_pdf_matched.pdf'.format(component),bbox_inches='tight',pad_inches=0)
		plt.show()
		plt.close()
'''
Plot HWI histogram for individual model
'''
if False:
	hwi_obs=cf.read('/home/users/lguo/cssp_aero/analysis_script/data/hwi_obs.nc')[0]
	for mno in range(17):
		hwi_cmip=cf.read('/home/users/lguo/cssp_aero/analysis_script/data/hwi_historical_{:02d}.nc'.format(mno))[0]

		kwargs=dict(bins=50,hist_kws={'alpha':.6},kde_kws={'linewidth':2})
		sns.distplot(hwi_obs,color="dodgerblue",label="ERA5",**kwargs)
		sns.distplot(hwi_cmip,color="orange",label="model",**kwargs)
		plt.ylim(0,1.)
		plt.xlim(-5,5)
		plt.gca().set(title='PDF of HWI in DJF (1979-2014)',xlabel=hwi_cmip.get_property("source_id"),ylabel='Density of Probability')
		plt.legend()
		plt.show()
#		plt.savefig(fig_nm,bbox_inches='tight',pad_inches=0)
		plt.close()
'''
Calculate HWI:
	for (not just matched) models and their ensemble members in ScenarioMIP
'''
if False:
	mip='ScenarioMIP'
	scenarios=['ssp126','ssp245','ssp370','ssp585']
	mip_var=['ta','ta','ua','ua','va']
	hwi_component=['hwit850','hwit250','hwiu500_1','hwiu500_2','hwiv850']
	mip_var_standardname=[]
	'''
	Read Historical
	'''
	historical_hwi_component_dict={}
	for j,component in enumerate(hwi_component):
		print("historical",component)
		os.system("ls /home/users/lguo/cssp_aero/data/cmip6/cmip/historical/{0}/{1}/{0}_Amon_*_historical_*_gn_185001-201412.nc".format(mip_var[j],component))
		historical_hwi_component_dict[component]=cf.read("/home/users/lguo/cssp_aero/data/cmip6/cmip/historical/{0}/{1}/{0}_Amon_*_historical_*_gn_185001-201412.nc".format(mip_var[j],component))
	'''
	Read Scenario
	'''
	scenario_dict_hwi_component_dict={}
	for i,scenario in enumerate(scenarios):
		scenario_dict_hwi_component_dict[scenario]={}
		for j,component in enumerate(hwi_component):
			print(scenario,component)
			os.system("ls /home/users/lguo/cssp_aero/data/cmip6/{0}/{1}/{2}/{3}/{2}_Amon_*_{1}_*_gn_201501-210012.nc".format(mip,scenario,mip_var[j],component))
			scenario_dict_hwi_component_dict[scenario][component]=cf.read("/home/users/lguo/cssp_aero/data/cmip6/{0}/{1}/{2}/{3}/{2}_Amon_*_{1}_*_gn_201501-210012.nc".format(mip,scenario,mip_var[j],component))
	'''
	Check availability among components. Choose the shortest component as baseline
	'''
	for scenario in scenarios:
		print(scenario,[len(scenario_dict_hwi_component_dict[scenario][component]) for component in hwi_component])
	base_component=["hwit850","hwiv850","hwiv850","hwit850","hwiv850"]
	'''
	Calculate HWI and concatenate all ensemble members
	'''
	for i,scenario in enumerate(scenarios):
		hwi_scenario_lst=[];hwi_delta_t_scenario_lst=[];hwi_delta_u_scenario_lst=[];hwi_delta_v_scenario_lst=[]
		for sample in scenario_dict_hwi_component_dict[scenario][base_component[i]]:
			try:
				print(scenario,sample.get_property("source_id"),sample.get_property("variant_label"))
				t850=scenario_dict_hwi_component_dict[scenario]["hwit850"].select("source_id="+sample.get_property("source_id")).select_field("variant_label="+sample.get_property("variant_label")).subspace(T=cf.djf())
				t250=scenario_dict_hwi_component_dict[scenario]["hwit250"].select("source_id="+sample.get_property("source_id")).select_field("variant_label="+sample.get_property("variant_label")).subspace(T=cf.djf())
				u500_1=scenario_dict_hwi_component_dict[scenario]["hwiu500_1"].select("source_id="+sample.get_property("source_id")).select_field("variant_label="+sample.get_property("variant_label")).subspace(T=cf.djf())
				u500_2=scenario_dict_hwi_component_dict[scenario]["hwiu500_2"].select("source_id="+sample.get_property("source_id")).select_field("variant_label="+sample.get_property("variant_label")).subspace(T=cf.djf())
				v850=scenario_dict_hwi_component_dict[scenario]["hwiv850"].select("source_id="+sample.get_property("source_id")).select_field("variant_label="+sample.get_property("variant_label")).subspace(T=cf.djf())

				print("historical",sample.get_property("source_id"),sample.get_property("variant_label"))
				t850_historical=historical_hwi_component_dict["hwit850"].select("source_id="+sample.get_property("source_id")).select_field("variant_label="+sample.get_property("variant_label")).subspace(T=cf.wi(cf.dt("1979-01-01"),cf.dt("2014-12-30"))).subspace(T=cf.djf())
				t250_historical=historical_hwi_component_dict["hwit250"].select("source_id="+sample.get_property("source_id")).select_field("variant_label="+sample.get_property("variant_label")).subspace(T=cf.wi(cf.dt("1979-01-01"),cf.dt("2014-12-30"))).subspace(T=cf.djf())
				u500_1_historical=historical_hwi_component_dict["hwiu500_1"].select("source_id="+sample.get_property("source_id")).select_field("variant_label="+sample.get_property("variant_label")).subspace(T=cf.wi(cf.dt("1979-01-01"),cf.dt("2014-12-30"))).subspace(T=cf.djf())
				u500_2_historical=historical_hwi_component_dict["hwiu500_2"].select("source_id="+sample.get_property("source_id")).select_field("variant_label="+sample.get_property("variant_label")).subspace(T=cf.wi(cf.dt("1979-01-01"),cf.dt("2014-12-30"))).subspace(T=cf.djf())
				v850_historical=historical_hwi_component_dict["hwiv850"].select("source_id="+sample.get_property("source_id")).select_field("variant_label="+sample.get_property("variant_label")).subspace(T=cf.wi(cf.dt("1979-01-01"),cf.dt("2014-12-30"))).subspace(T=cf.djf())


				hwi,delta_t,delta_u,delta_v=loc_func.cal_hwi_scenario(t850,t850_historical,t250,t250_historical,u500_1,u500_1_historical,u500_2,u500_2_historical,v850,v850_historical)
				hwi_scenario_lst.append(hwi);hwi_delta_t_scenario_lst.append(delta_t);hwi_delta_u_scenario_lst.append(delta_u);hwi_delta_v_scenario_lst.append(delta_v)
			except ValueError:
				print("... One component in either Scenario or Historical is not available .....")
				continue
			'''
			Save for AOD comparison
			'''
			dirout="/home/users/lguo/cssp_aero/data/cmip6/ScenarioMIP/{}/hwi".format(scenario)
			if not os.path.isdir(dirout): os.makedirs(dirout)

			fonm='{}/hwi_{}_{}_{}_gn_197901-210012.nc'.format(dirout,t850.get_property("source_id"),scenario,t850.get_property("variant_label"))
			print(fonm)
			cf.write(hwi,fonm)

		hwi_scenario_lst_mm=loc_func.fldlst2lstmm(hwi_scenario_lst);hwi_delta_t_scenario_lst_mm=loc_func.fldlst2lstmm(hwi_delta_t_scenario_lst);hwi_delta_u_scenario_lst_mm=loc_func.fldlst2lstmm(hwi_delta_u_scenario_lst);hwi_delta_v_scenario_lst_mm=loc_func.fldlst2lstmm(hwi_delta_v_scenario_lst)

		for im,mm in enumerate(hwi_scenario_lst_mm):
			mm_normalised=(mm-mm.collapse('T: mean'))/mm.collapse('T: sd')
			cf.write(mm_normalised,'/home/users/lguo/cssp_aero/analysis_script/data/hwi_historical-{}_{:02d}.nc'.format(scenario,im))
			'''
			To put model name into the concatenated .nc file. model name needs to be record
			'''
			if im == 0:
				f=open("/home/users/lguo/cssp_aero/analysis_script/data/hwi_historical-{}_mdlnm.txt".format(scenario),"w")
				f.write('"{}"s'.format(mm.get_property("source_id")))
			else:
				f=open("/home/users/lguo/cssp_aero/analysis_script/data/hwi_historical-{}_mdlnm.txt".format(scenario),"a")
				f.write(',"{}"s'.format(mm.get_property("source_id")))
		for im,mm in enumerate(hwi_delta_t_scenario_lst_mm):
			mm_normalised=(mm-mm.collapse('T: mean'))/mm.collapse('T: sd')
			cf.write(mm_normalised,'/home/users/lguo/cssp_aero/analysis_script/data/hwi_delta_t_historical-{}_{:02d}.nc'.format(scenario,im))
		for im,mm in enumerate(hwi_delta_u_scenario_lst_mm):
			mm_normalised=(mm-mm.collapse('T: mean'))/mm.collapse('T: sd')
			cf.write(mm_normalised,'/home/users/lguo/cssp_aero/analysis_script/data/hwi_delta_u_historical-{}_{:02d}.nc'.format(scenario,im))
		for im,mm in enumerate(hwi_delta_v_scenario_lst_mm):
			mm_normalised=(mm-mm.collapse('T: mean'))/mm.collapse('T: sd')
			cf.write(mm_normalised,'/home/users/lguo/cssp_aero/analysis_script/data/hwi_delta_v_historical-{}_{:02d}.nc'.format(scenario,im))
'''
Visualise HWI changes in ScenarioMIP
'''
if False:
	scenarios=scenarios[1:]
	var_titles=["HWI anomaly","HWI $\Delta T$ anomaly","HWI $\Delta U$ anomaly","HWI $\Delta V$ anomaly"]

	for ivar,var in enumerate(["hwi","hwi_delta_t","hwi_delta_u","hwi_delta_v"]):
		hwi_all={}
		hwi_scenario_ano={};hwi_scenario_ano_data={} # it is a dictionary of scenarios; each entry contains a 2D list with 0d as time slice and 1d as mm of HWIano.
		hwi_scenario={};hwi_scenario_data={}
		hwi_historical={};hwi_historical_data={}
		for scenario in scenarios:
			hwi_all[scenario]=cf.read('/home/users/lguo/cssp_aero/analysis_script/data/mm_{}_historical-{}.nc'.format(var,scenario))[0]

			hwi_scenario_ano[scenario],hwi_scenario[scenario],hwi_historical[scenario]=loc_func.hwi_anomaly3(hwi_all[scenario])

			hwi_scenario_ano_data[scenario]=loc_func.fldlst2list_2level(hwi_scenario_ano[scenario])
			hwi_scenario_data[scenario]=loc_func.fldlst2list_2level(hwi_scenario[scenario])
			hwi_historical_data[scenario]=loc_func.fldlst2list_2level(hwi_historical[scenario])
		'''
		Make boxplot
		'''
		boxcolors=['#1E9684','#1D3354','#EADD3D','#F21111','#840B22']
		boxcolors=boxcolors[1:]
		medianprops=dict(linewidth=2,color='k')
		box={}
		for i,scenario in enumerate(scenarios):
			boxprops=dict(color=boxcolors[i],facecolor="w",linewidth=2,alpha=.7)
			pos=[-.5+i*.25, 1.5+i*.25, 3.5+i*.25, 5.5+i*.25]
			box[scenario]=plt.boxplot(hwi_scenario_ano_data[scenario],whis=(5,95),widths=.15,positions=pos,showcaps=True,showfliers=False,patch_artist=True,medianprops=medianprops,boxprops=boxprops)
			'''
			Conduct ttest between future and historical
			'''
			for ip,po in enumerate(pos):
				sig=loc_func.ttest4mean(np.nanmean(hwi_scenario_data[scenario][ip]),np.nanmean(hwi_historical_data[scenario][ip]),np.nanstd(hwi_scenario_data[scenario][ip]),np.nanstd(hwi_historical_data[scenario][ip]),len(hwi_scenario_data[scenario][ip]),len(hwi_historical_data[scenario][ip]))
				print(scenario,ip,sig)
				if sig<.1:
					print("plotting ...")
					boxprops=dict(color=boxcolors[i],facecolor=boxcolors[i],linewidth=2,alpha=.7)
					plt.boxplot(hwi_scenario_ano_data[scenario][ip],whis=(5,95),widths=.15,positions=[pos[ip]],showcaps=True,showfliers=False,patch_artist=True,medianprops=medianprops,boxprops=boxprops)
		'''
		Add individual model
		'''
		markers=[".",",","o","v","^","<",">","1","2","3","4","8","s","p","P","*","h","H","+","x","X","D","d"]
		model_lst=[]
		model_lgd=[]
		model_mkr={};imkr=0
		model_clr={};iclr=0
		for i,snr in enumerate(scenarios):
			pos=[-.5+i*.25, 1.5+i*.25, 3.5+i*.25, 5.5+i*.25]
			for j,hwia_period in enumerate(hwi_scenario_ano[snr]):
				for hwia in hwia_period:
					mdlnm=hwia.coordinates("source_id").value().array[0].decode("utf-8")
					print(snr,mdlnm,hwia.array[0])
				#	if mdlnm in model_lst:
				#		plt.scatter(pos[j],hwia.array[0],s=15,c=model_clr[mdlnm],marker=model_mkr[mdlnm])
				#	else:
				#		model_lst.append(mdlnm)
				#		model_mkr[mdlnm]=markers[imkr];imkr+=1
				#		'''
				#		Generate RBG color randomly
				#		'''
				#		mkr_rbg=np.random.choice(range(256),size=3)
				#		'''
				#		Convert RBG to hex color, as RBG sequence is indistinguishable from an array of values to be colormapped
				#		'''
				#		model_clr[mdlnm]='#{:02x}{:02x}{:02x}'.format(mkr_rbg[0],mkr_rbg[1],mkr_rbg[2])
				#		model_lgd.append(plt.scatter(pos[j],hwia.array[0],s=15,c=model_clr[mdlnm],marker=model_mkr[mdlnm]))
					'''
					Read in recorded colours
					'''
					with open('/home/users/lguo/cssp_aero/analysis_script/data/marker_colour.pickle','rb') as f:
						model_clr=pickle.load(f)
	
					if mdlnm in model_lst:
						plt.scatter(pos[j],hwia.array[0],s=15,c=model_clr[mdlnm],marker=model_mkr[mdlnm])
					else:
						try:
							colour_key=model_clr[mdlnm]
						except KeyError:
							print("{} is not in Model List ...".format(mdlnm))
							continue
						model_lst.append(mdlnm)
						model_mkr[mdlnm]=markers[imkr];imkr+=1
						model_lgd.append(plt.scatter(pos[j],hwia.array[0],s=15,c=model_clr[mdlnm],marker=model_mkr[mdlnm]))
		'''
		Make double legends:
			1. build legend0;
			2. plot legend1;
			3. plot legend0;
		'''
		scenarios_str=['SSP 1-1.9','SSP 1-2.6','SSP 2-4.5','SSP 3-7.0','SSP 5-8.5']
		scenarios_str=scenarios_str[1:]
		legend0=plt.legend([box[snr]["boxes"][0] for snr in scenarios],scenarios_str,bbox_to_anchor=(.02,.81,.05,.12),loc='upper left',fontsize=15)
		plt.legend(model_lgd, model_lst, bbox_to_anchor=(0,.95,1,0.3), loc="lower left", mode="expand", ncol=4, fontsize='xx-small')
		if ivar==0:
			plt.gca().add_artist(legend0)
		'''
		x&y axes
		'''
		plt.xticks([0,2,4,6],['2025-2034','2035-2044','2045-2054','2090-2099'])
		plt.ylim(-1,2.5)
		plt.xlim(-1,7)
		plt.plot([-1,7],[0,0],c='k',lw=1,ls=':')
		plt.ylabel(var_titles[ivar],fontsize=20)
		plt.gca().tick_params(labelsize="15",direction="in",bottom=True,left=True)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().spines['top'].set_visible(False)
		'''
		output
		'''
		plt.savefig('./plot/hwi_scenariomip/fig_{}_scenario_anomaly.pdf'.format(var),bbox_inches='tight',pad_inches=0)
	#	plt.show()
		plt.close()
'''
Visualise HWI changes in ScenarioMIP
	Single model
'''
if False:
	scenarios=scenarios[1:]

	model="CanESM5"
	var="hwi"

	hwi_all={};mm_all={};
	hwi_scenario_ano={};hwi_scenario_ano_data={};mm_scenario_ano={};mm_scenario_ano_data={};
	hwi_scenario={};hwi_scenario_data={};mm_scenario={};mm_scenario_data={};
	hwi_historical={};hwi_historical_data={};mm_historical={};mm_historical_data={};
	for scenario in scenarios:
		hwi_all[scenario]=loc_func.aggregate_fieldlist(cf.read("{0}/ScenarioMIP/{1}/hwi/hwi_{2}_{1}_*_gn_197901-210012.nc".format(cmipdir,scenario,model)))

		hwi_scenario_ano[scenario],hwi_scenario[scenario],hwi_historical[scenario],=loc_func.hwi_anomaly3(hwi_all[scenario])

		hwi_scenario_ano_data[scenario]=loc_func.fldlst2list_2level(hwi_scenario_ano[scenario])
		hwi_scenario_data[scenario]=loc_func.fldlst2list_2level(hwi_scenario[scenario])
		hwi_historical_data[scenario]=loc_func.fldlst2list_2level(hwi_historical[scenario])

		mm_all[scenario]=cf.read('/home/users/lguo/cssp_aero/analysis_script/data/mm_{}_historical-{}.nc'.format(var,scenario))[0]

		mm_scenario_ano[scenario],mm_scenario[scenario],mm_historical[scenario],=loc_func.hwi_anomaly3(mm_all[scenario])

		mm_scenario_ano_data[scenario]=loc_func.fldlst2list_2level(mm_scenario_ano[scenario])
		mm_scenario_data[scenario]=loc_func.fldlst2list_2level(mm_scenario[scenario])
		mm_historical_data[scenario]=loc_func.fldlst2list_2level(mm_historical[scenario])
	'''
	Make boxplot
	'''
	boxcolors=['#1E9684','#1D3354','#EADD3D','#F21111','#840B22']
	boxcolors=boxcolors[1:]
	medianprops=dict(linewidth=2,color='k')
	box={}
	for i,scenario in enumerate(scenarios):
		boxprops=dict(color=boxcolors[i],facecolor="w",linewidth=2,alpha=.7)
		pos=[-.5+i*.25, 1.5+i*.25, 3.5+i*.25, 5.5+i*.25]
		box[scenario]=plt.boxplot(hwi_scenario_ano_data[scenario],positions=pos,whis=(5,95),widths=.15,showcaps=True,showfliers=False,patch_artist=True,medianprops=medianprops,boxprops=boxprops)
		'''
		Conduct ftest
		'''
		for ip,po in enumerate(pos):
			sig=loc_func.ftest(np.nanstd(mm_scenario_data[scenario][ip]),np.nanstd(hwi_scenario_data[scenario][ip]),len(mm_scenario_data[scenario][ip]),len(hwi_scenario_data[scenario][ip]))
			print(scenario,ip,sig)
			if sig<.1:
				print("plotting ...")
				boxprops=dict(color=boxcolors[i],facecolor=boxcolors[i],linewidth=2,alpha=.7)
				plt.boxplot(hwi_scenario_ano_data[scenario][ip],whis=(5,95),positions=[pos[ip]],widths=.15,showcaps=True,showfliers=False,patch_artist=True,medianprops=medianprops,boxprops=boxprops)
	'''
	Add individual model
	'''
	for i,scenario in enumerate(scenarios):
		pos=[-.5+i*.25, 1.5+i*.25, 3.5+i*.25, 5.5+i*.25]
		for j,hwia_period in enumerate(hwi_scenario_ano[scenario]):
			for hwia in hwia_period:
				plt.scatter(pos[j],hwia.array[0],s=15,c="k",marker=".",alpha=.6)
	'''
	Legend
	'''
	scenarios_str=['SSP 1-1.9','SSP 1-2.6','SSP 2-4.5','SSP 3-7.0','SSP 5-8.5']
	scenarios_str=scenarios_str[1:]
	plt.legend([box[scenario]["boxes"][0] for scenario in scenarios],scenarios_str,title=model,bbox_to_anchor=(.02,.89,.05,.12),loc='upper left')
	'''
	x&y axes
	'''
	plt.xticks([0,2,4,6],['2025-2034','2035-2044','2045-2054','2090-2099'])
	plt.ylim(-1,2.5)
	plt.xlim(-1,7)
	plt.gca().spines['right'].set_visible(False)
	plt.gca().spines['top'].set_visible(False)
	plt.gca().yaxis.set_ticks_position('left')
	plt.gca().xaxis.set_ticks_position('bottom')
	plt.gca().set(ylabel='{} Anomaly'.format(var.upper()))
	'''
	output
	'''
	plt.savefig('./plot/hwi_scenariomip/fig_{}_scenario_anomaly_{}.pdf'.format(var,model),bbox_inches='tight',pad_inches=0)
	plt.show()
	plt.close()
'''
Calculate Wang-Chen Index from obs for 1979-2014 DJF
'''
if False:
	sib_obs=cf.read('/home/users/lguo/cssp_aero/data/obs/psl_wci_sib_djf19792014.nc').select_field("air_pressure_at_mean_sea_level")
	al_obs=cf.read('/home/users/lguo/cssp_aero/data/obs/psl_wci_al_djf19792014.nc').select_field("air_pressure_at_mean_sea_level")
	mar_obs=cf.read('/home/users/lguo/cssp_aero/data/obs/psl_wci_mar_djf19792014.nc').select_field("air_pressure_at_mean_sea_level")
	wci_obs=loc_func.cal_wci(sib_obs,al_obs,mar_obs)
	cf.write(wci_obs,'/home/users/lguo/cssp_aero/analysis_script/data/wci_obs.nc')
'''
Calculate WCI for each matched model and each ensemble member in CMIP historical
'''
if False:
	mip='cmip'
	mip_var="psl"
	wci_component=["sib","al","mar"]
	'''
	Read in matched model name
	'''
	with open("./data/hwi_match_cmip-scenario.pickle","rb") as fp:
		model_lst=pickle.load(fp)
	'''
	Create a 2-level dictionary with level 1 HWI_component and level 2 model_name
	'''
	wci_component_dir_model_dir={}
	for i,component in enumerate(wci_component):
		print("===== "+component+" =====")
		var_fldlst=cf.read('/home/users/lguo/cssp_aero/data/cmip6/{0}/historical/{1}/{2}/{1}*_185001-201412.nc'.format(mip,mip_var,component))
		wci_component_dir_model_dir[component]={}
		for model in model_lst:
			print(model)
			wci_component_dir_model_dir[component][model]=var_fldlst.select('source_id={}'.format(model))
	'''
	Calculate HWI and concatenate all ensemble members
	'''
	wci_lst=[];wci_sib_lst=[];wci_al_lst=[];wci_mar_lst=[]
	for i,model in enumerate(model_lst):
		for j,sample in enumerate(wci_component_dir_model_dir["sib"][model]):
			print(model,sample.get_property("variant_label"))
			sib=wci_component_dir_model_dir["sib"][model].select_field("variant_label="+sample.get_property("variant_label"))[1548:].subspace(T=cf.djf())
			al=wci_component_dir_model_dir["al"][model].select_field("variant_label="+sample.get_property("variant_label"))[1548:].subspace(T=cf.djf())
			mar=wci_component_dir_model_dir["mar"][model].select_field("variant_label="+sample.get_property("variant_label"))[1548:].subspace(T=cf.djf())
			wci,wci_sib,wci_al,wci_mar=loc_func.cal_wci(sib,al,mar)
			'''
			Save for AOD comparison
			'''
			dirout="/home/users/lguo/cssp_aero/data/cmip6/cmip/historical/wci"
			if not os.path.isdir(dirout):
				os.makedirs(dirout)

			fonm='{}/wci_{}_historical_{}_gn_djf19792014.nc'.format(dirout,model,wci.get_property("variant_label"))
			print(fonm)
			cf.write(wci,fonm)

			wci_lst.append(wci);wci_sib_lst.append(wci_sib);wci_al_lst.append(wci_al);wci_mar_lst.append(wci_mar)
	wci_lst_mm=loc_func.fldlst2lstmm(wci_lst)
	wci_sib_lst_mm=loc_func.fldlst2lstmm(wci_sib_lst)
	wci_al_lst_mm=loc_func.fldlst2lstmm(wci_al_lst)
	wci_mar_lst_mm=loc_func.fldlst2lstmm(wci_mar_lst)
#	'''
#	Write out wci_conc for ploting fig02
#	'''
#	with open('/home/users/lguo/cssp_aero/analysis_script/data/wci_historical_concatenated.npy','wb') as f:
#		np.save(f,wci_conc)
	'''
	Concatenation need to carry out 
	'''
	for i,mdlnm in enumerate(model_lst):
		wci_lst_mm_normalised=(wci_lst_mm[i]-wci_lst_mm[i].collapse('T: mean'))/wci_lst_mm[i].collapse('T: sd')
		cf.write(wci_lst_mm_normalised,'/home/users/lguo/cssp_aero/analysis_script/data/wci_historical_{:02d}.nc'.format(i))
		if i==0:
			wci_conc=wci_lst_mm_normalised.array
		else:
			wci_conc=np.concatenate((wci_conc,wci_lst_mm_normalised.array))
		wci_sib_lst_mm_normalised=(wci_sib_lst_mm[i]-wci_sib_lst_mm[i].collapse('T: mean'))/wci_sib_lst_mm[i].collapse('T: sd')
		cf.write(wci_sib_lst_mm_normalised,'/home/users/lguo/cssp_aero/analysis_script/data/wci_sib_historical_{:02d}.nc'.format(i))
		wci_al_lst_mm_normalised=(wci_al_lst_mm[i]-wci_al_lst_mm[i].collapse('T: mean'))/wci_al_lst_mm[i].collapse('T: sd')
		cf.write(wci_al_lst_mm_normalised,'/home/users/lguo/cssp_aero/analysis_script/data/wci_al_historical_{:02d}.nc'.format(i))
		wci_mar_lst_mm_normalised=(wci_mar_lst_mm[i]-wci_mar_lst_mm[i].collapse('T: mean'))/wci_mar_lst_mm[i].collapse('T: sd')
		cf.write(wci_mar_lst_mm_normalised,'/home/users/lguo/cssp_aero/analysis_script/data/wci_mar_historical_{:02d}.nc'.format(i))
	'''
	Write out wci_conc for ploting
	'''
	with open('/home/users/lguo/cssp_aero/analysis_script/data/wci_historical_concatenated.npy','wb') as f:
		np.save(f,wci_conc)
'''
Plot WCI histogram for obs and CMIP historical concatenation
'''
if False:
	wci_obs=cf.read('/home/users/lguo/cssp_aero/analysis_script/data/wci_obs.nc')[0]
	with open('/home/users/lguo/cssp_aero/analysis_script/data/wci_historical_concatenated.npy','rb') as f:
		wci_cmip=np.load(f)
	'''
	Kolmogorov-Smirnov test: H0: Two samples are drawn from the same distribution.
	'''
	d,pvalue=stats.ks_2samp(wci_obs.array,wci_cmip)
	print("The P-value of Kolmogorov-Smirnov test is:",pvalue)

	kwargs=dict(bins=50,hist_kws={'alpha':.6},kde_kws={'linewidth':2})
	sns.distplot(wci_obs*-1.,color="dodgerblue",label="ERA5",**kwargs)
	sns.distplot(wci_cmip*-1.,color="orange",label="CMIP6",**kwargs)
	plt.ylim(0,1.)
	plt.xlim(-5,5)
	plt.xlabel(r"WCI$^*$",fontsize=15)
	plt.ylabel("Probability Density Function",fontsize=15)
	plt.gca().tick_params(labelsize="15",direction="out",bottom=True,left=True)
	plt.legend(fontsize=15)
	plt.savefig('./plot/hwi_scenariomip/fig_wci_historical_pdf_matched.pdf',bbox_inches='tight',pad_inches=0)
	plt.show()
	plt.close()
'''
Calculate WCI 
	for (not just matched) models and samples in ScenarioMIP
'''
if False:
	mip='ScenarioMIP'
	scenarios=scenarios[1:]
	mip_var="psl"
	wci_component=["sib","al","mar"]
	'''
	Read historical
	'''
	historical_wci_component_dict={}
	for j,component in enumerate(wci_component):
		print("historical",component)
		historical_wci_component_dict[component]=cf.read("/home/users/lguo/cssp_aero/data/cmip6/cmip/historical/{0}/{1}/{0}_Amon_*_historical_*_gn_185001-201412.nc".format(mip_var,component))
	'''
	Read Scenario
	'''
	scenario_dict_wci_component_dict={}
	for i,scenario in enumerate(scenarios):
		scenario_dict_wci_component_dict[scenario]={}
		for j,component in enumerate(wci_component):
			print(scenario,component)
			scenario_dict_wci_component_dict[scenario][component]=cf.read('/home/users/lguo/cssp_aero/data/cmip6/{0}/{1}/{2}/{3}/{2}_Amon_*_{1}_*_gn_201501-210012.nc'.format(mip,scenario,mip_var,component))
	'''
	Flatten 2D FieldList to 1D
	'''
	for i,scenario in enumerate(scenarios):
		wci_scenario_lst=[];wci_sib_scenario_lst=[];wci_al_scenario_lst=[];wci_mar_scenario_lst=[]
		for sample in scenario_dict_wci_component_dict[scenario]["sib"]:
			try:
				print(scenario,sample.get_property("source_id"),sample.get_property("variant_label"))
				sib=scenario_dict_wci_component_dict[scenario]["sib"].select("source_id="+sample.get_property("source_id")).select_field("variant_label="+sample.get_property("variant_label")).subspace(T=cf.djf())
				al=scenario_dict_wci_component_dict[scenario]["al"].select("source_id="+sample.get_property("source_id")).select_field("variant_label="+sample.get_property("variant_label")).subspace(T=cf.djf())
				mar=scenario_dict_wci_component_dict[scenario]["mar"].select("source_id="+sample.get_property("source_id")).select_field("variant_label="+sample.get_property("variant_label")).subspace(T=cf.djf())

				print("historical",sample.get_property("source_id"),sample.get_property("variant_label"))
				sib_historical=historical_wci_component_dict["sib"].select("source_id="+sample.get_property("source_id")).select_field("variant_label="+sample.get_property("variant_label")).subspace(T=cf.wi(cf.dt("1979-01-01"),cf.dt("2014-12-30"))).subspace(T=cf.djf())
				al_historical=historical_wci_component_dict["al"].select("source_id="+sample.get_property("source_id")).select_field("variant_label="+sample.get_property("variant_label")).subspace(T=cf.wi(cf.dt("1979-01-01"),cf.dt("2014-12-30"))).subspace(T=cf.djf())
				mar_historical=historical_wci_component_dict["mar"].select("source_id="+sample.get_property("source_id")).select_field("variant_label="+sample.get_property("variant_label")).subspace(T=cf.wi(cf.dt("1979-01-01"),cf.dt("2014-12-30"))).subspace(T=cf.djf())

				wci,wci_sib,wci_al,wci_mar=loc_func.cal_wci_scenario(sib,sib_historical,al,al_historical,mar,mar_historical)
				wci_scenario_lst.append(wci);wci_sib_scenario_lst.append(wci_sib);wci_al_scenario_lst.append(wci_al);wci_mar_scenario_lst.append(wci_mar)
			except ValueError:
				print("... One component in either Scenario or Historical is not available .....")
				continue
			'''
			Save for AOD comparison
			'''
			dirout="/home/users/lguo/cssp_aero/data/cmip6/ScenarioMIP/{}/wci".format(scenario)
			if not os.path.isdir(dirout): os.makedirs(dirout)

			fonm='{}/wci_{}_{}_{}_gn_197901-210012.nc'.format(dirout,sib.get_property("source_id"),scenario,sib.get_property("variant_label"))
			print(fonm)
			cf.write(wci,fonm)

		wci_scenario_lst_mm=loc_func.fldlst2lstmm(wci_scenario_lst);wci_sib_scenario_lst_mm=loc_func.fldlst2lstmm(wci_sib_scenario_lst);wci_al_scenario_lst_mm=loc_func.fldlst2lstmm(wci_al_scenario_lst);wci_mar_scenario_lst_mm=loc_func.fldlst2lstmm(wci_mar_scenario_lst)
		'''
		Write out wci_lst_mm for ploting 
		'''
		for im,mm in enumerate(wci_scenario_lst_mm):
			mm_normalised=(mm-mm.collapse('T: mean'))/mm.collapse('T: sd')
			cf.write(mm_normalised,'/home/users/lguo/cssp_aero/analysis_script/data/wci_historical-{}_{:02d}.nc'.format(scenario,im))
			'''
			To put model name into the concatenated .nc file. model name needs to be record
			'''
			if im == 0:
				f=open("/home/users/lguo/cssp_aero/analysis_script/data/wci_historical-{}_mdlnm.txt".format(scenario),"w")
				f.write('"{}"s'.format(mm.get_property("source_id")))
			else:
				f=open("/home/users/lguo/cssp_aero/analysis_script/data/wci_historical-{}_mdlnm.txt".format(scenario),"a")
				f.write(',"{}"s'.format(mm.get_property("source_id")))
		for im,mm in enumerate(wci_sib_scenario_lst_mm):
			mm_normalised=(mm-mm.collapse('T: mean'))/mm.collapse('T: sd')
			cf.write(mm_normalised,'/home/users/lguo/cssp_aero/analysis_script/data/wci_sib_historical-{}_{:02d}.nc'.format(scenario,im))
		for im,mm in enumerate(wci_al_scenario_lst_mm):
			mm_normalised=(mm-mm.collapse('T: mean'))/mm.collapse('T: sd')
			cf.write(mm_normalised,'/home/users/lguo/cssp_aero/analysis_script/data/wci_al_historical-{}_{:02d}.nc'.format(scenario,im))
		for im,mm in enumerate(wci_mar_scenario_lst_mm):
			mm_normalised=(mm-mm.collapse('T: mean'))/mm.collapse('T: sd')
			cf.write(mm_normalised,'/home/users/lguo/cssp_aero/analysis_script/data/wci_mar_historical-{}_{:02d}.nc'.format(scenario,im))
'''
Visualise WCI changes in ScenarioMIP
'''
if False:
	scenarios=scenarios[1:]
	record_marker_colour=False
	var_titles=[r"WCI* anomaly",r"-1$\times$SLP$_{sib}$ anomaly",r"SLP$_{al}$ anomaly",r"SLP$_{mar}$ anomaly"]

#	for ivar,var in enumerate(["wci","wci_sib","wci_al","wci_mar"]):
	for ivar,var in enumerate(["wci"]):
		if var=="wci": isflip=True;ishalf=False
		if var=="wci_sib": isflip=True;ishalf=False
		if var=="wci_al": isflip=False;ishalf=True
		if var=="wci_mar": isflip=False;ishalf=True

		wci_all={}
		wci_scenario_ano={};wci_scenario_ano_data={} # it is a dictionary of scenarios; each entry contains a 2D list with 0d as time slice and 1d as mm of HWIano.
		wci_scenario={};wci_scenario_data={}
		wci_historical={};wci_historical_data={}
		for scenario in scenarios:
			wci_all[scenario]=cf.read('/home/users/lguo/cssp_aero/analysis_script/data/mm_{}_historical-{}.nc'.format(var,scenario))[0]

			wci_scenario_ano[scenario],wci_scenario[scenario],wci_historical[scenario]=loc_func.hwi_anomaly3(wci_all[scenario],flip=isflip,half=ishalf)

			wci_scenario_ano_data[scenario]=loc_func.fldlst2list_2level(wci_scenario_ano[scenario])
			wci_scenario_data[scenario]=loc_func.fldlst2list_2level(wci_scenario[scenario])
			wci_historical_data[scenario]=loc_func.fldlst2list_2level(wci_historical[scenario])
		'''
		Make boxplot
		'''
		boxcolors=['#1E9684','#1D3354','#EADD3D','#F21111','#840B22']
		boxcolors=boxcolors[1:]
		medianprops=dict(linewidth=2,color='k')
		box={}
		for i,scenario in enumerate(scenarios):
			boxprops=dict(color=boxcolors[i],facecolor="w",linewidth=2,alpha=.7)
			pos=[-.5+i*.25, 1.5+i*.25, 3.5+i*.25, 5.5+i*.25]
			box[scenario]=plt.boxplot(wci_scenario_ano_data[scenario],whis=(5,95),widths=.15,positions=pos,showcaps=True,showfliers=False,patch_artist=True,medianprops=medianprops,boxprops=boxprops)
			'''
			Conduct ttest
			'''
			for ip,po in enumerate(pos):
				sig=loc_func.ttest4mean(np.nanmean(wci_scenario_data[scenario][ip]),np.nanmean(wci_historical_data[scenario][ip]),np.nanstd(wci_scenario_data[scenario][ip]),np.nanstd(wci_historical_data[scenario][ip]),len(wci_scenario_data[scenario][ip]),len(wci_historical_data[scenario][ip]))
				print(scenario,ip,sig)
				if sig<.1:
					print("plotting ...")
					boxprops=dict(color=boxcolors[i],facecolor=boxcolors[i],linewidth=2,alpha=.7)
					plt.boxplot(wci_scenario_ano_data[scenario][ip],whis=(5,95),widths=.15,positions=[pos[ip]],showcaps=True,showfliers=False,patch_artist=True,medianprops=medianprops,boxprops=boxprops)
		'''
		Add individual model
		'''
		markers=[".",",","o","v","^","<",">","1","2","3","4","8","s","p","P","*","h","H","+","x","X","D","d"]
		model_lst=[]
		model_lgd=[]
		model_mkr={};imkr=0
		model_clr={};iclr=0
		for i,scenario in enumerate(scenarios):
			pos=[-.5+i*.25, 1.5+i*.25, 3.5+i*.25, 5.5+i*.25]
			for j,wcia_period in enumerate(wci_scenario_ano[scenario]):
				for wcia in wcia_period:
					mdlnm=wcia.coordinates("source_id").value().array[0].decode("utf-8")
					print(scenario,mdlnm,wcia.array[0])
					'''
					True for recreating marker colors
					'''
					if record_marker_colour:
						if mdlnm in model_lst:
							plt.scatter(pos[j],wcia.array[0],s=15,c=model_clr[mdlnm],marker=model_mkr[mdlnm])
						else:
							model_lst.append(mdlnm)
							model_mkr[mdlnm]=markers[imkr];imkr+=1
							'''
							Generate RBG color randomly
							'''
							mkr_rbg=np.random.choice(range(256),size=3)
							'''
							Convert RBG to hex color, as RBG sequence is indistinguishable from an array of values to be colormapped
							'''
							model_clr[mdlnm]='#{:02x}{:02x}{:02x}'.format(mkr_rbg[0],mkr_rbg[1],mkr_rbg[2])
							model_lgd.append(plt.scatter(pos[j],wcia.array[0],s=15,c=model_clr[mdlnm],marker=model_mkr[mdlnm]))
					else:
						'''
						Read in recorded colours
						'''
						with open('/home/users/lguo/cssp_aero/analysis_script/data/marker_colour.pickle','rb') as f:
							model_clr=pickle.load(f)
	
						if mdlnm in model_lst:
							plt.scatter(pos[j],wcia.array[0],s=15,c=model_clr[mdlnm],marker=model_mkr[mdlnm])
						else:
							try:
								colour_key=model_clr[mdlnm]
							except KeyError:
								print("{} is not in Model List ...".format(mdlnm))
								continue
							model_lst.append(mdlnm)
							model_mkr[mdlnm]=markers[imkr];imkr+=1
							model_lgd.append(plt.scatter(pos[j],wcia.array[0],s=15,c=model_clr[mdlnm],marker=model_mkr[mdlnm]))
	
		'''
		Record maker colours for future use
		'''
		if record_marker_colour:
			with open('/home/users/lguo/cssp_aero/analysis_script/data/marker_colour.pickle','wb') as f:
				pickle.dump(model_clr,f)
		'''
		Make double legends:
			1. build legend0;
			2. plot legend1;
			3. plot legend0;
		'''
		scenarios_str=['SSP 1-1.9','SSP 1-2.6','SSP 2-4.5','SSP 3-7.0','SSP 5-8.5']
		scenarios_str=scenarios_str[1:]
		legend0=plt.legend([box[scenario]["boxes"][0] for scenario in scenarios],scenarios_str,bbox_to_anchor=(.02,.81,.05,.12),loc='upper left',fontsize=15)
		plt.legend(model_lgd,model_lst,bbox_to_anchor=(0,.95,1,0.3),loc="lower left",mode="expand",ncol=4,fontsize='xx-small')
		if ivar==0:
			plt.gca().add_artist(legend0)
		'''
		x&y axes
		'''
		plt.xticks([0,2,4,6],['2025-2034','2035-2044','2045-2054','2090-2099'])
		plt.ylim(-1,2.5)
		plt.xlim(-1,7)
		plt.plot([-1,7],[0,0],c='k',lw=1,ls=':')
		plt.ylabel(var_titles[ivar],fontsize=20)
		plt.gca().tick_params(labelsize="15",direction="in",bottom=True,left=True)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().spines['top'].set_visible(False)
		'''
		output
		'''
		plt.savefig('./plot/hwi_scenariomip/fig_{}_scenario_anomaly.pdf'.format(var),bbox_inches='tight',pad_inches=0)
	#	plt.show()
		plt.close()
'''
Visualise WCI changes in ScenarioMIP
	Single model
'''
if False:
	scenarios=scenarios[1:]

	model="CanESM5"
	var="wci";isflip=True;ishalf=False
#	var="wci_sib";isflip=True;ishalf=False
#	var="wci_al";isflip=False;ishalf=True
#	var="wci_mar";isflip=False;ishalf=True

	wci_all={};mm_all={};
	wci_scenario_ano={};wci_scenario_ano_data={};mm_scenario_ano={};mm_scenario_ano_data={};
	wci_scenario={};wci_scenario_data={};mm_scenario={};mm_scenario_data={};
	wci_historical={};wci_historical_data={};mm_historical={};mm_historical_data={};
	for scenario in scenarios:
		wci_all[scenario]=loc_func.aggregate_fieldlist(cf.read("{0}/ScenarioMIP/{1}/wci/wci_{2}_{1}_*_gn_197901-210012.nc".format(cmipdir,scenario,model)))

		wci_scenario_ano[scenario],wci_scenario[scenario],wci_historical[scenario],=loc_func.hwi_anomaly3(wci_all[scenario],flip=isflip,half=ishalf)

		wci_scenario_ano_data[scenario]=loc_func.fldlst2list_2level(wci_scenario_ano[scenario])
		wci_scenario_data[scenario]=loc_func.fldlst2list_2level(wci_scenario[scenario])
		wci_historical_data[scenario]=loc_func.fldlst2list_2level(wci_historical[scenario])

		mm_all[scenario]=cf.read('/home/users/lguo/cssp_aero/analysis_script/data/mm_{}_historical-{}.nc'.format(var,scenario))[0]

		mm_scenario_ano[scenario],mm_scenario[scenario],mm_historical[scenario]=loc_func.hwi_anomaly3(mm_all[scenario],flip=isflip,half=ishalf)

		mm_scenario_ano_data[scenario]=loc_func.fldlst2list_2level(mm_scenario_ano[scenario])
		mm_scenario_data[scenario]=loc_func.fldlst2list_2level(mm_scenario[scenario])
		mm_historical_data[scenario]=loc_func.fldlst2list_2level(mm_historical[scenario])
	'''
	Make boxplot
	'''
	boxcolors=['#1E9684','#1D3354','#EADD3D','#F21111','#840B22']
	boxcolors=boxcolors[1:]
	medianprops=dict(linewidth=2,color='k')
	box={}
	for i,scenario in enumerate(scenarios):
		boxprops=dict(color=boxcolors[i],facecolor="w",linewidth=2,alpha=.7)
		pos=[-.5+i*.25, 1.5+i*.25, 3.5+i*.25, 5.5+i*.25]
		box[scenario]=plt.boxplot(wci_scenario_ano_data[scenario],positions=pos,whis=(5,95),widths=.15,showcaps=True,showfliers=False,patch_artist=True,boxprops=boxprops,medianprops=medianprops)
		'''
		Conduct ftest
		'''
		for ip,po in enumerate(pos):
			sig=loc_func.ftest(np.nanstd(mm_scenario_data[scenario][ip]),np.nanstd(wci_scenario_data[scenario][ip]),len(mm_scenario_data[scenario][ip]),len(wci_scenario_data[scenario][ip]))
			print(scenario,ip,sig)
			if sig<.1:
				print("plotting ...")
				boxprops=dict(color=boxcolors[i],facecolor=boxcolors[i],linewidth=2,alpha=.7)
				plt.boxplot(wci_scenario_ano_data[scenario][ip],whis=(5,95),positions=[pos[ip]],widths=.15,showcaps=True,showfliers=False,patch_artist=True,medianprops=medianprops,boxprops=boxprops)
	'''
	Add individual model
	'''
	for i,scenario in enumerate(scenarios):
		pos=[-.5+i*.25, 1.5+i*.25, 3.5+i*.25, 5.5+i*.25]
		for j,wcia_period in enumerate(wci_scenario_ano[scenario]):
			for wcia in wcia_period:
				plt.scatter(pos[j],wcia.array[0],s=15,c="k",marker=".",alpha=.6)
	'''
	Legend
	'''
	scenarios_str=['SSP 1-1.9','SSP 1-2.6','SSP 2-4.5','SSP 3-7.0','SSP 5-8.5']
	scenarios_str=scenarios_str[1:]
	plt.legend([box[snr]["boxes"][0] for snr in scenarios],scenarios_str,title=model,bbox_to_anchor=(.02,.89,.05,.12),loc='upper left')
	'''
	x&y axes
	'''
	plt.xticks([0,2,4,6],['2025-2034','2035-2044','2045-2054','2090-2099'])
	plt.ylim(-1,2.5)
	plt.xlim(-1,7)
	plt.gca().spines['right'].set_visible(False)
	plt.gca().spines['top'].set_visible(False)
	plt.gca().yaxis.set_ticks_position('left')
	plt.gca().xaxis.set_ticks_position('bottom')
	plt.gca().set(ylabel=r'-1$\times${} Anomaly'.format(var.upper()))
	'''
	output
	'''
	plt.savefig('./plot/hwi_scenariomip/fig_{}_scenario_anomaly_{}.pdf'.format(var,model),bbox_inches='tight',pad_inches=0)
	plt.show()
	plt.close()
'''
Plot HWI components between historical and obs: ua@500hPa, va@850hPa and ta@250&850hPa
'''
if False:
	mip_var=["ua","va","ta","ta"][:2]
	plev=["500","850","850","250"]
	'''
	Set domain
	'''
	lati=-25;latx=65;loni=60;lonx=180
	'''
	for indicating HWI domains
	'''
	box_xpt={'ua500':[[110.,110.,137.5,137.5,110.],[110.,110.,137.5,137.5,110.]]
		,'va850':[115.,115.,130.,130.,115.]
		,'ta850':[112.5,112.5,132.5,132.5,112.5]
		,'ta250':[122.5,122.5,137.5,137.5,122.5]}
	box_ypt={'ua500':[[22.7,37.5,37.5,22.7,22.7],[42.5,52.5,52.5,42.5,42.5]]
		,'va850':[30.,47.5,47.5,30.,30.]
		,'ta850':[32.5,45,45,32.5,32.5]
		,'ta250':[37.5,45,45,37.5,37.5]}
	'''
	Read obs:
		HWI components are read into a dictionary;
	'''
	hwi_comp_obs={}
	for i,var in enumerate(mip_var):
		hwi_comp_obs[var+plev[i]]=cf.read('{}/{}{}.mon.interp.nc'.format(obsdir,var,plev[i]))[0][:432].subspace(T=cf.djf()).collapse('mean','T').squeeze()
	'''
	Plot obs
	'''
	cfp.setvars(axis_label_fontsize=20,colorbar_fontsize=20)
	if True:
		for i,var in enumerate(mip_var):
			cfp.gopen()
			'''
			for publication, colour scale need to be consistent between obs and mmm
			'''
			if var=="ua":
				cfp.levs(min=-8,max=36,step=2)
			if var=="va":
				cfp.levs(min=-7,max=7,step=1)
			if var=="ta" and plev[i]=="850":
				cfp.levs(min=246,max=296,step=2)
			if var=="ta" and plev[i]=="250":
				cfp.levs(min=210,max=232,step=2)
			cfp.mapset(lonmin=loni,lonmax=lonx,latmin=lati,latmax=latx)
			cfp.con(hwi_comp_obs[var+plev[i]],lines=False,colorbar_title="m/s")
			cfp.plotvars.mymap.tick_params(direction="out",bottom=True,left=True,color="k")
			'''
			add boxes
			'''
			if var=='ua':
				cfp.plotvars.mymap.plot(box_xpt[var+plev[i]][0],box_ypt[var+plev[i]][0],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
				cfp.plotvars.mymap.plot(box_xpt[var+plev[i]][1],box_ypt[var+plev[i]][1],linewidth=2.,color='red',transform=ccrs.PlateCarree())
			else:
				cfp.plotvars.mymap.plot(box_xpt[var+plev[i]],box_ypt[var+plev[i]],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
			cfp.plotvars.master_plot.savefig('./plot/hwi_scenariomip/obs_{}{}_19792014djf.pdf'.format(var,plev[i]),bbox_inches='tight',pad_inches=0)
			cfp.gclose()
	'''
	Read mm and make mmm:
		Since cf2&3 are not compatible and 3 has issue of aggregation, this procedure is carried out in a stand along cf2 script: hwi_cmip_comp.py
		mm are made in var_hwi_trunk.py
	'''
	hwi_comp_cmip_mean={}
	hwi_comp_cmip_sd={}
	for i,var in enumerate(mip_var):
		hwi_comp_cmip_mean[var+plev[i]]=cf.read('{0}/cmip/historical/{1}/{2}hPa/mm_{1}_Amon_historical_mm_gn_djf.nc'.format(cmipdir,var,plev[i]))[0].collapse('mean','sample').squeeze()
		hwi_comp_cmip_sd[var+plev[i]]  =cf.read('{0}/cmip/historical/{1}/{2}hPa/mm_{1}_Amon_historical_mm_gn_djf.nc'.format(cmipdir,var,plev[i]))[0].collapse('sd'  ,'sample').squeeze()
	'''
	Plot cmip
	'''
	if True:
		for i,var in enumerate(mip_var):
			cfp.gopen()
			'''
			for publication, colour scale need to be consistent between obs and mmm
			'''
			if var=="ua":
				cfp.levs(min=-8,max=36,step=2)
			if var=="va":
				cfp.levs(min=-7,max=7,step=1)
			if var=="ta" and plev[i]=="850":
				cfp.levs(min=246,max=296,step=2)
			if var=="ta" and plev[i]=="250":
				cfp.levs(min=210,max=232,step=2)
			cfp.mapset(lonmin=loni,lonmax=lonx,latmin=lati,latmax=latx)
			cfp.con(hwi_comp_cmip_mean[var+plev[i]],lines=False,colorbar_title="m/s")
			cfp.plotvars.mymap.tick_params(direction="out",bottom=True,left=True,color="k")
			'''
			add boxes
			'''
			if var=='ua':
				cfp.plotvars.mymap.plot(box_xpt[var+plev[i]][0],box_ypt[var+plev[i]][0],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
				cfp.plotvars.mymap.plot(box_xpt[var+plev[i]][1],box_ypt[var+plev[i]][1],linewidth=2.,color='red',transform=ccrs.PlateCarree())
			else:
				cfp.plotvars.mymap.plot(box_xpt[var+plev[i]],box_ypt[var+plev[i]],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
			cfp.plotvars.master_plot.savefig('./plot/hwi_scenariomip/cmip_{}{}_19792014djf.pdf'.format(var,plev[i]),bbox_inches='tight',pad_inches=0)
			cfp.gclose(view=False)
	'''
	Plot difference between cmip mmm and obs
	'''
	if True:
		hwi_comp_diff={}
		hwi_comp_diff_pvalue={}
		for i,var in enumerate(mip_var):
			hwi_comp_diff[var+plev[i]]=hwi_comp_cmip_mean[var+plev[i]]-hwi_comp_obs[var+plev[i]]
			hwi_comp_diff_pvalue[var+plev[i]]=loc_func.ttest4mean_1sample(hwi_comp_cmip_mean[var+plev[i]],hwi_comp_obs[var+plev[i]],hwi_comp_cmip_sd[var+plev[i]],17)
		for i,var in enumerate(mip_var):
			cfp.gopen()
			if var=="ua":
				cfp.levs(min=-4,max=4,step=.5)
			if var=="va":
				cfp.levs(min=-4,max=4,step=1)
			if var=="ta" and plev[i]=="850":
				cfp.levs(min=-7,max=4,step=1)
			if var=="ta" and plev[i]=="250":
				cfp.levs(min=-8,max=1,step=.5)
			cfp.mapset(lonmin=loni,lonmax=lonx,latmin=lati,latmax=latx)
			cfp.con(hwi_comp_diff[var+plev[i]],lines=False,colorbar_title="m/s")
			cfp.stipple(hwi_comp_diff_pvalue[var+plev[i]],min=1e-100,max=.10,size=5,color='k',pts=100,marker='.')
			cfp.plotvars.mymap.tick_params(direction="out",bottom=True,left=True,color="k")
			'''
			add boxes
			'''
			if var=='ua':
				cfp.plotvars.mymap.plot(box_xpt[var+plev[i]][0],box_ypt[var+plev[i]][0],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
				cfp.plotvars.mymap.plot(box_xpt[var+plev[i]][1],box_ypt[var+plev[i]][1],linewidth=2.,color='red',transform=ccrs.PlateCarree())
			else:
				cfp.plotvars.mymap.plot(box_xpt[var+plev[i]],box_ypt[var+plev[i]],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
			cfp.plotvars.master_plot.savefig('./plot/hwi_scenariomip/diff_{}{}_19792014djf.pdf'.format(var,plev[i]),bbox_inches='tight',pad_inches=0)
			cfp.gclose(view=False)
'''
Plot HWI components between scenarios and historical
'''
if False:
	mip_var=["ua","va"]
	plev=["500","850"]
	'''
	Set domain
	'''
	lati=-25;latx=65;loni=60;lonx=180
	'''
	for indicating HWI domains
	'''
	box_xpt={'ua500':[[110.,110.,137.5,137.5,110.],[110.,110.,137.5,137.5,110.]]
		,'va850':[115.,115.,130.,130.,115.]
		,'ta850':[112.5,112.5,132.5,132.5,112.5]
		,'ta250':[122.5,122.5,137.5,137.5,122.5]}
	box_ypt={'ua500':[[22.7,37.5,37.5,22.7,22.7],[42.5,52.5,52.5,42.5,42.5]]
		,'va850':[30.,47.5,47.5,30.,30.]
		,'ta850':[32.5,45,45,32.5,32.5]
		,'ta250':[37.5,45,45,37.5,37.5]}
	'''
	Loop through variables
	'''
	for i,var in enumerate(mip_var):
		'''
		Read historical
		'''
		fld_historical=cf.read('{0}/cmip/historical/{1}/{2}hPa/mm_{1}_Amon_historical_mm_gn_djf.nc'.format(cmipdir,var,plev[i]))[0]
		'''
		Loop through Scenarios
		'''
		for scenario in scenarios[1:]:
			for period in ["2025-2034","2035-2044","2045-2054","2090-2099"]:
				'''
				Read Scenario
				'''
				fld_scenario=cf.read('{0}/ScenarioMIP/{1}/{2}/{3}hPa/{4}/mm_{2}_Amon_{1}_mm_gn_djf.nc'.format(cmipdir,scenario,var,plev[i],period))[0]
				var_historical_mean,var_historical_sd,var_scenario_mean,var_scenario_sd,nsample=loc_func.matched_model_mean(fld_historical,fld_scenario)
				'''
				Plot difference
				'''
				cfp.setvars(axis_label_fontsize=20,colorbar_fontsize=20)
				if True:
					var_diff=var_scenario_mean-var_historical_mean
					var_diff_pvalue=loc_func.ttest4mean(var_scenario_mean,var_historical_mean,var_scenario_sd,var_historical_sd,nsample,nsample)

					cfp.gopen(file="./plot/hwi_scenariomip/diff_{1}{2}_{0}_{3}djf.pdf".format(scenario,var,plev[i],period))
					cfp.mapset(lonmin=loni,lonmax=lonx,latmin=lati,latmax=latx)
					cfp.levs(min=-4,max=4,step=.5)
					cfp.con(var_diff,lines=False,title="{} {}: {}".format(var,scenario,period),colorbar_title='m/s')
					cfp.stipple(var_diff_pvalue,min=1e-100,max=.10,size=5,color='k',pts=100,marker='.')
					cfp.plotvars.mymap.tick_params(direction="out",bottom=True,left=True,color="k")
					'''
					add boxes
					'''
					if var=='ua':
						cfp.plotvars.mymap.plot(box_xpt[var+plev[i]][0],box_ypt[var+plev[i]][0],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
						cfp.plotvars.mymap.plot(box_xpt[var+plev[i]][1],box_ypt[var+plev[i]][1],linewidth=2.,color='red',transform=ccrs.PlateCarree())
					else:
						cfp.plotvars.mymap.plot(box_xpt[var+plev[i]],box_ypt[var+plev[i]],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
					cfp.gclose()
'''
Plot psl between historical and obs
'''
if False:
	mip_var="psl"
	'''
	Set domain
	'''
	lati=-25;latx=65;loni=60;lonx=180
	'''
	for indicating HWI domains
	'''
	box_xpt={'sib':[70.,70.,120.,120.,70.],'al':[140.,140.,170.,170.,140.],'mar':[110.,110.,160.,160.,110.]}
	box_ypt={'sib':[40.,60.,60.,40.,40.],'al':[30.,50.,50.,30.,30.],'mar':[-20.,10.,10.,-20.,-20.]}
	'''
	Read obs:
	'''
	psl_obs=cf.read('{}/{}.mon.interp.nc'.format(obsdir,mip_var))[0][:432].subspace(T=cf.djf()).collapse('mean','T').squeeze()
	'''
	Plot obs
	'''
	cfp.setvars(axis_label_fontsize=20,colorbar_fontsize=20)
	if True:
		cfp.gopen()
		'''
		for publication, colour scale need to be consistent between obs and mmm
		'''
		cfp.mapset(lonmin=loni,lonmax=lonx,latmin=lati,latmax=latx)
		cfp.levs(min=982,max=1034,step=2)
		cfp.con(psl_obs/100.,lines=False,colorbar_title='hPa')
		cfp.plotvars.mymap.tick_params(direction="out",bottom=True,left=True,color="k")
		'''
		add boxes
		'''
		cfp.plotvars.mymap.plot(box_xpt["sib"],box_ypt["sib"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
		cfp.plotvars.mymap.plot(box_xpt["al"],box_ypt["al"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
		cfp.plotvars.mymap.plot(box_xpt["mar"],box_ypt["mar"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
		cfp.plotvars.master_plot.savefig('./plot/hwi_scenariomip/obs_{}_19792014djf.pdf'.format(mip_var),bbox_inches='tight',pad_inches=0)
		cfp.gclose(view=False)
	'''
	Read mm and make mmm:
		Since cf2&3 are not compatible and 3 has issue of aggregation, this procedure is carried out in a stand along cf2 script: hwi_cmip_comp.py
		mm are made in var_hwi_trunk.py
	'''
	psl_cmip_mean=cf.read('{0}/cmip/historical/{1}/slice/1979-2014/mm_{1}_Amon_historical_mm_gn_djf.nc'.format(cmipdir,mip_var))[0].collapse('mean','sample').squeeze()
	psl_cmip_sd  =cf.read('{0}/cmip/historical/{1}/slice/1979-2014/mm_{1}_Amon_historical_mm_gn_djf.nc'.format(cmipdir,mip_var))[0].collapse('sd'  ,'sample').squeeze()
	'''
	Plot cmip
	'''
	if True:
		cfp.gopen()
		cfp.mapset(lonmin=loni,lonmax=lonx,latmin=lati,latmax=latx)
		cfp.con(psl_cmip_mean/100.,lines=False,colorbar_title='hPa')
		cfp.plotvars.mymap.tick_params(direction="out",bottom=True,left=True,color="k")
		'''
		add boxes
		'''
		cfp.plotvars.mymap.plot(box_xpt["sib"],box_ypt["sib"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
		cfp.plotvars.mymap.plot(box_xpt["al"],box_ypt["al"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
		cfp.plotvars.mymap.plot(box_xpt["mar"],box_ypt["mar"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
		cfp.plotvars.master_plot.savefig('./plot/hwi_scenariomip/cmip_{}_19792014djf.pdf'.format(mip_var),bbox_inches='tight',pad_inches=0)
		cfp.gclose(view=False)
	'''
	Plot difference between cmip mmm and obs
	'''
	if True:
		psl_diff=psl_cmip_mean-psl_obs
		psl_diff_pvalue=loc_func.ttest4mean_1sample(psl_cmip_mean,psl_obs,psl_cmip_sd,17)
		cfp.gopen()
		cfp.levs()
		cfp.mapset(lonmin=loni,lonmax=lonx,latmin=lati,latmax=latx)
		cfp.con(psl_diff/100.,lines=False,colorbar_title='hPa')
		cfp.stipple(psl_diff_pvalue,min=1e-100,max=.10,size=5,color='k',pts=100,marker='.')
		cfp.plotvars.mymap.tick_params(direction="out",bottom=True,left=True,color="k")
		'''
		add boxes
		'''
		cfp.plotvars.mymap.plot(box_xpt["sib"],box_ypt["sib"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
		cfp.plotvars.mymap.plot(box_xpt["al"],box_ypt["al"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
		cfp.plotvars.mymap.plot(box_xpt["mar"],box_ypt["mar"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
		cfp.plotvars.master_plot.savefig('./plot/hwi_scenariomip/diff_{}_19792014djf.pdf'.format(mip_var),bbox_inches='tight',pad_inches=0)
		cfp.gclose(view=False)
'''
Plot psl between scenario and historical
'''
if False:
	mip_var="psl"
	'''
	Set domain
	'''
	lati=-25;latx=65;loni=60;lonx=180
	'''
	for indicating WCI domains
	'''
	box_xpt={'sib':[70.,70.,120.,120.,70.],'al':[140.,140.,170.,170.,140.],'mar':[110.,110.,160.,160.,110.]}
	box_ypt={'sib':[40.,60.,60.,40.,40.],'al':[30.,50.,50.,30.,30.],'mar':[-20.,10.,10.,-20.,-20.]}
	'''
	Read historical
	'''
	fld_historical=cf.read('{0}/cmip/historical/{1}/slice/1979-2014/mm_{1}_Amon_historical_mm_gn_djf.nc'.format(cmipdir,mip_var))[0]
	'''
	Loop through Scenarios
	'''
	for scenario in scenarios[1:]:
		for period in ["2025-2034","2035-2044","2045-2054","2090-2099"]:
			'''
			Read Scenario
			'''
			fld_scenario=cf.read('{0}/ScenarioMIP/{1}/{2}/slice/{3}/mm_{2}_Amon_{1}_mm_gn_djf.nc'.format(cmipdir,scenario,mip_var,period))[0]
			psl_historical_mean,psl_historical_sd,psl_scenario_mean,psl_scenario_sd,nsample=loc_func.matched_model_mean(fld_historical,fld_scenario)

			cfp.setvars(axis_label_fontsize=20,colorbar_fontsize=20)
			'''
			Plot historical
			'''
			if False:
				cfp.gopen()
				cfp.mapset(lonmin=loni,lonmax=lonx,latmin=lati,latmax=latx)
				cfp.con(psl_historical_mean/100.,lines=False,colorbar_title='{} (hPa)'.format(mip_var))
				'''
				add boxes
				'''
				cfp.plotvars.mymap.plot(box_xpt["sib"],box_ypt["sib"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
				cfp.plotvars.mymap.plot(box_xpt["al"],box_ypt["al"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
				cfp.plotvars.mymap.plot(box_xpt["mar"],box_ypt["mar"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
		#		cfp.plotvars.master_plot.savefig('./plot/hwi_scenariomip/cmip_{}_19792014djf.pdf'.format(mip_var),bbox_inches='tight',pad_inches=0)
				cfp.gclose()
			'''
			Plot Scenario
			'''
			if False:
				cfp.gopen()
				cfp.mapset(lonmin=loni,lonmax=lonx,latmin=lati,latmax=latx)
				cfp.con(psl_scenario_mean/100.,lines=False,colorbar_title='{} (hPa)'.format(mip_var))
				'''
				add boxes
				'''
				cfp.plotvars.mymap.plot(box_xpt["sib"],box_ypt["sib"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
				cfp.plotvars.mymap.plot(box_xpt["al"],box_ypt["al"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
				cfp.plotvars.mymap.plot(box_xpt["mar"],box_ypt["mar"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
		#		cfp.plotvars.master_plot.savefig('./plot/hwi_scenariomip/cmip_{}_19792014djf.pdf'.format(mip_var),bbox_inches='tight',pad_inches=0)
				cfp.gclose()
			'''
			Plot difference between scenario and historical
			'''
			if True:
				psl_diff=psl_scenario_mean-psl_historical_mean
				psl_diff_pvalue=loc_func.ttest4mean(psl_scenario_mean,psl_historical_mean,psl_scenario_sd,psl_historical_sd,nsample,nsample)
				cfp.gopen(file="./plot/hwi_scenariomip/diff_{1}_{0}_{2}djf.pdf".format(scenario,mip_var,period))
				cfp.mapset(lonmin=loni,lonmax=lonx,latmin=lati,latmax=latx)
				cfp.levs()
				cfp.con(psl_diff/100.,lines=False,title="{}: {}".format(scenario,period),colorbar_title='hPa')
				cfp.stipple(psl_diff_pvalue,min=1e-100,max=.10,size=5,color='k',pts=100,marker='.')
				cfp.plotvars.mymap.tick_params(direction="out",bottom=True,left=True,color="k")
				'''
				add boxes
				'''
				cfp.plotvars.mymap.plot(box_xpt["sib"],box_ypt["sib"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
				cfp.plotvars.mymap.plot(box_xpt["al"],box_ypt["al"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
				cfp.plotvars.mymap.plot(box_xpt["mar"],box_ypt["mar"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
				cfp.gclose()
			'''
			Plot the median of the difference
			'''
			if False:
				psl_diff_median=loc_func.matched_model_median(fld_historical,fld_scenario)

				cfp.gopen()
				cfp.levs()
				cfp.con(psl_diff_median/100.,lines=False,title="{}: {}".format(scenario,period),colorbar_title='{} (hPa)'.format(mip_var))
				'''
				add boxes
				'''
				cfp.plotvars.mymap.plot(box_xpt["sib"],box_ypt["sib"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
				cfp.plotvars.mymap.plot(box_xpt["al"],box_ypt["al"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
				cfp.plotvars.mymap.plot(box_xpt["mar"],box_ypt["mar"],linewidth=2.,color='blue',transform=ccrs.PlateCarree())
				cfp.plotvars.master_plot.savefig('./plot/hwi_scenariomip/diff_{}_{}_{}djf.pdf'.format(mip_var,scenario,period),bbox_inches='tight',pad_inches=0)
				cfp.gclose(view=False)
'''
Plot ta cross-section along 112.5E~137.5E
	Between Historical and Obs
'''
if False:
	mip_var="ta"
	box_xpt={'850':[32.5,45],'250':[37.5,45]}
	box_ypt={'850':[850,850],'250':[250,250]}
	'''
	Read obs:
		HWI components are read into a dictionary;
	'''
	hwi_comp_obs=cf.read("../data/obs/ta_pressure_mon_197901201912_interp.nc").select_field("air_temperature").subspace(T=cf.djf()).collapse("mean","T").subspace(latitude=cf.wi(0,70),longitude=cf.wi(112.5,137.5)).collapse("mean","X")
	'''
	Plot obs
	'''
	cfp.setvars(axis_label_fontsize=20,colorbar_fontsize=20)
	if True:
		cfp.gopen()
		'''
		for publication, colour scale need to be consistent between obs and mmm
		'''
		cfp.levs(min=200,max=300,step=10)
		cfp.con(hwi_comp_obs,lines=False,colorbar_title="K")
		cfp.plotvars.plot.tick_params(direction="out",bottom=True,left=True,color="k")
		cfp.plotvars.plot.set_ylabel("")
		cfp.plotvars.plot.set_xlabel("")
		'''
		add boxes
		'''
		cfp.plotvars.plot.plot(box_xpt["850"],box_ypt["850"],linewidth=2.,color='black')
		cfp.plotvars.plot.plot(box_xpt["250"],box_ypt["250"],linewidth=2.,color='black')
		cfp.plotvars.master_plot.savefig('./plot/hwi_scenariomip/obs_ta_pressure_19792014djf.pdf',bbox_inches='tight',pad_inches=0)
		cfp.gclose(view=False)
	'''
	Read mm and make mmm:
	'''
	hwi_comp_cmip_mean=cf.read("{0}/cmip/historical/{1}/multi-hPa/1979-2014/mm_{1}_Amon_historical_mm_gn_djf.nc".format(cmipdir,mip_var)).select_field("var_mean").collapse('mean','sample').subspace(latitude=cf.wi(0,70),longitude=cf.wi(112.5,137.5)).collapse("mean","X")
	hwi_comp_cmip_sd  =cf.read("{0}/cmip/historical/{1}/multi-hPa/1979-2014/mm_{1}_Amon_historical_mm_gn_djf.nc".format(cmipdir,mip_var)).select_field("var_mean").collapse('sd'  ,'sample').subspace(latitude=cf.wi(0,70),longitude=cf.wi(112.5,137.5)).collapse("mean","X")
	'''
	Plot cmip
	'''
	if True:
		cfp.gopen()
		cfp.levs(min=200,max=300,step=10)
		cfp.setvars(text_fontsize=20)
		cfp.con(hwi_comp_cmip_mean,lines=False,colorbar_title="K")
		cfp.plotvars.plot.tick_params(direction="out",bottom=True,left=True,color="k")
		cfp.plotvars.plot.set_ylabel("")
		cfp.plotvars.plot.set_xlabel("")
		'''
		add boxes
		'''
		cfp.plotvars.plot.plot(box_xpt["850"],box_ypt["850"],linewidth=2.,color='black')
		cfp.plotvars.plot.plot(box_xpt["250"],box_ypt["250"],linewidth=2.,color='black')
		cfp.plotvars.master_plot.savefig('./plot/hwi_scenariomip/cmip_ta_pressure_19792014djf.pdf',bbox_inches='tight',pad_inches=0)
		cfp.gclose(view=False)

	'''
	Plot difference between cmip mmm and obs
	'''
	if True:
		hwi_comp_obs.coordinate("Z").set_property("standard_name","air_pressure")
		hwi_comp_diff=hwi_comp_cmip_mean-hwi_comp_obs
		hwi_comp_diff_pvalue=loc_func.ttest4mean_1sample(hwi_comp_cmip_mean,hwi_comp_obs,hwi_comp_cmip_sd,17)

		cfp.gopen()
		cfp.levs()
		cfp.con(hwi_comp_diff,lines=False,colorbar_title="K")
		cfp.stipple(hwi_comp_diff_pvalue,min=1e-100,max=.10,size=5,color='k',pts=100,marker='.')
		cfp.plotvars.plot.tick_params(direction="out",bottom=True,left=True,color="k")
		cfp.plotvars.plot.set_ylabel("")
		cfp.plotvars.plot.set_xlabel("")
		'''
		add boxes
		'''
		cfp.plotvars.plot.plot(box_xpt["850"],box_ypt["850"],linewidth=2.,color='black')
		cfp.plotvars.plot.plot(box_xpt["250"],box_ypt["250"],linewidth=2.,color='black')
		cfp.plotvars.master_plot.savefig('./plot/hwi_scenariomip/diff_ta_pressure_19792014djf.pdf',bbox_inches='tight',pad_inches=0)
		cfp.gclose(view=False)
'''
Plot ta cross-section along 112.5E~137.5E
	Between Scenario and Historical
'''
if False:
	mip_var="ta"
	box_xpt={'850':[32.5,45],'250':[37.5,45]}
	box_ypt={'850':[850,850],'250':[250,250]}
	'''
	Read Historical
	'''
	fld_historical=cf.read("{0}/cmip/historical/{1}/multi-hPa/mm_{1}_Amon_historical_mm_gn_djf.nc".format(cmipdir,mip_var)).select_field("var_mean").subspace(latitude=cf.wi(0,70),longitude=cf.wi(112.5,137.5)).collapse("mean","X")
	'''
	Read Scenarios
	'''
	for scenario in scenarios[1:]:
		for period in ["2025-2034","2035-2044","2045-2054","2090-2099"]:
			'''
			Read Scenario
			'''
			fld_scenario=cf.read("{0}/ScenarioMIP/{1}/{2}/multi-hPa/{3}/mm_{2}_Amon_{1}_mm_gn_djf.nc".format(cmipdir,scenario,mip_var,period)).select_field("var_mean").subspace(latitude=cf.wi(0,70),longitude=cf.wi(112.5,137.5)).collapse("mean","X")
			fld_historical_mean,fld_historical_sd,fld_scenario_mean,fld_scenario_sd,nsample=loc_func.matched_model_mean(fld_historical,fld_scenario)
			'''
			Plot difference between cmip mmm and obs
			'''
			cfp.setvars(axis_label_fontsize=20,colorbar_fontsize=17)
			if True:
				fld_historical_mean.coordinate("Z").set_property("standard_name","air_pressure")
				fld_diff=fld_scenario_mean-fld_historical_mean
				fld_diff_pvalue=loc_func.ttest4mean(fld_scenario_mean,fld_historical_mean,fld_scenario_sd,fld_historical_sd,nsample,nsample)

				cfp.gopen(file="./plot/hwi_scenariomip/diff_{}_pressure_{}_{}djf.pdf".format(mip_var,scenario,period))
				cfp.levs(min=-5.,max=5.,step=.5)
				cfp.con(fld_diff,lines=False,title="{}: {}".format(scenario,period),colorbar_title="K")
				cfp.stipple(fld_diff_pvalue,min=1e-100,max=.01,size=5,color='k',pts=100,marker='.')
				cfp.plotvars.plot.tick_params(direction="out",bottom=True,left=True,color="k")
				cfp.plotvars.plot.set_ylabel("")
				cfp.plotvars.plot.set_xlabel("")
				'''
				add boxes
				'''
				cfp.plotvars.plot.plot(box_xpt["850"],box_ypt["850"],linewidth=2.,color='blue')
				cfp.plotvars.plot.plot(box_xpt["250"],box_ypt["250"],linewidth=2.,color='blue')
				cfp.gclose()
'''
Investigate AOD:
	1. With mm calculated from matched ensemble members
	2. hwi used here need to be from 1979-2100
	3. model standard deviation are calculated for signal-to-noise ratio
	4. Modify to fit single model
'''
if False:
	model="*"#"CanESM5"#

	result_historical={};result_scenario={}
	for scenario in scenarios[1:]:
		'''
		Present-day
			Find available source_ids and variant_labels from od550aer
		'''
		source_ids    =[x.split("_")[3] for x in glob.glob("/home/users/lguo/cssp_aero/data/cmip6/cmip/historical/od550aer/beijing/od550aer_AERmon_{}_historical_*_gn_185001-201412.nc".format(model))]
		variant_labels=[x.split("_")[5] for x in glob.glob("/home/users/lguo/cssp_aero/data/cmip6/cmip/historical/od550aer/beijing/od550aer_AERmon_{}_historical_*_gn_185001-201412.nc".format(model))]
		'''
		Present-day DJF mean
			Two periods are calculated: 1979-2014 & 2005-2014;
			Three means are calculated: all, HWI>1 and HWI<=0;
		'''
		od550aer_historical={}    ;od550aer_historical["long"]=[]    ;od550aer_historical["short"]=[];
		od550aer_historical_ge1={};od550aer_historical_ge1["long"]=[];od550aer_historical_ge1["short"]=[];
		od550aer_historical_le0={};od550aer_historical_le0["long"]=[];od550aer_historical_le0["short"]=[];
		od550aer_historical_ge2={};od550aer_historical_ge2["long"]=[];od550aer_historical_ge2["short"]=[];
		sources_n_variants=[] # used for matching historical with scenario
		for source_id,variant_label in zip(source_ids,variant_labels):
			try:
				print(scenario,"Historical",source_id,variant_label)
				finm="{}/cmip/historical/od550aer/beijing/od550aer_AERmon_{}_historical_{}_gn_185001-201412.nc".format(cmipdir,source_id,variant_label);print(finm)
				od550aer=cf.read(finm).select_field("atmosphere_optical_thickness_due_to_ambient_aerosol_particles")
				od550aer_long =od550aer.subspace(T=cf.wi(cf.dt("1979-01-01"),cf.dt("2014-12-30"))).subspace(T=cf.djf())
				od550aer_short=od550aer.subspace(T=cf.wi(cf.dt("2005-01-01"),cf.dt("2014-12-30"))).subspace(T=cf.djf())
				finm="{0}/ScenarioMIP/{1}/hwi/hwi_{2}_{1}_{3}_gn_197901-210012.nc".format(cmipdir,scenario,source_id,variant_label);print(finm)
				hwi=cf.read(finm).select_field("hwi")
				hwi_long =hwi.subspace(T=cf.wi(cf.dt("1979-01-01"),cf.dt("2014-12-30")))
				hwi_short=hwi.subspace(T=cf.wi(cf.dt("2005-01-01"),cf.dt("2014-12-30")))
				'''
				For climatology
				'''
				od550aer_historical["long"].append(od550aer_long.collapse("mean","T"))
				od550aer_historical["short"].append(od550aer_short.collapse("mean","T"))
				'''
				For severe haze months
				'''
				od550aer_long_ge1=od550aer_long.where(hwi_long<1,cf.masked)
				od550aer_historical_ge1["long"].append(od550aer_long_ge1.collapse("mean","T"))
				od550aer_short_ge1=od550aer_short.where(hwi_short<1,cf.masked)
				od550aer_historical_ge1["short"].append(od550aer_short_ge1.collapse("mean","T"))
				'''
				For non-haze months
				'''
				od550aer_long_le0=od550aer_long.where(hwi_long>0,cf.masked)
				od550aer_historical_le0["long"].append(od550aer_long_le0.collapse("mean","T"))
				od550aer_short_le0=od550aer_short.where(hwi_short>0,cf.masked)
				od550aer_historical_le0["short"].append(od550aer_short_le0.collapse("mean","T"))
				'''
				For severe haze months
				'''
				od550aer_long_ge2=od550aer_long.where(hwi_long<2,cf.masked)
				od550aer_historical_ge2["long"].append(od550aer_long_ge2.collapse("mean","T"))
				od550aer_short_ge2=od550aer_short.where(hwi_short<2,cf.masked)
				od550aer_historical_ge2["short"].append(od550aer_short_ge2.collapse("mean","T"))
				'''
				Store sourceid/variant_label
				'''
				sources_n_variants.append("/".join([source_id,variant_label]))
				
			except:
				print("{}_{} does not available in hwi ...".format(source_id,variant_label))
		'''
		Future scenarios
		'''
		source_ids    =[x.split("_")[3] for x in glob.glob("/home/users/lguo/cssp_aero/data/cmip6/ScenarioMIP/{0}/od550aer/beijing/od550aer_AERmon_{1}_{0}_*_gn_201501-210012.nc".format(scenario,model))]
		variant_labels=[x.split("_")[5] for x in glob.glob("/home/users/lguo/cssp_aero/data/cmip6/ScenarioMIP/{0}/od550aer/beijing/od550aer_AERmon_{1}_{0}_*_gn_201501-210012.nc".format(scenario,model))]
		'''
		Future DJF mean
			Four periods are calcualted: 30, 40 50 & 95;
			Three means are calculated: all, HWI>1 and HWI<=0;
		'''
		od550aer_scenario_historical={}    ;od550aer_scenario_historical["long"]=[]    ;od550aer_scenario_historical["short"]=[]
		od550aer_scenario_historical_ge1={};od550aer_scenario_historical_ge1["long"]=[];od550aer_scenario_historical_ge1["short"]=[]
		od550aer_scenario_historical_le0={};od550aer_scenario_historical_le0["long"]=[];od550aer_scenario_historical_le0["short"]=[]
		od550aer_scenario_historical_ge2={};od550aer_scenario_historical_ge2["long"]=[];od550aer_scenario_historical_ge2["short"]=[]

		od550aer_scenario={}    ;od550aer_scenario["30"]=[]    ;od550aer_scenario["40"]=[]    ;od550aer_scenario["50"]=[]    ;od550aer_scenario["95"]=[]
		od550aer_scenario_ge1={};od550aer_scenario_ge1["30"]=[];od550aer_scenario_ge1["40"]=[];od550aer_scenario_ge1["50"]=[];od550aer_scenario_ge1["95"]=[]
		od550aer_scenario_le0={};od550aer_scenario_le0["30"]=[];od550aer_scenario_le0["40"]=[];od550aer_scenario_le0["50"]=[];od550aer_scenario_le0["95"]=[]
		od550aer_scenario_ge2={};od550aer_scenario_ge2["30"]=[];od550aer_scenario_ge2["40"]=[];od550aer_scenario_ge2["50"]=[];od550aer_scenario_ge2["95"]=[]
		for source_id,variant_label in zip(source_ids,variant_labels):
			try:
				print(scenario,scenario,source_id,variant_label)
				'''
				Match to Present-day od550aer with availability in future scenarios
				'''
				try:
					index_historical=sources_n_variants.index("/".join([source_id,variant_label]))
					print(index_historical)
					od550aer_scenario_historical["long"].append(od550aer_historical["long"][index_historical])
					od550aer_scenario_historical["short"].append(od550aer_historical["short"][index_historical])
					od550aer_scenario_historical_ge1["long"].append(od550aer_historical_ge1["long"][index_historical])
					od550aer_scenario_historical_ge1["short"].append(od550aer_historical_ge1["short"][index_historical])
					od550aer_scenario_historical_le0["long"].append(od550aer_historical_le0["long"][index_historical])
					od550aer_scenario_historical_le0["short"].append(od550aer_historical_le0["short"][index_historical])
					od550aer_scenario_historical_ge2["long"].append(od550aer_historical_ge2["long"][index_historical])
					od550aer_scenario_historical_ge2["short"].append(od550aer_historical_ge2["short"][index_historical])
				except:
					print("... Does not exist in Historical ...")
					continue

				print("... Continue processing Scenario ...")
				finm="{0}/ScenarioMIP/{1}/od550aer/beijing/od550aer_AERmon_{2}_{1}_{3}_gn_201501-210012.nc".format(cmipdir,scenario,source_id,variant_label);print(finm)
				od550aer=cf.read(finm).select_field("atmosphere_optical_thickness_due_to_ambient_aerosol_particles").subspace(T=cf.djf())
				finm="{0}/ScenarioMIP/{1}/hwi/hwi_{2}_{1}_{3}_gn_197901-210012.nc".format(cmipdir,scenario,source_id,variant_label);print(finm)
				hwi=cf.read(finm).select_field("hwi").subspace(T=cf.wi(cf.dt("2015-01-01"),cf.dt("2100-12-30")))
				'''
				For climatology
				'''
				print("-----> {} {} <-----".format(od550aer.get_property("source_id"),od550aer.get_property("variant_label")))
				od550aer_scenario["30"].append(od550aer.subspace(T=cf.wi(cf.dt("2025-01-01"),cf.dt("2034-12-30"))).collapse("mean","T"))
				od550aer_scenario["40"].append(od550aer.subspace(T=cf.wi(cf.dt("2035-01-01"),cf.dt("2044-12-30"))).collapse("mean","T"))
				od550aer_scenario["50"].append(od550aer.subspace(T=cf.wi(cf.dt("2045-01-01"),cf.dt("2054-12-30"))).collapse("mean","T"))
				od550aer_scenario["95"].append(od550aer.subspace(T=cf.wi(cf.dt("2090-01-01"),cf.dt("2099-12-30"))).collapse("mean","T"))
				'''
				For severe haze months
				'''
				od550aer_ge1=od550aer.where(hwi<1,cf.masked)
				od550aer_scenario_ge1["30"].append(od550aer_ge1.subspace(T=cf.wi(cf.dt("2025-01-01"),cf.dt("2034-12-30"))).collapse("mean","T"))
				od550aer_scenario_ge1["40"].append(od550aer_ge1.subspace(T=cf.wi(cf.dt("2035-01-01"),cf.dt("2044-12-30"))).collapse("mean","T"))
				od550aer_scenario_ge1["50"].append(od550aer_ge1.subspace(T=cf.wi(cf.dt("2045-01-01"),cf.dt("2054-12-30"))).collapse("mean","T"))
				od550aer_scenario_ge1["95"].append(od550aer_ge1.subspace(T=cf.wi(cf.dt("2090-01-01"),cf.dt("2099-12-30"))).collapse("mean","T"))
				'''
				For non-haze months
				'''
				od550aer_le0=od550aer.where(hwi>0,cf.masked)
				od550aer_scenario_le0["30"].append(od550aer_le0.subspace(T=cf.wi(cf.dt("2025-01-01"),cf.dt("2034-12-30"))).collapse("mean","T"))
				od550aer_scenario_le0["40"].append(od550aer_le0.subspace(T=cf.wi(cf.dt("2035-01-01"),cf.dt("2044-12-30"))).collapse("mean","T"))
				od550aer_scenario_le0["50"].append(od550aer_le0.subspace(T=cf.wi(cf.dt("2045-01-01"),cf.dt("2054-12-30"))).collapse("mean","T"))
				od550aer_scenario_le0["95"].append(od550aer_le0.subspace(T=cf.wi(cf.dt("2090-01-01"),cf.dt("2099-12-30"))).collapse("mean","T"))
				'''
				For severe haze months
				'''
				od550aer_ge2=od550aer.where(hwi<2,cf.masked)
				od550aer_scenario_ge2["30"].append(od550aer_ge2.subspace(T=cf.wi(cf.dt("2025-01-01"),cf.dt("2034-12-30"))).collapse("mean","T"))
				od550aer_scenario_ge2["40"].append(od550aer_ge2.subspace(T=cf.wi(cf.dt("2035-01-01"),cf.dt("2044-12-30"))).collapse("mean","T"))
				od550aer_scenario_ge2["50"].append(od550aer_ge2.subspace(T=cf.wi(cf.dt("2045-01-01"),cf.dt("2054-12-30"))).collapse("mean","T"))
				od550aer_scenario_ge2["95"].append(od550aer_ge2.subspace(T=cf.wi(cf.dt("2090-01-01"),cf.dt("2099-12-30"))).collapse("mean","T"))
			except FileNotFoundError:
				print("{}_{}_{} does not exist in od550aer ...".format(scenario,source_id,variant_label))
		'''
		Calculate model mean for aforementioned categories
		'''
		if model=="*":
			od550aer_scenario_historical["long"]=loc_func.fldlst2lstmm_new(od550aer_scenario_historical["long"])
			od550aer_scenario_historical["short"]=loc_func.fldlst2lstmm_new(od550aer_scenario_historical["short"])
			od550aer_scenario_historical_ge1["long"]=loc_func.fldlst2lstmm_new(od550aer_scenario_historical_ge1["long"])
			od550aer_scenario_historical_ge1["short"]=loc_func.fldlst2lstmm_new(od550aer_scenario_historical_ge1["short"])
			od550aer_scenario_historical_le0["long"]=loc_func.fldlst2lstmm_new(od550aer_scenario_historical_le0["long"])
			od550aer_scenario_historical_le0["short"]=loc_func.fldlst2lstmm_new(od550aer_scenario_historical_le0["short"])
			od550aer_scenario_historical_ge2["long"]=loc_func.fldlst2lstmm_new(od550aer_scenario_historical_ge2["long"])
			od550aer_scenario_historical_ge2["short"]=loc_func.fldlst2lstmm_new(od550aer_scenario_historical_ge2["short"])
	
			od550aer_scenario["30"]=loc_func.fldlst2lstmm_new(od550aer_scenario["30"])
			od550aer_scenario["40"]=loc_func.fldlst2lstmm_new(od550aer_scenario["40"])
			od550aer_scenario["50"]=loc_func.fldlst2lstmm_new(od550aer_scenario["50"])
			od550aer_scenario["95"]=loc_func.fldlst2lstmm_new(od550aer_scenario["95"])
	
			od550aer_scenario_ge1["30"]=loc_func.fldlst2lstmm_new(od550aer_scenario_ge1["30"])
			od550aer_scenario_ge1["40"]=loc_func.fldlst2lstmm_new(od550aer_scenario_ge1["40"])
			od550aer_scenario_ge1["50"]=loc_func.fldlst2lstmm_new(od550aer_scenario_ge1["50"])
			od550aer_scenario_ge1["95"]=loc_func.fldlst2lstmm_new(od550aer_scenario_ge1["95"])
	
			od550aer_scenario_le0["30"]=loc_func.fldlst2lstmm_new(od550aer_scenario_le0["30"])
			od550aer_scenario_le0["40"]=loc_func.fldlst2lstmm_new(od550aer_scenario_le0["40"])
			od550aer_scenario_le0["50"]=loc_func.fldlst2lstmm_new(od550aer_scenario_le0["50"])
			od550aer_scenario_le0["95"]=loc_func.fldlst2lstmm_new(od550aer_scenario_le0["95"])
	
			od550aer_scenario_ge2["30"]=loc_func.fldlst2lstmm_new(od550aer_scenario_ge2["30"])
			od550aer_scenario_ge2["40"]=loc_func.fldlst2lstmm_new(od550aer_scenario_ge2["40"])
			od550aer_scenario_ge2["50"]=loc_func.fldlst2lstmm_new(od550aer_scenario_ge2["50"])
			od550aer_scenario_ge2["95"]=loc_func.fldlst2lstmm_new(od550aer_scenario_ge2["95"])
		else:
			od550aer_scenario_historical["long"]=od550aer_scenario_historical["long"]
			od550aer_scenario_historical["short"]=od550aer_scenario_historical["short"]
			od550aer_scenario_historical_ge1["long"]=od550aer_scenario_historical_ge1["long"]
			od550aer_scenario_historical_ge1["short"]=od550aer_scenario_historical_ge1["short"]
			od550aer_scenario_historical_le0["long"]=od550aer_scenario_historical_le0["long"]
			od550aer_scenario_historical_le0["short"]=od550aer_scenario_historical_le0["short"]
			od550aer_scenario_historical_ge2["long"]=od550aer_scenario_historical_ge2["long"]
			od550aer_scenario_historical_ge2["short"]=od550aer_scenario_historical_ge2["short"]
	
			od550aer_scenario["30"]=od550aer_scenario["30"]
			od550aer_scenario["40"]=od550aer_scenario["40"]
			od550aer_scenario["50"]=od550aer_scenario["50"]
			od550aer_scenario["95"]=od550aer_scenario["95"]
	
			od550aer_scenario_ge1["30"]=od550aer_scenario_ge1["30"]
			od550aer_scenario_ge1["40"]=od550aer_scenario_ge1["40"]
			od550aer_scenario_ge1["50"]=od550aer_scenario_ge1["50"]
			od550aer_scenario_ge1["95"]=od550aer_scenario_ge1["95"]
	
			od550aer_scenario_le0["30"]=od550aer_scenario_le0["30"]
			od550aer_scenario_le0["40"]=od550aer_scenario_le0["40"]
			od550aer_scenario_le0["50"]=od550aer_scenario_le0["50"]
			od550aer_scenario_le0["95"]=od550aer_scenario_le0["95"]
	
			od550aer_scenario_ge2["30"]=od550aer_scenario_ge2["30"]
			od550aer_scenario_ge2["40"]=od550aer_scenario_ge2["40"]
			od550aer_scenario_ge2["50"]=od550aer_scenario_ge2["50"]
			od550aer_scenario_ge2["95"]=od550aer_scenario_ge2["95"]
		'''
		Calculate multi-model mean and standara deviation
		'''
		result_historical[scenario]={}
		result_historical[scenario]["cli"]={}
		result_historical[scenario]["cli"]["long"]=[loc_func.mean_fldlst(od550aer_scenario_historical["long"])]
		result_historical[scenario]["cli"]["long"].append(loc_func.median_fldlst(od550aer_scenario_historical["long"]))
		result_historical[scenario]["cli"]["long"].append(loc_func.stdev_fldlst(od550aer_scenario_historical["long"]))
		result_historical[scenario]["cli"]["long"].append(loc_func.quartile_fldlst(od550aer_scenario_historical["long"]))
		result_historical[scenario]["cli"]["long"].append(od550aer_scenario_historical["long"])
		result_historical[scenario]["cli"]["long"].append(loc_func.tolist_fldlst(od550aer_scenario_historical["long"]))
		result_historical[scenario]["cli"]["short"]=[loc_func.mean_fldlst(od550aer_scenario_historical["short"])]
		result_historical[scenario]["cli"]["short"]=[loc_func.median_fldlst(od550aer_scenario_historical["short"])]
		result_historical[scenario]["cli"]["short"].append(loc_func.stdev_fldlst(od550aer_scenario_historical["short"]))
		result_historical[scenario]["cli"]["short"].append(loc_func.quartile_fldlst(od550aer_scenario_historical["short"]))
		result_historical[scenario]["cli"]["short"].append(loc_func.tolist_fldlst(od550aer_scenario_historical["short"]))
		result_historical[scenario]["cli"]["short"].append(od550aer_scenario_historical["short"])
		result_historical[scenario]["ge1"]={}
		result_historical[scenario]["ge1"]["long"]=[loc_func.mean_fldlst(od550aer_scenario_historical_ge1["long"])]
		result_historical[scenario]["ge1"]["long"].append(loc_func.median_fldlst(od550aer_scenario_historical_ge1["long"]))
		result_historical[scenario]["ge1"]["long"].append(loc_func.stdev_fldlst(od550aer_scenario_historical_ge1["long"]))
		result_historical[scenario]["ge1"]["long"].append(loc_func.quartile_fldlst(od550aer_scenario_historical_ge1["long"]))
		result_historical[scenario]["ge1"]["long"].append(loc_func.tolist_fldlst(od550aer_scenario_historical_ge1["long"]))
		result_historical[scenario]["ge1"]["long"].append(od550aer_scenario_historical_ge1["long"])
		result_historical[scenario]["ge1"]["short"]=[loc_func.mean_fldlst(od550aer_scenario_historical_ge1["short"])]
		result_historical[scenario]["ge1"]["short"].append(loc_func.median_fldlst(od550aer_scenario_historical_ge1["short"]))
		result_historical[scenario]["ge1"]["short"].append(loc_func.stdev_fldlst(od550aer_scenario_historical_ge1["short"]))
		result_historical[scenario]["ge1"]["short"].append(loc_func.quartile_fldlst(od550aer_scenario_historical_ge1["short"]))
		result_historical[scenario]["ge1"]["short"].append(loc_func.tolist_fldlst(od550aer_scenario_historical_ge1["short"]))
		result_historical[scenario]["ge1"]["short"].append(od550aer_scenario_historical_ge1["short"])
		result_historical[scenario]["le0"]={}
		result_historical[scenario]["le0"]["long"] =[loc_func.mean_fldlst(od550aer_scenario_historical_le0["long"])]
		result_historical[scenario]["le0"]["long"].append(loc_func.median_fldlst(od550aer_scenario_historical_le0["long"]))
		result_historical[scenario]["le0"]["long"].append(loc_func.stdev_fldlst(od550aer_scenario_historical_le0["long"]))
		result_historical[scenario]["le0"]["long"].append(loc_func.quartile_fldlst(od550aer_scenario_historical_le0["long"]))
		result_historical[scenario]["le0"]["long"].append(loc_func.tolist_fldlst(od550aer_scenario_historical_le0["long"]))
		result_historical[scenario]["le0"]["long"].append(od550aer_scenario_historical_le0["long"])
		result_historical[scenario]["le0"]["short"]=[loc_func.mean_fldlst(od550aer_scenario_historical_le0["short"])]
		result_historical[scenario]["le0"]["short"].append(loc_func.median_fldlst(od550aer_scenario_historical_le0["short"]))
		result_historical[scenario]["le0"]["short"].append(loc_func.stdev_fldlst(od550aer_scenario_historical_le0["short"]))
		result_historical[scenario]["le0"]["short"].append(loc_func.quartile_fldlst(od550aer_scenario_historical_le0["short"]))
		result_historical[scenario]["le0"]["short"].append(loc_func.tolist_fldlst(od550aer_scenario_historical_le0["short"]))
		result_historical[scenario]["le0"]["short"].append(od550aer_scenario_historical_le0["short"])
		result_historical[scenario]["ge2"]={}
		result_historical[scenario]["ge2"]["long"]=[loc_func.mean_fldlst(od550aer_scenario_historical_ge2["long"])]
		result_historical[scenario]["ge2"]["long"].append(loc_func.median_fldlst(od550aer_scenario_historical_ge2["long"]))
		result_historical[scenario]["ge2"]["long"].append(loc_func.stdev_fldlst(od550aer_scenario_historical_ge2["long"]))
		result_historical[scenario]["ge2"]["long"].append(loc_func.quartile_fldlst(od550aer_scenario_historical_ge2["long"]))
		result_historical[scenario]["ge2"]["long"].append(loc_func.tolist_fldlst(od550aer_scenario_historical_ge2["long"]))
		result_historical[scenario]["ge2"]["long"].append(od550aer_scenario_historical_ge2["long"])
		result_historical[scenario]["ge2"]["short"]=[loc_func.mean_fldlst(od550aer_scenario_historical_ge2["short"])]
		result_historical[scenario]["ge2"]["short"].append(loc_func.median_fldlst(od550aer_scenario_historical_ge2["short"]))
		result_historical[scenario]["ge2"]["short"].append(loc_func.stdev_fldlst(od550aer_scenario_historical_ge2["short"]))
		result_historical[scenario]["ge2"]["short"].append(loc_func.quartile_fldlst(od550aer_scenario_historical_ge2["short"]))
		result_historical[scenario]["ge2"]["short"].append(loc_func.tolist_fldlst(od550aer_scenario_historical_ge2["short"]))
		result_historical[scenario]["ge2"]["short"].append(od550aer_scenario_historical_ge2["short"])

		result_scenario[scenario]={}
		result_scenario[scenario]["30"]={};
		result_scenario[scenario]["30"]["cli"]=[loc_func.mean_fldlst(od550aer_scenario["30"])]
		result_scenario[scenario]["30"]["cli"].append(loc_func.median_fldlst(od550aer_scenario["30"]))
		result_scenario[scenario]["30"]["cli"].append(loc_func.stdev_fldlst(od550aer_scenario["30"]))
		result_scenario[scenario]["30"]["cli"].append(loc_func.quartile_fldlst(od550aer_scenario["30"]))
		result_scenario[scenario]["30"]["cli"].append(loc_func.tolist_fldlst(od550aer_scenario["30"]))
		result_scenario[scenario]["30"]["cli"].append(od550aer_scenario["30"])
		result_scenario[scenario]["30"]["ge1"]=[loc_func.mean_fldlst(od550aer_scenario_ge1["30"])]
		result_scenario[scenario]["30"]["ge1"].append(loc_func.median_fldlst(od550aer_scenario_ge1["30"]))
		result_scenario[scenario]["30"]["ge1"].append(loc_func.stdev_fldlst(od550aer_scenario_ge1["30"]))
		result_scenario[scenario]["30"]["ge1"].append(loc_func.quartile_fldlst(od550aer_scenario_ge1["30"]))
		result_scenario[scenario]["30"]["ge1"].append(loc_func.tolist_fldlst(od550aer_scenario_ge1["30"]))
		result_scenario[scenario]["30"]["ge1"].append(od550aer_scenario_ge1["30"])
		result_scenario[scenario]["30"]["le0"]=[loc_func.mean_fldlst(od550aer_scenario_le0["30"])]
		result_scenario[scenario]["30"]["le0"].append(loc_func.median_fldlst(od550aer_scenario_le0["30"]))
		result_scenario[scenario]["30"]["le0"].append(loc_func.stdev_fldlst(od550aer_scenario_le0["30"]))
		result_scenario[scenario]["30"]["le0"].append(loc_func.quartile_fldlst(od550aer_scenario_le0["30"]))
		result_scenario[scenario]["30"]["le0"].append(loc_func.tolist_fldlst(od550aer_scenario_le0["30"]))
		result_scenario[scenario]["30"]["le0"].append(od550aer_scenario_le0["30"])
		result_scenario[scenario]["30"]["ge2"]=[loc_func.mean_fldlst(od550aer_scenario_ge2["30"])]
		result_scenario[scenario]["30"]["ge2"].append(loc_func.median_fldlst(od550aer_scenario_ge2["30"]))
		result_scenario[scenario]["30"]["ge2"].append(loc_func.stdev_fldlst(od550aer_scenario_ge2["30"]))
		result_scenario[scenario]["30"]["ge2"].append(loc_func.quartile_fldlst(od550aer_scenario_ge2["30"]))
		result_scenario[scenario]["30"]["ge2"].append(loc_func.tolist_fldlst(od550aer_scenario_ge2["30"]))
		result_scenario[scenario]["30"]["ge2"].append(od550aer_scenario_ge2["30"])
		result_scenario[scenario]["40"]={};
		result_scenario[scenario]["40"]["cli"]=[loc_func.mean_fldlst(od550aer_scenario["40"])]
		result_scenario[scenario]["40"]["cli"].append(loc_func.median_fldlst(od550aer_scenario["40"]))
		result_scenario[scenario]["40"]["cli"].append(loc_func.stdev_fldlst(od550aer_scenario["40"]))
		result_scenario[scenario]["40"]["cli"].append(loc_func.quartile_fldlst(od550aer_scenario["40"]))
		result_scenario[scenario]["40"]["cli"].append(loc_func.tolist_fldlst(od550aer_scenario["40"]))
		result_scenario[scenario]["40"]["cli"].append(od550aer_scenario["40"])
		result_scenario[scenario]["40"]["ge1"]=[loc_func.mean_fldlst(od550aer_scenario_ge1["40"])]
		result_scenario[scenario]["40"]["ge1"].append(loc_func.median_fldlst(od550aer_scenario_ge1["40"]))
		result_scenario[scenario]["40"]["ge1"].append(loc_func.stdev_fldlst(od550aer_scenario_ge1["40"]))
		result_scenario[scenario]["40"]["ge1"].append(loc_func.quartile_fldlst(od550aer_scenario_ge1["40"]))
		result_scenario[scenario]["40"]["ge1"].append(loc_func.tolist_fldlst(od550aer_scenario_ge1["40"]))
		result_scenario[scenario]["40"]["ge1"].append(od550aer_scenario_ge1["40"])
		result_scenario[scenario]["40"]["le0"]=[loc_func.mean_fldlst(od550aer_scenario_le0["40"])]
		result_scenario[scenario]["40"]["le0"].append(loc_func.median_fldlst(od550aer_scenario_le0["40"]))
		result_scenario[scenario]["40"]["le0"].append(loc_func.stdev_fldlst(od550aer_scenario_le0["40"]))
		result_scenario[scenario]["40"]["le0"].append(loc_func.quartile_fldlst(od550aer_scenario_le0["40"]))
		result_scenario[scenario]["40"]["le0"].append(loc_func.tolist_fldlst(od550aer_scenario_le0["40"]))
		result_scenario[scenario]["40"]["le0"].append(od550aer_scenario_le0["40"])
		result_scenario[scenario]["40"]["ge2"]=[loc_func.mean_fldlst(od550aer_scenario_ge2["40"])]
		result_scenario[scenario]["40"]["ge2"].append(loc_func.median_fldlst(od550aer_scenario_ge2["40"]))
		result_scenario[scenario]["40"]["ge2"].append(loc_func.stdev_fldlst(od550aer_scenario_ge2["40"]))
		result_scenario[scenario]["40"]["ge2"].append(loc_func.quartile_fldlst(od550aer_scenario_ge2["40"]))
		result_scenario[scenario]["40"]["ge2"].append(loc_func.tolist_fldlst(od550aer_scenario_ge2["40"]))
		result_scenario[scenario]["40"]["ge2"].append(od550aer_scenario_ge2["40"])
		result_scenario[scenario]["50"]={};
		result_scenario[scenario]["50"]["cli"]=[loc_func.mean_fldlst(od550aer_scenario["50"])]
		result_scenario[scenario]["50"]["cli"].append(loc_func.median_fldlst(od550aer_scenario["50"]))
		result_scenario[scenario]["50"]["cli"].append(loc_func.stdev_fldlst(od550aer_scenario["50"]))
		result_scenario[scenario]["50"]["cli"].append(loc_func.quartile_fldlst(od550aer_scenario["50"]))
		result_scenario[scenario]["50"]["cli"].append(loc_func.tolist_fldlst(od550aer_scenario["50"]))
		result_scenario[scenario]["50"]["cli"].append(od550aer_scenario["50"])
		result_scenario[scenario]["50"]["ge1"]=[loc_func.mean_fldlst(od550aer_scenario_ge1["50"])]
		result_scenario[scenario]["50"]["ge1"].append(loc_func.median_fldlst(od550aer_scenario_ge1["50"]))
		result_scenario[scenario]["50"]["ge1"].append(loc_func.stdev_fldlst(od550aer_scenario_ge1["50"]))
		result_scenario[scenario]["50"]["ge1"].append(loc_func.quartile_fldlst(od550aer_scenario_ge1["50"]))
		result_scenario[scenario]["50"]["ge1"].append(loc_func.tolist_fldlst(od550aer_scenario_ge1["50"]))
		result_scenario[scenario]["50"]["ge1"].append(od550aer_scenario_ge1["50"])
		result_scenario[scenario]["50"]["le0"]=[loc_func.mean_fldlst(od550aer_scenario_le0["50"])]
		result_scenario[scenario]["50"]["le0"].append(loc_func.median_fldlst(od550aer_scenario_le0["50"]))
		result_scenario[scenario]["50"]["le0"].append(loc_func.stdev_fldlst(od550aer_scenario_le0["50"]))
		result_scenario[scenario]["50"]["le0"].append(loc_func.quartile_fldlst(od550aer_scenario_le0["50"]))
		result_scenario[scenario]["50"]["le0"].append(loc_func.tolist_fldlst(od550aer_scenario_le0["50"]))
		result_scenario[scenario]["50"]["le0"].append(od550aer_scenario_le0["50"])
		result_scenario[scenario]["50"]["ge2"]=[loc_func.mean_fldlst(od550aer_scenario_ge2["50"])]
		result_scenario[scenario]["50"]["ge2"].append(loc_func.median_fldlst(od550aer_scenario_ge2["50"]))
		result_scenario[scenario]["50"]["ge2"].append(loc_func.stdev_fldlst(od550aer_scenario_ge2["50"]))
		result_scenario[scenario]["50"]["ge2"].append(loc_func.quartile_fldlst(od550aer_scenario_ge2["50"]))
		result_scenario[scenario]["50"]["ge2"].append(loc_func.tolist_fldlst(od550aer_scenario_ge2["50"]))
		result_scenario[scenario]["50"]["ge2"].append(od550aer_scenario_ge2["50"])
		result_scenario[scenario]["95"]={};
		result_scenario[scenario]["95"]["cli"]=[loc_func.mean_fldlst(od550aer_scenario["95"])]
		result_scenario[scenario]["95"]["cli"].append(loc_func.median_fldlst(od550aer_scenario["95"]))
		result_scenario[scenario]["95"]["cli"].append(loc_func.stdev_fldlst(od550aer_scenario["95"]))
		result_scenario[scenario]["95"]["cli"].append(loc_func.quartile_fldlst(od550aer_scenario["95"]))
		result_scenario[scenario]["95"]["cli"].append(loc_func.tolist_fldlst(od550aer_scenario["95"]))
		result_scenario[scenario]["95"]["cli"].append(od550aer_scenario["95"])
		result_scenario[scenario]["95"]["ge1"]=[loc_func.mean_fldlst(od550aer_scenario_ge1["95"])]
		result_scenario[scenario]["95"]["ge1"].append(loc_func.median_fldlst(od550aer_scenario_ge1["95"]))
		result_scenario[scenario]["95"]["ge1"].append(loc_func.stdev_fldlst(od550aer_scenario_ge1["95"]))
		result_scenario[scenario]["95"]["ge1"].append(loc_func.quartile_fldlst(od550aer_scenario_ge1["95"]))
		result_scenario[scenario]["95"]["ge1"].append(loc_func.tolist_fldlst(od550aer_scenario_ge1["95"]))
		result_scenario[scenario]["95"]["ge1"].append(od550aer_scenario_ge1["95"])
		result_scenario[scenario]["95"]["le0"]=[loc_func.mean_fldlst(od550aer_scenario_le0["95"])]
		result_scenario[scenario]["95"]["le0"].append(loc_func.median_fldlst(od550aer_scenario_le0["95"]))
		result_scenario[scenario]["95"]["le0"].append(loc_func.stdev_fldlst(od550aer_scenario_le0["95"]))
		result_scenario[scenario]["95"]["le0"].append(loc_func.quartile_fldlst(od550aer_scenario_le0["95"]))
		result_scenario[scenario]["95"]["le0"].append(loc_func.tolist_fldlst(od550aer_scenario_le0["95"]))
		result_scenario[scenario]["95"]["le0"].append(od550aer_scenario_le0["95"])
		result_scenario[scenario]["95"]["ge2"]=[loc_func.mean_fldlst(od550aer_scenario_ge2["95"])]
		result_scenario[scenario]["95"]["ge2"].append(loc_func.median_fldlst(od550aer_scenario_ge2["95"]))
		result_scenario[scenario]["95"]["ge2"].append(loc_func.stdev_fldlst(od550aer_scenario_ge2["95"]))
		result_scenario[scenario]["95"]["ge2"].append(loc_func.quartile_fldlst(od550aer_scenario_ge2["95"]))
		result_scenario[scenario]["95"]["ge2"].append(loc_func.tolist_fldlst(od550aer_scenario_ge2["95"]))
		result_scenario[scenario]["95"]["ge2"].append(od550aer_scenario_ge2["95"])
	'''
	Write out
	'''
	if model=="*": fonm="./data/od550aer_ano.pickle"
	else: fonm="./data/od550aer_ano_{}.pickle".format(model)
	with open(fonm,"wb") as fp:
		pickle.dump(result_historical,fp)
		pickle.dump(result_scenario,fp)
'''
Summarise AOD anomalies
	Illustrated into seperated plots
'''
if False:
	hist_bs="ssp370"
	comp_bs="le0"
	period_bs="long"
	mean="mean"
	model=""#"CanESM5"#

	scenarios=scenarios[1:]
	boxcolors=['#1E9684','#1D3354','#EADD3D','#F21111','#840B22']
	boxcolors=boxcolors[1:]
	legend_labels=["Hist.","SSP1-2.6","SSP2-4.5","SSP3-7.0","SSP5-8.5"]
	'''
	Read
	'''
	fp = open("./data/od550aer_ano.pickle","rb") if model=="" else open("./data/od550aer_ano_{}.pickle".format(model),"rb")	
	aod_historical=pickle.load(fp)
	aod_scenario=pickle.load(fp)
	'''
	Multi-information for f-test
	'''
	if model=="CanESM5": 
		fp=open("./data/od550aer_ano.pickle","rb")
		mm_historical=pickle.load(fp)
		mm_scenario=pickle.load(fp)
	'''
	Plot AOD anomaly against PD - relative changes
	'''
	if False:
		diff_scenario={};sig_scenario={}
		if model=="":
			for scenario in scenarios:
				source_ids=[]
				for fld in aod_historical[scenario]["ge1"][period_bs][5]:
					source_ids.append(fld.get_property("source_id"))

				diff_scenario[scenario]=[];sig_scenario[scenario]=[]
				for period in ["30","40","50","95"]:
					print(scenario,period)
					tmp_lst=[]
					for fld in aod_scenario[scenario][period]["ge1"][5]:
						source_id=fld.get_property("source_id")
						idx=source_ids.index(source_id)
						tmp_lst.append(((fld-aod_historical[scenario]["ge1"][period_bs][5][idx])/aod_historical[scenario]["ge1"][period_bs][5][idx]*100.).array[0])
					diff_scenario[scenario].append(np.array(tmp_lst))
					'''
					T-test: test the significant difference between means 
					'''
					sig_scenario[scenario].append(loc_func.ttest4mean(aod_scenario[scenario][period]["ge1"][0],aod_historical[scenario][comp_bs][period_bs][0],aod_scenario[scenario][period]["ge1"][2],aod_historical[scenario][comp_bs][period_bs][2],12,12))
		if model=="CanESM5":
			for scenario in scenarios:
				diff_scenario[scenario]=[];sig_scenario[scenario]=[];
				source_ids=[]
				for fld in aod_historical[scenario]["ge1"][period_bs][5]:
					source_ids.append(fld.get_property("variant_label"))
				for period in ["30","40","50","95"]:
					print(scenario,period)
					tmp_lst=[]
					for fld in aod_scenario[scenario][period]["ge1"][5]:
						source_id=fld.get_property("variant_label")
						idx=source_ids.index(source_id)
						if fld.mask.array[0]:
							print(source_id,idx,fld.array[0])
						else:
							tmp_lst.append(((fld-aod_historical[scenario]["ge1"][period_bs][5][idx])/aod_historical[scenario]["ge1"][period_bs][5][idx]*100.).array[0])
					diff_scenario[scenario].append(np.array(tmp_lst))
					'''
					F-test: test the significant difference between variances
					'''
					sig_scenario[scenario].append(loc_func.ftest(aod_scenario[scenario][period]["ge1"][2],mm_scenario[scenario][period]["ge1"][2],10,10))
		box={}
		for i,scenario in enumerate(scenarios):
			medianprops=dict(linewidth=2,color='k')
			flierprops=dict(color=boxcolors[i],markerfacecolor=boxcolors[i],marker="o",markersize=4,alpha=.4)
			capprops=dict(linestyle="-",linewidth=1,color="k")
			pos=[.7+i*0.2, 2.6+i*0.2, 4.6+i*0.2, 6.6+i*0.2]
			for j,sig in enumerate(sig_scenario[scenario]):
				boxprops = dict(color=boxcolors[i],facecolor="w",linewidth=2,alpha=.4) if sig>.1 else dict(color=boxcolors[i],facecolor=boxcolors[i],linewidth=2,alpha=.4)
				plt.boxplot(diff_scenario[scenario][j],whis=(5,95),widths=.1,positions=[pos[j]],showcaps=True,showfliers=False,patch_artist=True,medianprops=medianprops,boxprops=boxprops,capprops=capprops)
		plt.ylim(-99,99)
		plt.xlim(0,8.)
		plt.plot([0,8.],[0,0],c='k',lw=1,ls=':')
		plt.xticks([1,3,5,7],['2025-2034','2035-2044','2045-2054','2090-2099'])
		plt.ylabel('AOD@550nm Anomaly (%)',fontsize=15)
		plt.gca().tick_params(labelsize=14,direction="out",bottom=True,left=True)
	
		if model=="": plt.savefig('./plot/hwi_scenariomip/fig_aod_scenario_pd_anomaly_{}.pdf'.format(period_bs),bbox_inches='tight',pad_inches=0)
		else: plt.savefig('./plot/hwi_scenariomip/fig_aod_scenario_pd_anomaly_{}_{}.pdf'.format(period_bs,model),bbox_inches='tight',pad_inches=0)
		plt.show()
		plt.close()
	'''
	Plot AOD anomaly against PD - relative changes without spread
	'''
	if False:
		diff_scenario={};sig_scenario={};
		if model=="":
			for scenario in scenarios:
				source_ids=[]
				for fld in aod_historical[scenario]["ge1"][period_bs][5]:
					source_ids.append(fld.get_property("source_id"))

				diff_scenario[scenario]=[];sig_scenario[scenario]=[]
				for period in ["30","40","50","95"]:
					tmp_lst=[]
					for fld in aod_scenario[scenario][period]["ge1"][5]:
						source_id=fld.get_property("source_id")
						idx=source_ids.index(source_id)
						tmp_lst.append(((fld-aod_historical[scenario]["ge1"][period_bs][5][idx])/aod_historical[scenario]["ge1"][period_bs][5][idx]*100.).array[0])
					diff_scenario[scenario].append(np.nanmean(np.array(tmp_lst)))
					'''
					T-test: test the significant difference between means 
					'''
					#sig_scenario[scenario].append(loc_func.ttest4mean(aod_scenario[scenario][period]["ge1"][0],aod_historical[scenario][comp_bs][period_bs][0],aod_scenario[scenario][period]["ge1"][2],aod_historical[scenario][comp_bs][period_bs][2],12,12))
					'''
					T-test: whether AOD changes in scenarios are significantly different from present-day.
						This is one-sample test against zero change.
					'''
					#sig_scenario[scenario].append(loc_func.ttest4mean_1sample(np.nanmean(tmp_lst),0,np.nanstd(tmp_lst),len(tmp_lst)))
		if model=="CanESM5":
			for scenario in scenarios:
				source_ids=[]
				for fld in aod_historical[scenario]["ge1"][period_bs][5]:
					source_ids.append(fld.get_property("variant_label"))

				diff_scenario[scenario]=[];sig_scenario[scenario]=[];
				for period in ["30","40","50","95"]:
					print(scenario,period)
					tmp_lst=[]
					for fld in aod_scenario[scenario][period]["ge1"][5]:
						source_id=fld.get_property("variant_label")
						idx=source_ids.index(source_id)
						if fld.mask.array[0]:
							print(source_id,idx,fld.array[0])
						else:
							tmp_lst.append(((fld-aod_historical[scenario]["ge1"][period_bs][5][idx])/aod_historical[scenario]["ge1"][period_bs][5][idx]*100.).array[0])
					diff_scenario[scenario].append(np.nanmean(np.array(tmp_lst)))
					'''
					F-test: test the significant difference between variances
					'''
					sig_scenario[scenario].append(loc_func.ftest(aod_scenario[scenario][period]["ge1"][2],mm_scenario[scenario][period]["ge1"][2],20,20))
		box={}
		pos=[1,3,5,7]
		for i,scenario in enumerate(scenarios):
			box[scenario]=plt.plot(pos,diff_scenario[scenario],marker='s',color=boxcolors[i],markersize=6,linestyle=":",linewidth=2,alpha=.8,markerfacecolor="none",markeredgewidth=2)
			if model=="CanESM5":
				for j,sig in enumerate(sig_scenario[scenario]):
					if sig<=.1:
						print(scenario,j,sig)
						plt.plot(pos[j],diff_scenario[scenario][j],marker="s",color=boxcolors[i],markersize=6,alpha=.8,markerfacecolor=boxcolors[i],markeredgewidth=0)
		plt.legend((box[scenario][0] for scenario in scenarios),legend_labels[1:],loc='upper right',fontsize=15,ncol=2,framealpha=.4)

		plt.xlim(0,8.);plt.ylim(-99,99)
		plt.plot([0,8.],[0,0],c='k',lw=1,ls=':')
		plt.xticks([1,3,5,7],['2025-2034','2035-2044','2045-2054','2090-2099'])
		plt.ylabel('AOD@550nm Anomaly (%)',fontsize=15)
		plt.gca().tick_params(labelsize=14,direction="out",bottom=True,left=True)
	
		plt.savefig('./plot/hwi_scenariomip/fig_aod_scenario_pd_anomaly_{}.pdf'.format(period_bs),bbox_inches='tight',pad_inches=0) if model=="" else plt.savefig('./plot/hwi_scenariomip/fig_aod_scenario_pd_anomaly_{}_{}.pdf'.format(period_bs,model),bbox_inches='tight',pad_inches=0)
		plt.show()
		plt.close()
	'''
	Plot AOD anomaly against comteporary - relative changes without spread
	'''
	if True:
		diff_historical=[];
		for i,fld in enumerate(aod_historical[hist_bs]["ge2"][period_bs][5]):
			diff_historical.append(((fld-aod_historical[hist_bs]["le0"][period_bs][5][i])/aod_historical[hist_bs]["le0"][period_bs][5][i]*100.).array[0])
		if model=="": sig_historical=loc_func.ttest4mean(aod_historical[hist_bs]["ge2"][period_bs][0],aod_historical[hist_bs][comp_bs][period_bs][0],aod_historical[hist_bs]["ge2"][period_bs][2],aod_historical[hist_bs][comp_bs][period_bs][2],12,12)
		if model=="CanESM5": sig_historical=loc_func.ftest(aod_historical[hist_bs]["ge2"][period_bs][2],mm_historical[hist_bs]["ge2"][period_bs][2],10,10)

		diff_scenario={};sig_scenario={};
		if model=="":
			for scenario in scenarios:
				source_ids=[]
				for fld in aod_historical[scenario]["ge1"][period_bs][5]:
					source_ids.append(fld.get_property("source_id"))

				diff_scenario[scenario]=[];sig_scenario[scenario]=[];
				for period in ["30","40","50","95"]:
					print(scenario,period)
					tmp_lst=[]
					for fld in aod_scenario[scenario][period]["ge2"][5]:
						try:
							source_id=fld.get_property("source_id")
							idx=source_ids.index(source_id)
							tmp_lst.append(((fld-aod_scenario[scenario][period]["le0"][5][idx])/aod_scenario[scenario][period]["le0"][5][idx]*100.).array[0])
						except ValueError:
							continue
					diff_scenario[scenario].append(np.nanmean(np.array(tmp_lst)))
				#	sig_scenario[scenario].append(loc_func.ttest4mean(aod_scenario[scenario][period]["ge2"][0],aod_scenario[scenario][period][comp_bs][0],aod_scenario[scenario][period]["ge2"][2],aod_scenario[scenario][period][comp_bs][2],12,12))
		if model=="CanESM5":
			for scenario in scenarios:
				source_ids=[]
				for fld in aod_historical[scenario]["ge2"][period_bs][5]:
					source_ids.append(fld.get_property("variant_label"))

				diff_scenario[scenario]=[];sig_scenario[scenario]=[];
				for period in ["30","40","50","95"]:
					print(scenario,period)
					tmp_lst=[]
					for fld in aod_scenario[scenario][period]["ge2"][5]:
						source_id=fld.get_property("variant_label")
						idx=source_ids.index(source_id)
						if fld.mask.array[0]:
							print(source_id,idx,fld.array[0])
						else:
							tmp_lst.append(((fld-aod_scenario[scenario][period]["le0"][5][idx])/aod_scenario[scenario][period]["le0"][5][idx]*100.).array[0])
					diff_scenario[scenario].append(np.nanmean(np.array(tmp_lst)))
					sig_scenario[scenario].append(loc_func.ftest(aod_scenario[scenario][period]["ge2"][2],mm_scenario[scenario][period]["ge2"][2],20,20))
		box={}
		box["hist"]=plt.plot(-1,np.nanmean(diff_historical),marker="s",color="k",markersize=6,linestyle="none",alpha=.8,markerfacecolor="none",markeredgewidth=2)
		if model=="CanESM5":
			plt.plot(-1,np.nanmean(diff_historical),marker="s",color="k",markersize=6,alpha=.8,markerfacecolor="k",markeredgewidth=0)

		pos=[1,3,5,7]
		for i,scenario in enumerate(scenarios):
			box[scenario]=plt.plot(pos,diff_scenario[scenario],marker='s',color=boxcolors[i],markersize=6,linestyle=":",linewidth=2,alpha=.8,markerfacecolor="none",markeredgewidth=2)
			if model=="CanESM5":
				for j,sig in enumerate(sig_scenario[scenario]):
					if sig<=.1:
						print(scenario,j,sig)
						plt.plot(pos[j],diff_scenario[scenario][j],marker="s",color=boxcolors[i],markersize=6,alpha=.8,markerfacecolor=boxcolors[i],markeredgewidth=0)
		scenarios.insert(0,"hist")
		plt.legend((box[scenario][0] for scenario in scenarios),legend_labels,loc='upper right',fontsize=15,ncol=2,framealpha=.4)

		plt.xlim(-1.5,8.);plt.ylim(0,99)
		plt.plot([-1.5,8.],[0,0],c='k',lw=1,ls=':')
		plt.xticks([-1,1,3,5,7],["1979-2014",'2025-2034','2035-2044','2045-2054','2090-2099'])
		plt.ylabel('AOD@550nm Anomaly (%)',fontsize=15)
		plt.gca().tick_params(labelsize=14,direction="out",bottom=True,left=True)
	
	#	plt.savefig('./plot/hwi_scenariomip/fig_aod_scenario_ct_anomaly_{}.pdf'.format(period_bs),bbox_inches='tight',pad_inches=0) if model=="" else plt.savefig('./plot/hwi_scenariomip/fig_aod_scenario_ct_anomaly_{}_{}.pdf'.format(period_bs,model),bbox_inches='tight',pad_inches=0)
		plt.show()
		plt.close()
	'''
	Plot AOD anomaly against comteporary - relative changes
	'''
	if False:
		diff_historical=[];
		for i,fld in enumerate(aod_historical[hist_bs]["ge1"][period_bs][5]):
			diff_historical.append(((fld-aod_historical[hist_bs]["le0"][period_bs][5][i])/aod_historical[hist_bs]["le0"][period_bs][5][i]*100.).array[0])
		if model=="": sig_historical=loc_func.ttest4mean(aod_historical[hist_bs]["ge1"][period_bs][0],aod_historical[hist_bs][comp_bs][period_bs][0],aod_historical[hist_bs]["ge1"][period_bs][2],aod_historical[hist_bs][comp_bs][period_bs][2],12,12)
		if model=="CanESM5": sig_historical=loc_func.ftest(aod_historical[hist_bs]["ge1"][period_bs][2],mm_historical[hist_bs]["ge1"][period_bs][2],10,10)

		diff_scenario={};sig_scenario={};
		if model=="":
			for scenario in scenarios:
				source_ids=[]
				for fld in aod_historical[scenario]["ge1"][period_bs][5]:
					source_ids.append(fld.get_property("source_id"))

				diff_scenario[scenario]=[];sig_scenario[scenario]=[];
				for period in ["30","40","50","95"]:
					print(scenario,period)
					tmp_lst=[]
					for fld in aod_scenario[scenario][period]["ge1"][5]:
						source_id=fld.get_property("source_id")
						idx=source_ids.index(source_id)
						tmp_lst.append(((fld-aod_scenario[scenario][period]["le0"][5][idx])/aod_scenario[scenario][period]["le0"][5][idx]*100.).array[0])
					diff_scenario[scenario].append(np.array(tmp_lst))
				#	sig_scenario[scenario].append(loc_func.ttest4mean(aod_scenario[scenario][period]["ge1"][0],aod_scenario[scenario][period][comp_bs][0],aod_scenario[scenario][period]["ge1"][2],aod_scenario[scenario][period][comp_bs][2],12,12))
		if model=="CanESM5":
			for scenario in scenarios:
				source_ids=[]
				for fld in aod_historical[scenario]["ge1"][period_bs][5]:
					source_ids.append(fld.get_property("variant_label"))

				diff_scenario[scenario]=[];sig_scenario[scenario]=[];
				for period in ["30","40","50","95"]:
					print(scenario,period)
					tmp_lst=[]
					for fld in aod_scenario[scenario][period]["ge1"][5]:
						source_id=fld.get_property("variant_label")
						idx=source_ids.index(source_id)
						if fld.mask.array[0]:
							print(source_id,idx,fld.array[0])
						else:
							tmp_lst.append(((fld-aod_scenario[scenario][period]["le0"][5][idx])/aod_scenario[scenario][period]["le0"][5][idx]*100.).array[0])
					diff_scenario[scenario].append(np.array(tmp_lst))
					sig_scenario[scenario].append(loc_func.ftest(aod_scenario[scenario][period]["ge1"][2],mm_scenario[scenario][period]["ge1"][2],20,20))
		box={}
		medianprops=dict(linewidth=2,color='k')
		boxprops = dict(color="k",facecolor="w",alpha=.4) if sig_historical>.1 else dict(color="k",facecolor="k",alpha=.4)
		capprops=dict(linestyle="-",linewidth=1,color="k")
		plt.boxplot(diff_historical,whis=(5,95),positions=[-1],widths=.1,showcaps=True,showfliers=False,patch_artist=True,medianprops=medianprops,boxprops=boxprops,capprops=capprops)

		for i,scenario in enumerate(scenarios):
			boxprops=dict(color=boxcolors[i],facecolor=boxcolors[i],alpha=.4)
			pos=[.7+i*0.2, 2.6+i*0.2, 4.6+i*0.2, 6.6+i*0.2]
			for j,sig in enumerate(sig_scenario[scenario]):
				boxprops = dict(color=boxcolors[i],facecolor="w",linewidth=2,alpha=.4) if sig>0.1 else dict(color=boxcolors[i],facecolor=boxcolors[i],linewidth=2,alpha=.4)
				plt.boxplot(diff_scenario[scenario][j],whis=(5,95),widths=.1,positions=[pos[j]],showcaps=True,showfliers=False,patch_artist=True,medianprops=medianprops,boxprops=boxprops,capprops=capprops)
		plt.ylim(0,149)
		plt.xlim(-1.5,8.)
		plt.plot([-1.5,8.],[0,0],c='k',lw=1,ls=':')
		plt.xticks([-1,1,3,5,7],["1979-2014",'2025-2034','2035-2044','2045-2054','2090-2099'])
		plt.ylabel('AOD@550nm Anomaly (%)',fontsize=15)
		plt.gca().tick_params(labelsize=14,direction="out",bottom=True,left=True)
	
		plt.savefig('./plot/hwi_scenariomip/fig_aod_scenario_ct_anomaly_{}.pdf'.format(period_bs),bbox_inches='tight',pad_inches=0) if model=="" else plt.savefig('./plot/hwi_scenariomip/fig_aod_scenario_ct_anomaly_{}_{}.pdf'.format(period_bs,model),bbox_inches='tight',pad_inches=0)
		plt.show()
		plt.close()
	'''
	Plot AOD
	'''
	if False:
		aod_scenario_cli_mean={};aod_scenario_cli_median={};aod_scenario_cli_std={};aod_scenario_cli_quartile={};aod_scenario_cli={}
		for scenario in scenarios:
			aod_scenario_cli_mean[scenario]=[];aod_scenario_cli_median[scenario]=[];aod_scenario_cli_std[scenario]=[];aod_scenario_cli_quartile[scenario]=[];aod_scenario_cli[scenario]=[]
			for period in ["30","40","50","95"]:
				aod_scenario_cli_mean[scenario].append(aod_scenario[scenario][period]["cli"][0])
				aod_scenario_cli_median[scenario].append(aod_scenario[scenario][period]["cli"][1])
				aod_scenario_cli_std[scenario].append(aod_scenario[scenario][period]["cli"][2])
				aod_scenario_cli_quartile[scenario].append(aod_scenario[scenario][period]["cli"][3])
				aod_scenario_cli[scenario].append(aod_scenario[scenario][period]["cli"][4])
			'''
			Swap axes for plt.errorbar
			'''
			aod_scenario_cli_quartile[scenario]=np.array(aod_scenario_cli_quartile[scenario]).swapaxes(0,1)
			aod_scenario_cli[scenario]=np.array(aod_scenario_cli[scenario]).swapaxes(0,1)

		eb_kwargs=dict(marker="s",color="k",markersize=6,linestyle="none",capsize=5,elinewidth=1,alpha=.8)
		em_kwargs=dict(marker="o",color="k",markersize=3,alpha=.5)
		box={}
		if mean=="mean":
			box["pd"]=plt.errorbar(-1.,aod_historical[hist_bs]["cli"][period_bs][0],yerr=aod_historical[hist_bs]["cli"][period_bs][2],**eb_kwargs)
		elif mean=="median":
			box["pd"]=plt.errorbar(-1.,aod_historical[hist_bs]["cli"][period_bs][1],yerr=aod_historical[hist_bs]["cli"][period_bs][3].reshape((2,1)),**eb_kwargs)
		for pt in aod_historical[hist_bs]["cli"][period_bs][4]:
			plt.plot(-1,pt,**em_kwargs)

		for i,scenario in enumerate(scenarios):
			pos=np.array([.5,2.5,4.5,6.5])+i*.2
			eb_kwargs["color"]=boxcolors[i];em_kwargs["color"]=boxcolors[i]
			if mean=="mean":
				box[scenario]=plt.errorbar(pos,aod_scenario_cli_mean[scenario],yerr=aod_scenario_cli_std[scenario],**eb_kwargs)
			elif mean=="median":
				box[scenario]=plt.errorbar(pos,aod_scenario_cli_median[scenario],yerr=aod_scenario_cli_quartile[scenario],**eb_kwargs)
			for ip,po in enumerate(pos):
				for pt in aod_scenario_cli[scenario][:,ip]:
					plt.plot(po,pt,**em_kwargs)

		if model=="": plt.ylim(.0,.65)
		if model=="CanESM5": plt.ylim(.22,.87)
		plt.xlim(-1.5,8.)
		if period_bs=="short": plt.xticks([-1,1,3,5,7],['2005-2014','2025-2034','2035-2044','2045-2054','2090-2099'])
		elif period_bs=="long": plt.xticks([-1,1,3,5,7],['1979-2014','2025-2034','2035-2044','2045-2054','2090-2099'])
		plt.ylabel('AOD@550nm',fontsize=15)
		plt.gca().tick_params(labelsize=14,direction="out",bottom=True,left=True)
		plt.legend([box[experiment] for experiment in legend_ct_handles],legend_labels,loc='upper right',fontsize=13,ncol=3,framealpha=.4)

		if model=="": plt.savefig('./plot/hwi_scenariomip/fig_aod_scenario_{}_{}.pdf'.format(mean,period_bs),bbox_inches='tight',pad_inches=0)
		else: plt.savefig('./plot/hwi_scenariomip/fig_aod_scenario_{}_{}_{}.pdf'.format(mean,period_bs,model),bbox_inches='tight',pad_inches=0)
		plt.show()
		plt.close()
	'''
	Plot AOD anomaly against PD
	'''
	if False:
		diff_historical=aod_historical[hist_bs]["ge1"][period_bs][0]-aod_historical[hist_bs][comp_bs][period_bs][0]
		sig_historical=loc_func.ttest4mean(aod_historical[hist_bs]["ge1"][period_bs][0],aod_historical[hist_bs][comp_bs][period_bs][0],aod_historical[hist_bs]["ge1"][period_bs][1],aod_historical[hist_bs][comp_bs][period_bs][1],36,36)

		diff_scenario_pd={};sig_scenario_pd={}
		for scenario in scenarios:
			diff_scenario_pd[scenario]=[];sig_scenario_pd[scenario]=[]
			for period in ["30","40","50","95"]:
				print(scenario,period)
				diff_scenario_pd[scenario].append(aod_scenario[scenario][period]["ge1"][0]-aod_historical[scenario][comp_bs][period_bs][0])
				sig_scenario_pd[scenario].append(loc_func.ttest4mean(aod_scenario[scenario][period]["ge1"][0],aod_historical[scenario][comp_bs][period_bs][0],aod_scenario[scenario][period]["ge1"][1],aod_historical[scenario][comp_bs][period_bs][1],10,36))
		box={}
		if sig_historical<=.12: kwargs=dict(markerfacecolor="k",markeredgewidth=0)
		else: kwargs=dict(markerfacecolor="none",markeredgewidth=2)
		box["pd"]=plt.plot(-1.,diff_historical,marker="s",color="k",markersize=6,alpha=.8,**kwargs)

		for i,scenario in enumerate(scenarios):
			pos=[1.,3.,5.,7.]
			box[scenario+"pd"]=plt.plot(pos,diff_scenario_pd[scenario],marker='s',color=boxcolors[i],markersize=6,linestyle=":",linewidth=2,alpha=.8,markerfacecolor="none",markeredgewidth=2)
			for ip,sig in enumerate(sig_scenario_pd[scenario]):
				print("Significant level:",scenario,ip,sig)
				if sig<=.12:
					print("In ...")
					plt.plot(pos[ip],diff_scenario_pd[scenario][ip],marker="s",color=boxcolors[i],markersize=6,alpha=.8,markerfacecolor=boxcolors[i],markeredgewidth=0)

		if model=="": plt.ylim(-.17,.13)
		if model=="CanESM5": plt.ylim(-.1,.4)
		plt.xlim(-1.5,8.)
		if period_bs=="short": plt.xticks([-1,1,3,5,7],['2005-2014','2025-2034','2035-2044','2045-2054','2090-2099'])
		if period_bs=="long": plt.xticks([-1,1,3,5,7],['1979-2014','2025-2034','2035-2044','2045-2054','2090-2099'])
		plt.ylabel('AOD@550nm Anomaly',fontsize=15)
		plt.gca().tick_params(labelsize=14,direction="out",bottom=True,left=True)
	
		if model=="": plt.savefig('./plot/hwi_scenariomip/fig_aod_scenario_pd_anomaly_{}.pdf'.format(period_bs),bbox_inches='tight',pad_inches=0)
		else: plt.savefig('./plot/hwi_scenariomip/fig_aod_scenario_pd_anomaly_{}_{}.pdf'.format(period_bs,model),bbox_inches='tight',pad_inches=0)
		plt.show()
		plt.close()
	'''
	Plot AOD anomaly against comtemporary baseline
	'''
	if False:
		diff_historical=aod_historical[hist_bs]["ge1"][period_bs][0]-aod_historical[hist_bs][comp_bs][period_bs][0]
		sig_historical=loc_func.ttest4mean(aod_historical[hist_bs]["ge1"][period_bs][0],aod_historical[hist_bs][comp_bs][period_bs][0],aod_historical[hist_bs]["ge1"][period_bs][1],aod_historical[hist_bs][comp_bs][period_bs][1],36,36)

		diff_scenario={};sig_scenario={}
		for scenario in scenarios:
			diff_scenario[scenario]=[];sig_scenario[scenario]=[]
			for period in ["30","40","50","95"]:
				print(scenario,period)
				diff_scenario[scenario].append(aod_scenario[scenario][period]["ge1"][0]-aod_scenario[scenario][period][comp_bs][0])
				sig_scenario[scenario].append(loc_func.ttest4mean(aod_scenario[scenario][period]["ge1"][0],aod_scenario[scenario][period][comp_bs][0],aod_scenario[scenario][period]["ge1"][1],aod_scenario[scenario][period][comp_bs][1],10,10))
		box={}
		if sig_historical<=.12: kwargs=dict(markerfacecolor="k",markeredgewidth=0)
		else: kwargs=dict(markerfacecolor="none",markeredgewidth=2)
		box["pd"]=plt.plot(-1.,diff_historical,marker="s",color="k",markersize=6,alpha=.8,**kwargs)

		for i,scenario in enumerate(scenarios):
			pos=[1.,3.,5.,7.]
			box[scenario]=plt.plot(pos,diff_scenario[scenario],marker='s',color=boxcolors[i],markersize=6,linestyle=":",linewidth=2,alpha=.8,markerfacecolor="none",markeredgewidth=2)
			for ip,sig in enumerate(sig_scenario[scenario]):
				print("Significant level:",scenario,ip,sig)
				if sig<=.12:
					print("In ...")
					plt.plot(pos[ip],diff_scenario[scenario][ip],marker="s",color=boxcolors[i],markersize=6,alpha=.8,markerfacecolor=boxcolors[i],markeredgewidth=0)

		if model=="": plt.ylim(-.17,.13)
		if model=="CanESM5": plt.ylim(-.1,.4)
		plt.xlim(-1.5,8.)
		if period_bs=="short": plt.xticks([-1,1,3,5,7],['2005-2014','2025-2034','2035-2044','2045-2054','2090-2099'])
		elif period_bs=="long": plt.xticks([-1,1,3,5,7],['1979-2014','2025-2034','2035-2044','2045-2054','2090-2099'])
		plt.ylabel('AOD@550nm Anomaly',fontsize=15)
		plt.gca().tick_params(labelsize=14,direction="out",bottom=True,left=True)
	
		if model=="": plt.savefig('./plot/hwi_scenariomip/fig_aod_scenario_ct_anomaly_{}.pdf'.format(period_bs),bbox_inches='tight',pad_inches=0)
		else: plt.savefig('./plot/hwi_scenariomip/fig_aod_scenario_ct_anomaly_{}_{}.pdf'.format(period_bs,model),bbox_inches='tight',pad_inches=0)
		plt.show()
		plt.close()
'''
AOD diagnosis for CanESM5 shows a rapid drop from present-day to near future. In this script, the time series of AOD over Beijing from 1850-2100 for all scenarios are investigate to find out whether this drop is reasonable?
'''
model="CanESM5"
if False:
	'''
	Sort CMIP/Historical
	'''
	source_ids    =[x.split("_")[3] for x in glob.glob("/home/users/lguo/cssp_aero/data/cmip6/cmip/historical/od550aer/beijing/od550aer_AERmon_{}_historical_*_gn_185001-201412.nc".format(model))]
	variant_labels=[x.split("_")[5] for x in glob.glob("/home/users/lguo/cssp_aero/data/cmip6/cmip/historical/od550aer/beijing/od550aer_AERmon_{}_historical_*_gn_185001-201412.nc".format(model))]
	'''
	Historical DJF
	'''
	od550aer_historical={}
	for source_id,variant_label in zip(source_ids,variant_labels):
		try:
			print("Historical",source_id,variant_label)

			od550aer_historical[variant_label]=cf.read("/home/users/lguo/cssp_aero/data/cmip6/cmip/historical/od550aer/beijing/od550aer_AERmon_{}_historical_{}_gn_185001-201412.nc".format(source_id,variant_label)).select_field("atmosphere_optical_thickness_due_to_ambient_aerosol_particles")
		except:
			print("{}_{} does not available in hwi ...".format(source_id,variant_label))
	'''
	Sort ScenarioMIP
	'''
	od550aer={}
	for scenario in scenarios:
		source_ids    =[x.split("_")[3] for x in glob.glob("/home/users/lguo/cssp_aero/data/cmip6/ScenarioMIP/{0}/od550aer/beijing/od550aer_AERmon_{1}_{0}_*_gn_201501-210012.nc".format(scenario,model))]
		variant_labels=[x.split("_")[5] for x in glob.glob("/home/users/lguo/cssp_aero/data/cmip6/ScenarioMIP/{0}/od550aer/beijing/od550aer_AERmon_{1}_{0}_*_gn_201501-210012.nc".format(scenario,model))]

		od550aer_scenario=[];od550aer_historical_scenario=[]
		for source_id,variant_label in zip(source_ids,variant_labels):
			try:
				print(scenario,source_id,variant_label)
				od550aer_scenario.append(cf.read("/home/users/lguo/cssp_aero/data/cmip6/ScenarioMIP/{0}/od550aer/beijing/od550aer_AERmon_{1}_{0}_{2}_gn_201501-210012.nc".format(scenario,source_id,variant_label)).select_field("atmosphere_optical_thickness_due_to_ambient_aerosol_particles"))
				od550aer_historical_scenario.append(od550aer_historical[variant_label])
			except:
				print("{}_{}_{} does not have corresponding historical run ...".format(scenario,source_id,variant_label))
		od550aer_scenario_mean=sum(od550aer_scenario)/len(od550aer_scenario)
		od550aer_historical_mean=sum(od550aer_historical_scenario)/len(od550aer_historical_scenario)
		'''
		Differences in CellMeasure and CellMethod prevent cf.aggregate().
		Delete differences
		'''
		od550aer_historical_mean.del_construct(identity='measure:area')
		od550aer_historical_mean.del_construct(identity="cellmethod1");od550aer_scenario_mean.del_construct(identity="cellmethod1")
		'''
		Aggregate time serieses
		'''
		od550aer[scenario]=cf.aggregate([od550aer_historical_mean,od550aer_scenario_mean])[0].collapse("mean","T",group=cf.djf())
	'''
	Temporary storage
	'''
	with open("/home/users/lguo/cssp_aero/analysis_script/data/od550aer_{}.pickle".format(model),"wb") as f:
		pickle.dump(od550aer,f)
if False:
	with open("/home/users/lguo/cssp_aero/analysis_script/data/od550aer_{}.pickle".format(model),"rb") as f:
		od550aer=pickle.load(f)
	'''
	Plot
	'''
	colors=['#1E9684','#1D3354','#EADD3D','#F21111','#840B22']
	scenario_labels=["SSP1-1.9","SSP1-2.6","SSP2-4.5","SSP3-7.0","SSP5-8.5"]
	cfp.gopen()
	'''
	Reset axis label font size
	'''
	cfp.setvars(axis_label_fontsize=15)
	for i,scenario in enumerate(scenarios):
		if i==(len(scenarios)-1):
			cfp.lineplot(od550aer[scenario].moving_window("mean",10,axis="T"),color=colors[i],linewidth=2,label=scenario_labels[i],ylabel="AOD@550nm",legend_location='upper right')
		else:
			cfp.lineplot(od550aer[scenario].moving_window("mean",10,axis="T"),color=colors[i],linewidth=2,label=scenario_labels[i])
	'''
	Retreat time dimention for referencing line
	'''
	time=od550aer["ssp119"].convert("time")
	time.subspace(T=cf.dt("2014-01-15")).array[0]
	xpts=[time.subspace(T=cf.dt("2014-01-15")).array[0],time.subspace(T=cf.dt("2014-01-15")).array[0]];ypts=[0,1]
#	cfp.lineplot(x=xpts,y=ypts,color="k",linewidth=1,linestyle=":")
	cfp.plotvars.plot.plot(xpts,ypts,lw=1,ls=":",color='k')
#	cfp.plotvars.master_plot.savefig('./plot/hwi_scenariomip/od550aer_{}.pdf'.format(model),bbox_inches='tight',pad_inches=0)
	cfp.gclose(view=False)
'''
Investigate PM2.5:
	With mm calculated from matched ensemble members
'''
if False:
	result_historical={};result_scenario={}
	for scenario in scenarios[1:]:
		'''
		Sort Historical
		'''
		source_ids    =[x.split("/")[12].split("_")[0] for x in glob.glob("{}/cmip/historical/mmr/*_*_monthly_mean_surf_PM2pt5_averaged_over_Beijing_area_for_CMIP6_historical_1850_2014.nc".format(cmipdir))]
		variant_labels=[x.split("/")[12].split("_")[1] for x in glob.glob("{}/cmip/historical/mmr/*_*_monthly_mean_surf_PM2pt5_averaged_over_Beijing_area_for_CMIP6_historical_1850_2014.nc".format(cmipdir))]
		'''
		Historical DJF
		'''
		od550aer_historical={};od550aer_historical["long"]=[];od550aer_historical["short"]=[];
		od550aer_historical_ge1={};od550aer_historical_ge1["long"]=[];od550aer_historical_ge1["short"]=[];
		od550aer_historical_le0={};od550aer_historical_le0["long"]=[];od550aer_historical_le0["short"]=[];
		sources_n_variants=[] # used for matching historical with scenario
		for source_id,variant_label in zip(source_ids,variant_labels):
			print(scenario,"Historical",source_id,variant_label)
			finm="{}/cmip/historical/mmr/{}_{}_monthly_mean_surf_PM2pt5_averaged_over_Beijing_area_for_CMIP6_historical_1850_2014.nc".format(cmipdir,source_id,variant_label)
			print(finm)
			od550aer=cf.read(finm).select_field("long_name=monthly mean surface PM2.5 concentration in historical")
			od550aer_long=od550aer.subspace(T=cf.wi(cf.dt("1979-01-01"),cf.dt("2014-12-30"))).subspace(T=cf.djf())
			od550aer_short=od550aer.subspace(T=cf.wi(cf.dt("2005-01-01"),cf.dt("2014-12-30"))).subspace(T=cf.djf())

			finm="{0}/ScenarioMIP/{1}/hwi/hwi_{2}_{1}_{3}_gn_197901-210012.nc".format(cmipdir,scenario,source_id,variant_label)
			print(finm)
			try:
				hwi=cf.read(finm).select_field("hwi")
			except FileNotFoundError:
				print("{}_{} does not available in hwi ...".format(source_id,variant_label))
			hwi_long=hwi.subspace(T=cf.wi(cf.dt("1979-01-01"),cf.dt("2014-12-30")))
			hwi_short=hwi.subspace(T=cf.wi(cf.dt("2005-01-01"),cf.dt("2014-12-30")))
			'''
			For climatology
			'''
			od550aer_historical["long"].append(od550aer_long.collapse("mean","T"))
			od550aer_historical["short"].append(od550aer_short.collapse("mean","T"))
			'''
			For severe haze months
			'''
			od550aer_long_ge1=od550aer_long.where(hwi_long<1,cf.masked)
			od550aer_historical_ge1["long"].append(od550aer_long_ge1.collapse("mean","T"))
			od550aer_short_ge1=od550aer_short.where(hwi_short<1,cf.masked)
			od550aer_historical_ge1["short"].append(od550aer_short_ge1.collapse("mean","T"))
			'''
			For non-haze months
			'''
			od550aer_long_le0=od550aer_long.where(hwi_long>0,cf.masked)
			od550aer_historical_le0["long"].append(od550aer_long_le0.collapse("mean","T"))
			od550aer_short_le0=od550aer_short.where(hwi_short>0,cf.masked)
			od550aer_historical_le0["short"].append(od550aer_short_le0.collapse("mean","T"))
			'''
			Store sourceid/variant_label
			'''
			sources_n_variants.append("/".join([source_id,variant_label]))
		'''
		Sort ScenarioMIP
		'''
		source_ids    =[x.split("/")[12].split("_")[0] for x in glob.glob("{0}/ScenarioMIP/{1}/mmr/*_*_monthly_mean_surf_PM2pt5_averaged_over_Beijing_area_for_CMIP6_{1}_2015_2100_all.nc".format(cmipdir,scenario))]
		variant_labels=[x.split("/")[12].split("_")[1] for x in glob.glob("{0}/ScenarioMIP/{1}/mmr/*_*_monthly_mean_surf_PM2pt5_averaged_over_Beijing_area_for_CMIP6_{1}_2015_2100_all.nc".format(cmipdir,scenario))]

		od550aer_scenario={};od550aer_scenario["30"]=[];od550aer_scenario["40"]=[];od550aer_scenario["50"]=[];od550aer_scenario["95"]=[]
		od550aer_scenario_ge1={};od550aer_scenario_ge1["30"]=[];od550aer_scenario_ge1["40"]=[];od550aer_scenario_ge1["50"]=[];od550aer_scenario_ge1["95"]=[]
		od550aer_scenario_le0={};od550aer_scenario_le0["30"]=[];od550aer_scenario_le0["40"]=[];od550aer_scenario_le0["50"]=[];od550aer_scenario_le0["95"]=[]

		od550aer_scenario_historical={};od550aer_scenario_historical["long"]=[];od550aer_scenario_historical["short"]=[]
		od550aer_scenario_historical_ge1={};od550aer_scenario_historical_ge1["long"]=[];od550aer_scenario_historical_ge1["short"]=[]
		od550aer_scenario_historical_le0={};od550aer_scenario_historical_le0["long"]=[];od550aer_scenario_historical_le0["short"]=[]
		for source_id,variant_label in zip(source_ids,variant_labels):
			print(scenario,scenario,source_id,variant_label)
			'''
			Match to historical
			'''
			try:
				index_historical=sources_n_variants.index("/".join([source_id,variant_label]))
				print(index_historical)
			except:
				print("... Does not exist in Historical ...")
				continue
			od550aer_scenario_historical["long"].append(od550aer_historical["long"][index_historical])
			od550aer_scenario_historical["short"].append(od550aer_historical["short"][index_historical])
			od550aer_scenario_historical_ge1["long"].append(od550aer_historical_ge1["long"][index_historical])
			od550aer_scenario_historical_ge1["short"].append(od550aer_historical_ge1["short"][index_historical])
			od550aer_scenario_historical_le0["long"].append(od550aer_historical_le0["long"][index_historical])
			od550aer_scenario_historical_le0["short"].append(od550aer_historical_le0["short"][index_historical])

			print("... Continue processing Scenario ...")
			finm="{0}/ScenarioMIP/{1}/mmr/{2}_{3}_monthly_mean_surf_PM2pt5_averaged_over_Beijing_area_for_CMIP6_{1}_2015_2100_all.nc".format(cmipdir,scenario,source_id,variant_label)
			print(finm)
			od550aer=cf.read(finm).select_field("long_name=monthly mean surface PM2.5 concentration in {}".format(scenario)).subspace(time=cf.djf())
			finm="{0}/ScenarioMIP/{1}/hwi/hwi_{2}_{1}_{3}_gn_197901-210012.nc".format(cmipdir,scenario,source_id,variant_label)
			print(finm)
			hwi=cf.read(finm).select_field("hwi").subspace(T=cf.wi(cf.dt("2015-01-01"),cf.dt("2100-12-30")))
			'''
			For climatology
			'''
			print("-----> {} {} <-----".format(hwi.get_property("source_id"),hwi.get_property("variant_label")))
			od550aer_scenario["30"].append(od550aer.subspace(time=cf.wi(cf.dt("2025-01-01"),cf.dt("2034-12-30"))).collapse("mean","time"))
			od550aer_scenario["40"].append(od550aer.subspace(time=cf.wi(cf.dt("2035-01-01"),cf.dt("2044-12-30"))).collapse("mean","time"))
			od550aer_scenario["50"].append(od550aer.subspace(time=cf.wi(cf.dt("2045-01-01"),cf.dt("2054-12-30"))).collapse("mean","time"))
			od550aer_scenario["95"].append(od550aer.subspace(time=cf.wi(cf.dt("2090-01-01"),cf.dt("2099-12-30"))).collapse("mean","time"))
			'''
			For severe haze months
			'''
			od550aer_ge1=od550aer.where(hwi<1,cf.masked)
			od550aer_scenario_ge1["30"].append(od550aer_ge1.subspace(time=cf.wi(cf.dt("2025-01-01"),cf.dt("2034-12-30"))).collapse("mean","time"))
			od550aer_scenario_ge1["40"].append(od550aer_ge1.subspace(time=cf.wi(cf.dt("2035-01-01"),cf.dt("2044-12-30"))).collapse("mean","time"))
			od550aer_scenario_ge1["50"].append(od550aer_ge1.subspace(time=cf.wi(cf.dt("2045-01-01"),cf.dt("2054-12-30"))).collapse("mean","time"))
			od550aer_scenario_ge1["95"].append(od550aer_ge1.subspace(time=cf.wi(cf.dt("2090-01-01"),cf.dt("2099-12-30"))).collapse("mean","time"))
			'''
			For non-haze months
			'''
			od550aer_le0=od550aer.where(hwi>0,cf.masked)
			od550aer_scenario_le0["30"].append(od550aer_le0.subspace(time=cf.wi(cf.dt("2025-01-01"),cf.dt("2034-12-30"))).collapse("mean","time"))
			od550aer_scenario_le0["40"].append(od550aer_le0.subspace(time=cf.wi(cf.dt("2035-01-01"),cf.dt("2044-12-30"))).collapse("mean","time"))
			od550aer_scenario_le0["50"].append(od550aer_le0.subspace(time=cf.wi(cf.dt("2045-01-01"),cf.dt("2054-12-30"))).collapse("mean","time"))
			od550aer_scenario_le0["95"].append(od550aer_le0.subspace(time=cf.wi(cf.dt("2090-01-01"),cf.dt("2099-12-30"))).collapse("mean","time"))
		'''
		Calculate model mean
		'''
		od550aer_scenario_historical["long"]=loc_func.fldlst2lstmm_new(od550aer_scenario_historical["long"])
		od550aer_scenario_historical["short"]=loc_func.fldlst2lstmm_new(od550aer_scenario_historical["short"])
		od550aer_scenario_historical_ge1["long"]=loc_func.fldlst2lstmm_new(od550aer_scenario_historical_ge1["long"])
		od550aer_scenario_historical_ge1["short"]=loc_func.fldlst2lstmm_new(od550aer_scenario_historical_ge1["short"])
		od550aer_scenario_historical_le0["long"]=loc_func.fldlst2lstmm_new(od550aer_scenario_historical_le0["long"])
		od550aer_scenario_historical_le0["short"]=loc_func.fldlst2lstmm_new(od550aer_scenario_historical_le0["short"])

		od550aer_scenario["30"]=loc_func.fldlst2lstmm_new(od550aer_scenario["30"])
		od550aer_scenario["40"]=loc_func.fldlst2lstmm_new(od550aer_scenario["40"])
		od550aer_scenario["50"]=loc_func.fldlst2lstmm_new(od550aer_scenario["50"])
		od550aer_scenario["95"]=loc_func.fldlst2lstmm_new(od550aer_scenario["95"])

		od550aer_scenario_ge1["30"]=loc_func.fldlst2lstmm_new(od550aer_scenario_ge1["30"])
		od550aer_scenario_ge1["40"]=loc_func.fldlst2lstmm_new(od550aer_scenario_ge1["40"])
		od550aer_scenario_ge1["50"]=loc_func.fldlst2lstmm_new(od550aer_scenario_ge1["50"])
		od550aer_scenario_ge1["95"]=loc_func.fldlst2lstmm_new(od550aer_scenario_ge1["95"])

		od550aer_scenario_le0["30"]=loc_func.fldlst2lstmm_new(od550aer_scenario_le0["30"])
		od550aer_scenario_le0["40"]=loc_func.fldlst2lstmm_new(od550aer_scenario_le0["40"])
		od550aer_scenario_le0["50"]=loc_func.fldlst2lstmm_new(od550aer_scenario_le0["50"])
		od550aer_scenario_le0["95"]=loc_func.fldlst2lstmm_new(od550aer_scenario_le0["95"])

		result_historical[scenario]={}
		result_historical[scenario]["cli"]={}
		result_historical[scenario]["cli"]["long"]=[loc_func.mean_fldlst(od550aer_scenario_historical["long"])]
		result_historical[scenario]["cli"]["long"].append(loc_func.median_fldlst(od550aer_scenario_historical["long"]))
		result_historical[scenario]["cli"]["long"].append(loc_func.stdev_fldlst(od550aer_scenario_historical["long"]))
		result_historical[scenario]["cli"]["long"].append(loc_func.quartile_fldlst(od550aer_scenario_historical["long"]))
		result_historical[scenario]["cli"]["long"].append(od550aer_scenario_historical["long"])
		result_historical[scenario]["cli"]["long"].append(loc_func.tolist_fldlst(od550aer_scenario_historical["long"]))
		result_historical[scenario]["cli"]["short"]=[loc_func.mean_fldlst(od550aer_scenario_historical["short"])]
		result_historical[scenario]["cli"]["short"]=[loc_func.median_fldlst(od550aer_scenario_historical["short"])]
		result_historical[scenario]["cli"]["short"].append(loc_func.stdev_fldlst(od550aer_scenario_historical["short"]))
		result_historical[scenario]["cli"]["short"].append(loc_func.quartile_fldlst(od550aer_scenario_historical["short"]))
		result_historical[scenario]["cli"]["short"].append(loc_func.tolist_fldlst(od550aer_scenario_historical["short"]))
		result_historical[scenario]["cli"]["short"].append(od550aer_scenario_historical["short"])
		result_historical[scenario]["ge1"]={}
		result_historical[scenario]["ge1"]["long"]=[loc_func.mean_fldlst(od550aer_scenario_historical_ge1["long"])]
		result_historical[scenario]["ge1"]["long"].append(loc_func.median_fldlst(od550aer_scenario_historical_ge1["long"]))
		result_historical[scenario]["ge1"]["long"].append(loc_func.stdev_fldlst(od550aer_scenario_historical_ge1["long"]))
		result_historical[scenario]["ge1"]["long"].append(loc_func.quartile_fldlst(od550aer_scenario_historical_ge1["long"]))
		result_historical[scenario]["ge1"]["long"].append(loc_func.tolist_fldlst(od550aer_scenario_historical_ge1["long"]))
		result_historical[scenario]["ge1"]["long"].append(od550aer_scenario_historical_ge1["long"])
		result_historical[scenario]["ge1"]["short"]=[loc_func.mean_fldlst(od550aer_scenario_historical_ge1["short"])]
		result_historical[scenario]["ge1"]["short"].append(loc_func.median_fldlst(od550aer_scenario_historical_ge1["short"]))
		result_historical[scenario]["ge1"]["short"].append(loc_func.stdev_fldlst(od550aer_scenario_historical_ge1["short"]))
		result_historical[scenario]["ge1"]["short"].append(loc_func.quartile_fldlst(od550aer_scenario_historical_ge1["short"]))
		result_historical[scenario]["ge1"]["short"].append(loc_func.tolist_fldlst(od550aer_scenario_historical_ge1["short"]))
		result_historical[scenario]["ge1"]["short"].append(od550aer_scenario_historical_ge1["short"])
		result_historical[scenario]["le0"]={}
		result_historical[scenario]["le0"]["long"] =[loc_func.mean_fldlst(od550aer_scenario_historical_le0["long"])]
		result_historical[scenario]["le0"]["long"].append(loc_func.median_fldlst(od550aer_scenario_historical_le0["long"]))
		result_historical[scenario]["le0"]["long"].append(loc_func.stdev_fldlst(od550aer_scenario_historical_le0["long"]))
		result_historical[scenario]["le0"]["long"].append(loc_func.quartile_fldlst(od550aer_scenario_historical_le0["long"]))
		result_historical[scenario]["le0"]["long"].append(loc_func.tolist_fldlst(od550aer_scenario_historical_le0["long"]))
		result_historical[scenario]["le0"]["long"].append(od550aer_scenario_historical_le0["long"])
		result_historical[scenario]["le0"]["short"]=[loc_func.mean_fldlst(od550aer_scenario_historical_le0["short"])]
		result_historical[scenario]["le0"]["short"].append(loc_func.median_fldlst(od550aer_scenario_historical_le0["short"]))
		result_historical[scenario]["le0"]["short"].append(loc_func.stdev_fldlst(od550aer_scenario_historical_le0["short"]))
		result_historical[scenario]["le0"]["short"].append(loc_func.quartile_fldlst(od550aer_scenario_historical_le0["short"]))
		result_historical[scenario]["le0"]["short"].append(loc_func.tolist_fldlst(od550aer_scenario_historical_le0["short"]))
		result_historical[scenario]["le0"]["short"].append(od550aer_scenario_historical_le0["short"])

		result_scenario[scenario]={}
		result_scenario[scenario]["30"]={};
		result_scenario[scenario]["30"]["cli"]=[loc_func.mean_fldlst(od550aer_scenario["30"])]
		result_scenario[scenario]["30"]["cli"].append(loc_func.median_fldlst(od550aer_scenario["30"]))
		result_scenario[scenario]["30"]["cli"].append(loc_func.stdev_fldlst(od550aer_scenario["30"]))
		result_scenario[scenario]["30"]["cli"].append(loc_func.quartile_fldlst(od550aer_scenario["30"]))
		result_scenario[scenario]["30"]["cli"].append(loc_func.tolist_fldlst(od550aer_scenario["30"]))
		result_scenario[scenario]["30"]["cli"].append(od550aer_scenario["30"])
		result_scenario[scenario]["30"]["ge1"]=[loc_func.mean_fldlst(od550aer_scenario_ge1["30"])]
		result_scenario[scenario]["30"]["ge1"].append(loc_func.median_fldlst(od550aer_scenario_ge1["30"]))
		result_scenario[scenario]["30"]["ge1"].append(loc_func.stdev_fldlst(od550aer_scenario_ge1["30"]))
		result_scenario[scenario]["30"]["ge1"].append(loc_func.quartile_fldlst(od550aer_scenario_ge1["30"]))
		result_scenario[scenario]["30"]["ge1"].append(loc_func.tolist_fldlst(od550aer_scenario_ge1["30"]))
		result_scenario[scenario]["30"]["ge1"].append(od550aer_scenario_ge1["30"])
		result_scenario[scenario]["30"]["le0"]=[loc_func.mean_fldlst(od550aer_scenario_le0["30"])]
		result_scenario[scenario]["30"]["le0"].append(loc_func.median_fldlst(od550aer_scenario_le0["30"]))
		result_scenario[scenario]["30"]["le0"].append(loc_func.stdev_fldlst(od550aer_scenario_le0["30"]))
		result_scenario[scenario]["30"]["le0"].append(loc_func.quartile_fldlst(od550aer_scenario_le0["30"]))
		result_scenario[scenario]["30"]["le0"].append(loc_func.tolist_fldlst(od550aer_scenario_le0["30"]))
		result_scenario[scenario]["30"]["le0"].append(od550aer_scenario_le0["30"])
		result_scenario[scenario]["40"]={};
		result_scenario[scenario]["40"]["cli"]=[loc_func.mean_fldlst(od550aer_scenario["40"])]
		result_scenario[scenario]["40"]["cli"].append(loc_func.median_fldlst(od550aer_scenario["40"]))
		result_scenario[scenario]["40"]["cli"].append(loc_func.stdev_fldlst(od550aer_scenario["40"]))
		result_scenario[scenario]["40"]["cli"].append(loc_func.quartile_fldlst(od550aer_scenario["40"]))
		result_scenario[scenario]["40"]["cli"].append(loc_func.tolist_fldlst(od550aer_scenario["40"]))
		result_scenario[scenario]["40"]["cli"].append(od550aer_scenario["40"])
		result_scenario[scenario]["40"]["ge1"]=[loc_func.mean_fldlst(od550aer_scenario_ge1["40"])]
		result_scenario[scenario]["40"]["ge1"].append(loc_func.median_fldlst(od550aer_scenario_ge1["40"]))
		result_scenario[scenario]["40"]["ge1"].append(loc_func.stdev_fldlst(od550aer_scenario_ge1["40"]))
		result_scenario[scenario]["40"]["ge1"].append(loc_func.quartile_fldlst(od550aer_scenario_ge1["40"]))
		result_scenario[scenario]["40"]["ge1"].append(loc_func.tolist_fldlst(od550aer_scenario_ge1["40"]))
		result_scenario[scenario]["40"]["ge1"].append(od550aer_scenario_ge1["40"])
		result_scenario[scenario]["40"]["le0"]=[loc_func.mean_fldlst(od550aer_scenario_le0["40"])]
		result_scenario[scenario]["40"]["le0"].append(loc_func.median_fldlst(od550aer_scenario_le0["40"]))
		result_scenario[scenario]["40"]["le0"].append(loc_func.stdev_fldlst(od550aer_scenario_le0["40"]))
		result_scenario[scenario]["40"]["le0"].append(loc_func.quartile_fldlst(od550aer_scenario_le0["40"]))
		result_scenario[scenario]["40"]["le0"].append(loc_func.tolist_fldlst(od550aer_scenario_le0["40"]))
		result_scenario[scenario]["40"]["le0"].append(od550aer_scenario_le0["40"])
		result_scenario[scenario]["50"]={};
		result_scenario[scenario]["50"]["cli"]=[loc_func.mean_fldlst(od550aer_scenario["50"])]
		result_scenario[scenario]["50"]["cli"].append(loc_func.median_fldlst(od550aer_scenario["50"]))
		result_scenario[scenario]["50"]["cli"].append(loc_func.stdev_fldlst(od550aer_scenario["50"]))
		result_scenario[scenario]["50"]["cli"].append(loc_func.quartile_fldlst(od550aer_scenario["50"]))
		result_scenario[scenario]["50"]["cli"].append(loc_func.tolist_fldlst(od550aer_scenario["50"]))
		result_scenario[scenario]["50"]["cli"].append(od550aer_scenario["50"])
		result_scenario[scenario]["50"]["ge1"]=[loc_func.mean_fldlst(od550aer_scenario_ge1["50"])]
		result_scenario[scenario]["50"]["ge1"].append(loc_func.median_fldlst(od550aer_scenario_ge1["50"]))
		result_scenario[scenario]["50"]["ge1"].append(loc_func.stdev_fldlst(od550aer_scenario_ge1["50"]))
		result_scenario[scenario]["50"]["ge1"].append(loc_func.quartile_fldlst(od550aer_scenario_ge1["50"]))
		result_scenario[scenario]["50"]["ge1"].append(loc_func.tolist_fldlst(od550aer_scenario_ge1["50"]))
		result_scenario[scenario]["50"]["ge1"].append(od550aer_scenario_ge1["50"])
		result_scenario[scenario]["50"]["le0"]=[loc_func.mean_fldlst(od550aer_scenario_le0["50"])]
		result_scenario[scenario]["50"]["le0"].append(loc_func.median_fldlst(od550aer_scenario_le0["50"]))
		result_scenario[scenario]["50"]["le0"].append(loc_func.stdev_fldlst(od550aer_scenario_le0["50"]))
		result_scenario[scenario]["50"]["le0"].append(loc_func.quartile_fldlst(od550aer_scenario_le0["50"]))
		result_scenario[scenario]["50"]["le0"].append(loc_func.tolist_fldlst(od550aer_scenario_le0["50"]))
		result_scenario[scenario]["50"]["le0"].append(od550aer_scenario_le0["50"])
		result_scenario[scenario]["95"]={};
		result_scenario[scenario]["95"]["cli"]=[loc_func.mean_fldlst(od550aer_scenario["95"])]
		result_scenario[scenario]["95"]["cli"].append(loc_func.median_fldlst(od550aer_scenario["95"]))
		result_scenario[scenario]["95"]["cli"].append(loc_func.stdev_fldlst(od550aer_scenario["95"]))
		result_scenario[scenario]["95"]["cli"].append(loc_func.quartile_fldlst(od550aer_scenario["95"]))
		result_scenario[scenario]["95"]["cli"].append(loc_func.tolist_fldlst(od550aer_scenario["95"]))
		result_scenario[scenario]["95"]["cli"].append(od550aer_scenario["95"])
		result_scenario[scenario]["95"]["ge1"]=[loc_func.mean_fldlst(od550aer_scenario_ge1["95"])]
		result_scenario[scenario]["95"]["ge1"].append(loc_func.median_fldlst(od550aer_scenario_ge1["95"]))
		result_scenario[scenario]["95"]["ge1"].append(loc_func.stdev_fldlst(od550aer_scenario_ge1["95"]))
		result_scenario[scenario]["95"]["ge1"].append(loc_func.quartile_fldlst(od550aer_scenario_ge1["95"]))
		result_scenario[scenario]["95"]["ge1"].append(loc_func.tolist_fldlst(od550aer_scenario_ge1["95"]))
		result_scenario[scenario]["95"]["ge1"].append(od550aer_scenario_ge1["95"])
		result_scenario[scenario]["95"]["le0"]=[loc_func.mean_fldlst(od550aer_scenario_le0["95"])]
		result_scenario[scenario]["95"]["le0"].append(loc_func.median_fldlst(od550aer_scenario_le0["95"]))
		result_scenario[scenario]["95"]["le0"].append(loc_func.stdev_fldlst(od550aer_scenario_le0["95"]))
		result_scenario[scenario]["95"]["le0"].append(loc_func.quartile_fldlst(od550aer_scenario_le0["95"]))
		result_scenario[scenario]["95"]["le0"].append(loc_func.tolist_fldlst(od550aer_scenario_le0["95"]))
		result_scenario[scenario]["95"]["le0"].append(od550aer_scenario_le0["95"])
	'''
	Write out
	'''
	with open("./data/pm25_ano.pickle","wb") as fp:
		pickle.dump(result_historical,fp)
		pickle.dump(result_scenario,fp)
'''
Summarise PM2.5 anomalies
'''
if False:
	hist_bs="ssp370"
	comp_bs="le0"
	period_bs="long"#"long"

	scenarios=scenarios[1:]
	boxcolors=['#1E9684','#1D3354','#EADD3D','#F21111','#840B22']
	boxcolors=boxcolors[1:]
	legend_ct_handles=["pd","ssp119","ssp126","ssp245","ssp370","ssp585"]
	legend_ct_handles.remove("ssp119")
	legend_pd_handles=["pd","ssp119pd","ssp126pd","ssp245pd","ssp370pd","ssp585pd"]
	legend_pd_handles.remove("ssp119pd")
	legend_labels=["Hist.","SSP1-1.9","SSP1-2.6","SSP2-4.5","SSP3-7.0","SSP5-8.5"]
	legend_labels.remove("SSP1-1.9")
	'''
	Read
	'''
	fp=open("./data/pm25_ano.pickle","rb")	
	aod_historical=pickle.load(fp)
	aod_scenario=pickle.load(fp)
#	'''
#	Anomaly
#	'''
#	diff_historical=aod_historical[hist_bs]["ge1"][period_bs]-aod_historical[hist_bs][comp_bs][period_bs]
#	diff_scenario={};diff_scenario_pd={};
#	for scenario in scenarios:
#		diff_scenario[scenario]=[];diff_scenario_pd[scenario]=[];
#		for period in ["30","40","50","95"]:
#			print(scenario,period)
#			print(aod_scenario[scenario][period]["ge1"],aod_scenario[scenario][period][comp_bs],aod_scenario[scenario][period]["ge1"]-aod_scenario[scenario][period][comp_bs])
#			diff_scenario[scenario].append(aod_scenario[scenario][period]["ge1"]-aod_scenario[scenario][period][comp_bs])
#			print(aod_scenario[scenario][period]["ge1"],aod_historical[scenario][comp_bs][period_bs],aod_scenario[scenario][period]["ge1"]-aod_historical[scenario][comp_bs][period_bs])
#			diff_scenario_pd[scenario].append(aod_scenario[scenario][period]["ge1"]-aod_historical[scenario][comp_bs][period_bs])
	'''
	Plot
	'''
	if False:
		aod_scenario_cli={}
		for scenario in scenarios:
			aod_scenario_cli[scenario]=[]
			for period in ["30","40","50","95"]:
				aod_scenario_cli[scenario].append(aod_scenario[scenario][period]["cli"])
		box={}
		box["pd"]=plt.scatter(-1.,aod_historical[hist_bs]["cli"][period_bs],c='k',s=40,marker='s',alpha=.8)
		for i,scenario in enumerate(scenarios):
			pos=[1.,3.,5.,7.]
			box[scenario]=plt.scatter(pos,aod_scenario_cli[scenario],c=boxcolors[i],s=40,marker='s',alpha=.8)
			plt.plot(pos,aod_scenario_cli[scenario],c=boxcolors[i],ls=':',lw=3,alpha=.8)
		plt.ylim(5,55)
		plt.xlim(-1.5,8.)
		if period_bs=="short": plt.xticks([-1,1,3,5,7],['2005-2014','2025-2034','2035-2044','2045-2054','2090-2099'])
		elif period_bs=="long": plt.xticks([-1,1,3,5,7],['1979-2014','2025-2034','2035-2044','2045-2054','2090-2099'])
		plt.plot([-1.5,8.],[10,10],c='k',lw=1,ls=':')
		plt.ylabel('PM2.5 ($mg/m^3$)',fontsize=15)
		plt.gca().tick_params(labelsize=14,direction="out",bottom=True,left=True)
		plt.legend([box[experiment] for experiment in legend_ct_handles],legend_labels,loc='upper left',ncol=3,fontsize=13,framealpha=0.4)

		plt.savefig('./plot/hwi_scenariomip/fig_pm25_scenario_{}.pdf'.format(period_bs),bbox_inches='tight',pad_inches=0)
		plt.show()
		plt.close()
	'''
	Plot anomaly against PD
	'''
	if False:
		box={}
		'''
		mark historical AODa
		'''
		box["pd"]=plt.scatter(-1.,diff_historical,c='k',s=40,marker='s',alpha=.8)
		'''
		add markers
		'''
		for i,scenario in enumerate(scenarios):
			pos=[1.,3.,5.,7.]
			box[scenario+"pd"]=plt.scatter(pos,diff_scenario_pd[scenario],c=boxcolors[i],s=40,marker='s',alpha=.8)
			plt.plot(pos,diff_scenario_pd[scenario],c=boxcolors[i],ls=':',lw=3,alpha=.8)
		if period_bs=="long":
			plt.ylim(-30,30)
		elif period_bs=="short":
			plt.ylim(-.18,.12)
		plt.xlim(-1.5,8.)
		if period_bs=="short":
			plt.xticks([-1,1,3,5,7],['2005-2014','2025-2034','2035-2044','2045-2054','2090-2099'])
		elif period_bs=="long":
			plt.xticks([-1,1,3,5,7],['1979-2014','2025-2034','2035-2044','2045-2054','2090-2099'])
		plt.ylabel('PM2.5 Anomaly ($mg/m^3$)',fontsize=15)
		plt.gca().tick_params(labelsize=14,direction="out",bottom=True,left=True)
	
		plt.savefig('./plot/hwi_scenariomip/fig_pm25_scenario_pd_anomaly_{}.pdf'.format(period_bs),bbox_inches='tight',pad_inches=0)
		plt.show()
		plt.close()
	'''
	Plot anomaly against comtemporary
	'''
	if False:
		box={}
		'''
		mark historical anomaly
		'''
		box["pd"]=plt.scatter(-1.,diff_historical,c='k',s=40,marker='s',alpha=.8)
		'''
		add markers
		'''
		for i,scenario in enumerate(scenarios):
			pos=[1.,3.,5.,7.]
			box[scenario]=plt.scatter(pos,diff_scenario[scenario],c=boxcolors[i],s=40,marker='s',alpha=.8)
			plt.plot(pos,diff_scenario[scenario],c=boxcolors[i],ls=':',lw=3,alpha=.8)
		plt.ylim(-30,30)
		plt.xlim(-1.5,8.)
		if period_bs=="short":
			plt.xticks([-1,1,3,5,7],['2005-2014','2025-2034','2035-2044','2045-2054','2090-2099'])
		elif period_bs=="long":
			plt.xticks([-1,1,3,5,7],['1979-2014','2025-2034','2035-2044','2045-2054','2090-2099'])
		plt.ylabel('PM2.5 Anomaly ($mg/m^3$)',fontsize=15)
		plt.gca().tick_params(labelsize=14,direction="out",bottom=True,left=True)
	
		plt.savefig('./plot/hwi_scenariomip/fig_pm25_scenario_ct_anomaly_{}.pdf'.format(period_bs),bbox_inches='tight',pad_inches=0)
		plt.show()
		plt.close()
	'''
	Plot anomaly against PD relative
	'''
	if True:
		diff_scenario={};sig_scenario={}
		for scenario in scenarios:
			source_ids=[]
			for fld in aod_historical[scenario]["ge1"][period_bs][5]:
				source_ids.append(fld.get_property("source_id"))

			diff_scenario[scenario]=[];sig_scenario[scenario]=[]
			for period in ["30","40","50","95"]:
				print(scenario,period)
				tmp_lst=[]
				for fld in aod_scenario[scenario][period]["ge1"][5]:
					source_id=fld.get_property("source_id")
					idx=source_ids.index(source_id)
					tmp_lst.append(((fld-aod_historical[scenario]["ge1"][period_bs][5][idx])/aod_historical[scenario]["ge1"][period_bs][5][idx]*100.).array[0])
				diff_scenario[scenario].append(np.nanmean(np.array(tmp_lst)))
				sig_scenario[scenario].append(loc_func.ttest4mean(aod_scenario[scenario][period]["ge1"][0],aod_historical[scenario][comp_bs][period_bs][0],aod_scenario[scenario][period]["ge1"][2],aod_historical[scenario][comp_bs][period_bs][2],6,6))
		box={}
		pos=[1,3,5,7]
		for i,scenario in enumerate(scenarios):
			print(i,scenario)
			box[scenario]=plt.plot(pos,diff_scenario[scenario],marker='s',color=boxcolors[i],markersize=6,linestyle=":",linewidth=2,alpha=.8,markerfacecolor="none",markeredgewidth=2)
		plt.legend((box[scenario][0] for scenario in scenarios),legend_labels[1:],loc='upper right',fontsize=15,ncol=2,framealpha=.4)

		plt.ylim(-99,99)
		plt.xlim(0,8.)
		plt.plot([0,8.],[0,0],c='k',lw=1,ls=':')
		plt.xticks([1,3,5,7],['2025-2034','2035-2044','2045-2054','2090-2099'])
		plt.ylabel(r"PM$_{2.5}$ Anomaly (%)",fontsize=15)
		plt.gca().tick_params(labelsize=14,direction="out",bottom=True,left=True)
	
		plt.savefig('./plot/hwi_scenariomip/fig_pm25_scenario_pd_anomaly_{}.pdf'.format(period_bs),bbox_inches='tight',pad_inches=0)
		plt.show()
		plt.close()
	'''
	Plot AOD anomaly against comteporary - relative changes
	'''
	if True:
		diff_historical=[];
		for i,fld in enumerate(aod_historical[hist_bs]["ge1"][period_bs][5]):
			diff_historical.append(((fld-aod_historical[hist_bs]["le0"][period_bs][5][i])/aod_historical[hist_bs]["le0"][period_bs][5][i]*100.).array[0])
		sig_historical=loc_func.ttest4mean(aod_historical[hist_bs]["ge1"][period_bs][0],aod_historical[hist_bs][comp_bs][period_bs][0],aod_historical[hist_bs]["ge1"][period_bs][2],aod_historical[hist_bs][comp_bs][period_bs][2],6,6)

		diff_scenario={};sig_scenario={};
		for scenario in scenarios:
			diff_scenario[scenario]=[];sig_scenario[scenario]=[];
			source_ids=[]
			for fld in aod_historical[scenario]["ge1"][period_bs][5]:
				source_ids.append(fld.get_property("source_id"))
			for period in ["30","40","50","95"]:
				print(scenario,period)
				tmp_lst=[]
				for fld in aod_scenario[scenario][period]["ge1"][5]:
					source_id=fld.get_property("source_id")
					idx=source_ids.index(source_id)
					tmp_lst.append(((fld-aod_scenario[scenario][period]["le0"][5][idx])/aod_scenario[scenario][period]["le0"][5][idx]*100.).array[0])
				diff_scenario[scenario].append(np.nanmean(np.array(tmp_lst)))
				sig_scenario[scenario].append(loc_func.ttest4mean(aod_scenario[scenario][period]["ge1"][0],aod_scenario[scenario][period][comp_bs][0],aod_scenario[scenario][period]["ge1"][2],aod_scenario[scenario][period][comp_bs][2],6,6))
		box={}
		box["hist"]=plt.plot(-1,np.nanmean(diff_historical),marker="s",color="k",markersize=6,linestyle="none",alpha=.8,markerfacecolor="none",markeredgewidth=2)
		pos=[1,3,5,7]
		for i,scenario in enumerate(scenarios):
			box[scenario]=plt.plot(pos,diff_scenario[scenario],marker='s',color=boxcolors[i],markersize=6,linestyle=":",linewidth=2,alpha=.8,markerfacecolor="none",markeredgewidth=2)
		scenarios.insert(0,"hist")
		plt.legend((box[scenario][0] for scenario in scenarios),legend_labels,loc='upper right',fontsize=15,ncol=2,framealpha=.4)

		plt.ylim(0,99)
		plt.xlim(-1.5,8.)
		plt.plot([-1.5,8.],[0,0],c='k',lw=1,ls=':')
		plt.xticks([-1,1,3,5,7],["1979-2014",'2025-2034','2035-2044','2045-2054','2090-2099'])
		plt.ylabel(r"PM$_{2.5}$ Anomaly (%)",fontsize=15)
		plt.gca().tick_params(labelsize=14,direction="out",bottom=True,left=True)
	
		plt.savefig('./plot/hwi_scenariomip/fig_pm25_scenario_ct_anomaly_{}.pdf'.format(period_bs),bbox_inches='tight',pad_inches=0)
		plt.show()
		plt.close()
