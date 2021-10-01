
import cf
import numpy as np
from scipy import stats
import statistics
import pdb
'''
Calculate HWI following ###
'''
def cal_hwi(t850,t850_mean,t250,t250_mean,u500_1,u500_1_mean,u500_2,u500_2_mean,v850,v850_mean):
	'''
	calculate anomalies
	'''
	t850=t850-t850_mean
	t250=t250-t250_mean
	u500_1=u500_1-u500_1_mean
	u500_2=u500_2-u500_2_mean
	v850=v850-v850_mean
	'''
	calculate shears
	'''
	delta_t=t850-t250
	delta_u=u500_2-u500_1
	'''
	standardise anomalies
	'''
#	delta_t=(delta_t-delta_t.collapse('mean','T'))/delta_t.collapse('sd','T')
#	delta_u=(delta_u-delta_u.collapse('mean','T'))/delta_u.collapse('sd','T')
#	v850=(v850-v850.collapse('mean','T'))/v850.collapse('sd','T')
	'''
	assemble
	'''
	hwi=delta_t+delta_u+v850
	'''
	standardise hwi
	'''
#	hwi=(hwi-hwi.collapse('T: mean'))/hwi.collapse('T: sd')

	hwi.set_property('standard_name','hwi')

	return hwi,delta_t,delta_u,v850

'''
Check aggregation in cal_hwi_scenario
'''
def check_aggregation(fldlst):
	if len(fldlst)!=1:
		print("cf.aggregate() fails")
		exit()
	else:
		fld=fldlst[0]
	return fld
'''
Calculate HWI:
	1. Concatenate time serieses from 197901-210012;
	2. Calculate HWI
'''
def cal_hwi_scenario(t850_scenario,t850_historical,t250_scenario,t250_historical,u500_1_scenario,u500_1_historical,u500_2_scenario,u500_2_historical,v850_scenario,v850_historical):
	for fld in [t850_scenario,t850_historical,t250_scenario,t250_historical,u500_1_scenario,u500_1_historical,u500_2_scenario,u500_2_historical,v850_scenario,v850_historical]:
		'''
		Revise Cell Method
			The descrepency of cell method is caused by cf python versions
		'''
		if fld.cell_method("cellmethod1").get_axes()!=("area",): fld.cell_method("cellmethod1").set_axes("area")
		'''
		Remove Cell Measure
			Caused by missing of external files
		'''
		if fld.has_construct(identity='measure:area'): fld.del_construct(identity='measure:area')
	t850=cf.aggregate([t850_historical,t850_scenario],axes='T',donotchecknonaggregatingaxes=True)
	t850=check_aggregation(t850)
	t250=cf.aggregate([t250_historical,t250_scenario],axes='T',donotchecknonaggregatingaxes=True)
	t250=check_aggregation(t250)
	u500_1=cf.aggregate([u500_1_historical,u500_1_scenario],axes='T',donotchecknonaggregatingaxes=True)
	u500_1=check_aggregation(u500_1)
	u500_2=cf.aggregate([u500_2_historical,u500_2_scenario],axes='T',donotchecknonaggregatingaxes=True)
	u500_2=check_aggregation(u500_2)
	v850=cf.aggregate([v850_historical,v850_scenario],axes='T',donotchecknonaggregatingaxes=True)
	v850=check_aggregation(v850)
	'''
	calculate anomalies
	'''
	t850=t850-t850.collapse("mean","T")
	t250=t250-t250.collapse("mean","T")
	u500_1=u500_1-u500_1.collapse("mean","T")
	u500_2=u500_2-u500_2.collapse("mean","T")
	v850=v850-v850.collapse("mean","T")
	'''
	calculate shears
	'''
	delta_t=t850-t250
	delta_u=u500_2-u500_1
	'''
	standardise anomalies
	'''
	delta_t=(delta_t-delta_t.collapse('mean','T'))/delta_t.collapse('sd','T')
	delta_u=(delta_u-delta_u.collapse('mean','T'))/delta_u.collapse('sd','T')
	v850=(v850-v850.collapse('mean','T'))/v850.collapse('sd','T')
	'''
	assemble
	'''
	hwi=delta_t+delta_u+v850
	'''
	standardise hwi
	'''
	hwi=(hwi-hwi.collapse('T: mean'))/hwi.collapse('T: sd')
	hwi.set_property('standard_name','hwi')

	return hwi,delta_t,delta_u,v850
'''
Calculate Wang-Chen index
'''
def cal_wci(sib,al,mar):
	sib_cc=(sib-sib.collapse('T: mean'))/sib.collapse('T: sd')
	al_cc=(al-al.collapse('T: mean'))/al.collapse('T: sd')
	mar_cc=(mar-mar.collapse('T: mean'))/mar.collapse('T: sd')

	wci=.5*(2.*sib_cc-al_cc-mar_cc)
	'''
	This is an extra step to bring comparison equal
	'''
	wci=(wci-wci.collapse('T: mean'))/wci.collapse('T: sd')

	wci.set_property("standard_name","wci")
	return wci,sib_cc,al_cc,mar_cc

def cal_wci_scenario(sib_scenario,sib_historical,al_scenario,al_historical,mar_scenario,mar_historical):
	sib=cf.aggregate([sib_historical,sib_scenario],axes='T',donotchecknonaggregatingaxes=True)
	sib=check_aggregation(sib)
	al=cf.aggregate([al_historical,al_scenario],axes='T',donotchecknonaggregatingaxes=True)
	al=check_aggregation(al)
	mar=cf.aggregate([mar_historical,mar_scenario],axes='T',donotchecknonaggregatingaxes=True)
	mar=check_aggregation(mar)

	sib_cc=(sib-sib.collapse('T: mean'))/sib.collapse('T: sd')
	al_cc=(al-al.collapse('T: mean'))/al.collapse('T: sd')
	mar_cc=(mar-mar.collapse('T: mean'))/mar.collapse('T: sd')

	wci=.5*(2.*sib_cc-al_cc-mar_cc)
	'''
	This is an extra step to bring comparison equal
	'''
	wci=(wci-wci.collapse('T: mean'))/wci.collapse('T: sd')

	wci.set_property("standard_name","wci")
	return wci,sib_cc,al_cc,mar_cc
'''
Sort sort a cf.FieldList, calculate model mean from ensemble, return a list of cf.Field
INPUT :
	fldlst : A cf.FieldList with cf.Field entries of model ensemble members.
OUTPUT :
	fldmm : A cf.Field with with an extra dimemtion, sample, representing model mean.
'''
def fldlst2lstmm(fldlst) :

	fldmm=[]; fldlst_names=[]

	# List fld names
	for fld in fldlst:
		fldlst_names.append(fld.get_property("source_id"))

	# Find unique model names from fldlst_new
	for i,fld_name in enumerate(fldlst_names):
		if i==0:
			print(fld_name)
			model_names=[fld_name]
		else:
			if fld_name!=fldlst_names[i-1]:
				print(fld_name)
				if (fld_name in model_names):
					print("{} exists!".format(fld_name))
				else:
					model_names.append(fld_name)

	# Calculate model mean, and put the mean into a list
	fldmm=[]
	for model_name in model_names:

		# Locate flds of the same source_id
		fld_idx=[i for i,fld_name in enumerate(fldlst_names) if fld_name==model_name]
		print(model_name,fld_idx)

		# Subset fldlst according to fld_idx
		fldlst_model=[]
		for i in fld_idx:
			fldlst_model.append(fldlst[i])

		# Model mean
		fldmm.append(sum(fldlst_model)/len(fldlst_model))

	return fldmm
'''
New function:
	- Only to deal single-value cf.Field;
	- Remove masked cf.Field;
	- Apply sum() to calculate List mean;
'''
def fldlst2lstmm_new(fldlst):
	'''
	Rebuild fldlst to contain only non missing value flds
	'''
	fldlst_new=[];fldlst_names=[]
	for i,fld in enumerate(fldlst):
		if fld.mask.array[0]:
			print("... Remove {}/{}: {}...".format(fld.get_property("source_id"),fld.get_property("variant_label"),fld.array[0]))
		else:
			fldlst_new.append(fld)
			fldlst_names.append(fld.get_property("source_id"))
	'''
	Find unique model names from fldlst_new
	'''
	for i,fld_name in enumerate(fldlst_names):
		if i==0:
			print(fld_name)
			model_names=[fld_name]
		else:
			if fld_name!=fldlst_names[i-1]:
				print(fld_name)
				if (fld_name in model_names):
					print("{} exists!".format(fld_name))
				else:
					model_names.append(fld_name)
	'''
	Model mean
	'''
	fldmm=[]
	for model_name in model_names:
		fld_idx=[i for i,fld_name in enumerate(fldlst_names) if fld_name==model_name]
		print(model_name,fld_idx)
		'''
		Subset fldlst according to fld_idx
		'''
		fldlst_model=[]
		for i in fld_idx:
			fldlst_model.append(fldlst_new[i])
		'''
		Model mean
		'''
		fldmm.append(sum(fldlst_model)/len(fldlst_model))
	return fldmm

'''
pack_field_2d:
	pack a 2d array into a cf.field
'''
def pack_field_2d(fld,data,**keys):
#	fld_dummy=cf.Field(properties=fld.properties())
#	for key,axis in fld.domain_axes.items():
#		print(key,axis)
#		fld_dummy.set_construct(axis,key=key)
#	fld_dummy.set_data(data,axes=fld.get_data_axes())
#	pdb.set_trace()
	fld_dummy=fld.copy()
#	fld_dummy.data=data
	fld_dummy.set_data(data,fld.get_data_axes())
	if 'standard_name' in keys:
		fld_dummy.set_property('standard_name',keys['standard_name'])
	return fld_dummy
'''
Inputs:
	x: 1d vector
	y: 3d cf.field
'''
def linregress_2d(x,var):
	slope=np.zeros(var.shape[1:])
	rvalue=np.zeros(var.shape[1:])
	pvalue=np.zeros(var.shape[1:])
	stderr=np.zeros(var.shape[1:])


	for j in range(var.shape[1]):
		for i in range(var.shape[2]):
			if not var[0,j,i].squeeze().mask.array:
				slope[j,i],interp,rvalue[j,i],pvalue[j,i],stderr[j,i]=stats.linregress(x,var[:,j,i].squeeze().array)
	
	slope_fld=pack_field_2d(var[0].squeeze(),cf.Data(slope,'1'),standard_name='slope')
	rvalue_fld=pack_field_2d(var[0].squeeze(),cf.Data(rvalue,'1'),standard_name='rvalue')
	pvalue_fld=pack_field_2d(var[0].squeeze(),cf.Data(pvalue,'1'),standard_name='pvalue')
	stderr_fld=pack_field_2d(var[0].squeeze(),cf.Data(stderr,var.units),standard_name='stderr')

	return slope_fld,rvalue_fld,pvalue_fld,stderr_fld
'''
match_list:
	1. Read in two file lists;
	2. Matches two lists according to the model name and the ensemble member;
Input: either the whole filename or model name and ensemble member;
Output: matched model name (lst_model) and ensemble member (lst_member);
'''
def match_list(**keys):
	if 'lst0' in keys:
		lst0 = keys['lst0']
		lst0_model = []
		lst0_member = []
		for i in range(len(lst0)):
			lst0_model.append(lst0[i].split("_")[2])
			lst0_member.append(lst0[i].split("_")[4])
	else:
		lst0_model = keys['lst0_model']
		lst0_member = keys['lst0_member']

	if 'lst1' in keys:
		lst1 = keys['lst1']
		lst1_model = []
		lst1_member = []
		for i in range(len(lst1)):
			lst1_model.append(lst1[i].split("_")[2])
			lst1_member.append(lst1[i].split("_")[4])
	else:
		lst1_model = keys['lst1_model']
		lst1_member = keys['lst1_member']

	lst_model = []
	lst_member = []
	for i in range(len(lst0_model)):
		for j in range(len(lst1_model)):
			if lst0_model[i] == lst1_model[j] and lst0_member[i] == lst1_member[j]:
				print(i)
				lst_model.append(lst0_model[i])
				lst_member.append(lst0_member[i])
	
	return lst_model,lst_member
'''
hwi_anomaly:
	calculate anomaly against historical from three periods: 2025-2034, 2035-2044 and 2045-2054
Input:
	scenario_flst: field list of hwi mm during 2015 to 2100;
	cmip_flst: field list of hwi mm during 1979 to 2014;
'''
def hwi_anomaly(scenario_flst,cmip_flst):
	hwi_ano=[[] for i in range(3)]
	for fld in scenario_flst:
		source_id=fld.get_property('source_id')
		try:
			cmip_mean=cmip_flst.select_field('source_id={}'.format(source_id)).collapse('T: mean')
			hwi_ano[0].append(fld[30:60].collapse('T: mean')-cmip_mean)
			hwi_ano[1].append(fld[61:90].collapse('T: mean')-cmip_mean)
			hwi_ano[2].append(fld[91:120].collapse('T: mean')-cmip_mean)
		except:
			print('{} cannot be found in Historical!'.format(source_id))
			pass
	return hwi_ano

def hwi_anomaly2(scenario_flst,cmip_flst,*,flip:bool=False,half:bool=False):
	if flip:
		w_flip=-1.
	else:
		w_flip=1.
	if half:
		w_half=.5
	else:
		w_half=1.


	mdllst_cmip=cmip_flst.coordinates("source_id").value().array[:]

	hwi_ano=[[] for i in range(4)]
	for fld in scenario_flst:
		mdlnm=fld.coordinates("source_id").value().array[0]
		if mdlnm in mdllst_cmip:
			cmip_mean=cmip_flst.subspace(source_id=mdlnm).collapse('T: mean')
			hwi_ano[0].append((fld[0,30:60].collapse('T: mean')-cmip_mean)*w_flip*w_half)
			hwi_ano[1].append((fld[0,60:90].collapse('T: mean')-cmip_mean)*w_flip*w_half)
			hwi_ano[2].append((fld[0,90:120].collapse('T: mean')-cmip_mean)*w_flip*w_half)
			hwi_ano[3].append((fld[0,225:255].collapse('T: mean')-cmip_mean)*w_flip*w_half)
		else:
			print('{} cannot be found in Historical!'.format(mdlnm))
			pass
	return hwi_ano

def hwi_anomaly3(scenario_flst,*,flip:bool=False,half:bool=False):
	if flip:
		w_flip=-1.
	else:
		w_flip=1.
	if half:
		w_half=.5
	else:
		w_half=1.

	hwi_ano=[[] for i in range(4)];hwi_snr=[[] for i in range(4)];hwi_his=[[] for i in range(4)];
	for fld in scenario_flst:
		hwi_historical=fld.subspace(T=cf.wi(cf.dt("1979-01-01"),cf.dt("2014-12-30"))).collapse('T: mean')

		hwi_ano[0].append((fld.subspace(T=cf.wi(cf.dt("2025-01-01"),cf.dt("2034-12-30"))).collapse('T: mean')-hwi_historical)*w_flip*w_half)
		hwi_ano[1].append((fld.subspace(T=cf.wi(cf.dt("2035-01-01"),cf.dt("2044-12-30"))).collapse('T: mean')-hwi_historical)*w_flip*w_half)
		hwi_ano[2].append((fld.subspace(T=cf.wi(cf.dt("2045-01-01"),cf.dt("2054-12-30"))).collapse('T: mean')-hwi_historical)*w_flip*w_half)
		hwi_ano[3].append((fld.subspace(T=cf.wi(cf.dt("2090-01-01"),cf.dt("2099-12-30"))).collapse('T: mean')-hwi_historical)*w_flip*w_half)

		hwi_snr[0].append(fld.subspace(T=cf.wi(cf.dt("2025-01-01"),cf.dt("2034-12-30"))).collapse('T: mean')*w_flip*w_half)
		hwi_snr[1].append(fld.subspace(T=cf.wi(cf.dt("2035-01-01"),cf.dt("2044-12-30"))).collapse('T: mean')*w_flip*w_half)
		hwi_snr[2].append(fld.subspace(T=cf.wi(cf.dt("2045-01-01"),cf.dt("2054-12-30"))).collapse('T: mean')*w_flip*w_half)
		hwi_snr[3].append(fld.subspace(T=cf.wi(cf.dt("2090-01-01"),cf.dt("2099-12-30"))).collapse('T: mean')*w_flip*w_half)

		hwi_his[0].append(hwi_historical*w_flip*w_half)
		hwi_his[1].append(hwi_historical*w_flip*w_half)
		hwi_his[2].append(hwi_historical*w_flip*w_half)
		hwi_his[3].append(hwi_historical*w_flip*w_half)

	return hwi_ano,hwi_snr,hwi_his
'''
fld2array:
	2D list of field; stript away metadata; keep data.
'''
def fldlst2list_2level(fldlst2d):
	hwi2d=[[] for i in range(len(fldlst2d))]
	for i,fldlst1d in enumerate(fldlst2d):
		for fld in fldlst1d:
			try:
				print(i,fld.coordinates("source_id").value().array[0],fld.array[0])
				hwi2d[i].append(fld.array[0,0])
			except ValueError:
				print(i,fld.get_property('source_id'),fld.squeeze().array)
				hwi2d[i].append(fld.array[0,0])
	return hwi2d
'''
aggregate_fieldlist:
	1. Copy T coordinate from the first field of the list (potnetially, this need to be fixed);
	2. Rebuild fieldlist with original properties and data, but identical T coordinate from the first field;
	3. Aggregate the new fieldlist with cf.aggregate;
	4. Build auxilary coordinates: experiment, institution, source, variant;
	5. Name the new aggregate coordinate as 'sample';
Input: a fieldlist
Output: an aggregated field
'''
#def aggregate_fieldlist(flst):
#	'''
#	keep t dimensional properties from the first field
#	'''
#	T0=flst[0].coord('T')
#	T0_axes=flst[0].get_data_axes('T')
#	'''
#	create a new fieldlist
#	'''
#	flst_cp=[];experiment_info=[];institution_info=[];source_info=[];variant_info=[]
#	for i,x in enumerate(flst):
#		xx=cf.Field(properties=x.properties())
#		'''
#		create domain axes explicitly
#		'''
#		for key,axis in x.domain_axes.items():
#			xx.set_construct(axis,key=key)
#		xx.set_data(x.data,axes=flst[0].get_data_axes())
#		xx.set_construct(T0,axes=T0_axes)
#		xx.set_property('sample',i)
#		xx.set_property('standard_name','precipitation_flux')
#		'''
#		Correct irregualte psl standard name
#		'''
#		if xx.get_property('standard_name')=='air_pressure_at_sea_level':
#			xx.set_property('standard_name','air_pressure_at_mean_sea_level')
#		flst_cp.append(xx)
#		institution_info.append(x.get_property('institution_id'))
#		experiment_info.append(x.get_property('experiment_id'))
#		source_info.append(x.get_property('source_id'))
#		variant_info.append(x.get_property('variant_label'))
#	'''
#	Since following aggregation code does not work, try mmm instead.
#	'''
#	mmm=0
#	for x in flst_cp:
#		mmm=mmm+x
#	mmm=mmm/float(len(flst_cp))
#
#	return mmm
#'''
#pack_field_2d:
#	pack a 2d array into a cf.field
#'''
#def pack_field_2d(fld,data,**keys):
#	fld_dummy=cf.Field(properties=fld.properties())
#
#	for key,axis in fld.domain_axes.items():
#		print(key,axis)
#		fld_dummy.set_construct(axis,key=key)
#	pdb.set_trace()
#	fld_dummy.set_construct(fld.coord('Y'),axes=fld.get_data_axes('Y'))
#	fld_dummy.set_construct(fld.coord('X'),axes=fld.get_data_axes('X'))
#	fld_dummy.set_data(data,fld.get_data_axes())
#	if 'standard_name' in keys:
#		fld_dummy.set_property('standard_name',keys['standard_name'])
#	return fld_dummy

'''
http://geog.uoregon.edu/GeogR/topics/ttest.pdf
ttest4mean_1sample:
	For one sample ttest
Input:
	mean, variance and size of sample for both samples
Output:
	p value packed into a field
'''
def ttest4mean_1sample(mean0, base, sd0, n0):

	n0 = float(n0)

	if isinstance(mean0, cf.field.Field) :
		tvalue = (mean0.array - base.array) / sd0
		df = n0-1
		pvalue = (1.-stats.t.cdf(abs(tvalue),df)) * 2.
		pvalue_fld = pack_field_2d(mean0, cf.Data(pvalue,'1'), standard_name='pvalue')

	if isinstance(mean0, float):
		tvalue = (mean0 - base) / sd0
		df = n0 - 1
		pvalue = (1.-stats.t.cdf(abs(tvalue),df)) * 2.
		pvalue_fld = pvalue

	return pvalue_fld
'''
ttest4mean:
	For two samples ttest with population variances are assumed to be equal.
Input:
	mean, variance and size of sample for both samples
Output:
	p value packed into a field
'''
def ttest4mean(mean0,mean1,sd0,sd1,n0,n1):
	n0=float(n0)
	n1=float(n1)
	'''
	For cf.Field
	'''
	if type(mean0)=="cf.field.Field":
		pooled_variance=(sd0.array*sd0.array*(n0-1)+sd1.array*sd1.array*(n1-1))/(n0+n1-2)
		variance=np.sqrt(pooled_variance)*np.sqrt(1/n0+1/n1)
		tvalue=(mean0.array-mean1.array)/variance
		df=n0+n1-2
		pvalue=(1.-stats.t.cdf(abs(tvalue),df))*2.
		pvalue_fld=pack_field_2d(mean0,cf.Data(pvalue,'1'),standard_name='pvalue')
	'''
	For a pair of value
	'''
	if isinstance(mean0,float):
		pooled_variance=(sd0*sd0*(n0-1)+sd1*sd1*(n1-1))/(n0+n1-2)
		variance=np.sqrt(pooled_variance)*np.sqrt(1/n0+1/n1)
		tvalue=(mean0-mean1)/variance
		df=n0+n1-2
		pvalue_fld=(1.-stats.t.cdf(abs(tvalue),df))*2.

	return pvalue_fld
'''
https://www.statology.org/f-test-python/
'''
def ftest(std0,std1,n0,n1):
	fvalue=(std0*std0)/(std1*std1)
	df0=n0-1;df1=n1-1;
	pvalue=1-stats.f.cdf(fvalue,df0,df1)
	return pvalue
'''
Match historcal to scemario
'''
def matched_model_mean(fld_historical,fld_scenario):
	mdlnm_lst_historical=[]
	for fld in fld_historical:
		mdlnm_lst_historical.append(fld.aux("source_info").array.tolist())

	mdlnm_lst_scenario=[]
	for fld in fld_scenario:
		mdlnm_lst_scenario.append(fld.aux("source_info").array.tolist())

	'''
	Match historical with Scenario
	'''
	idx_historical_matched=[]
	for i,mdlnm in enumerate(mdlnm_lst_historical):
		if mdlnm in mdlnm_lst_scenario:
			print(i,mdlnm)
			idx_historical_matched.append(i)
		
	var_scenario_mean=fld_scenario.collapse('mean','sample').squeeze()
	var_scenario_sd  =fld_scenario.collapse('sd'  ,'sample').squeeze()
	nsample=fld_scenario.shape[0]
		
	var_historical_mean=fld_historical[idx_historical_matched].collapse('mean','sample').squeeze()
	var_historical_sd  =fld_historical[idx_historical_matched].collapse('sd'  ,'sample').squeeze()

	return var_historical_mean,var_historical_sd,var_scenario_mean,var_scenario_sd,nsample

def matched_model_median(fld_historical,fld_scenario):
	mdlnm_lst_historical=[]
	for fld in fld_historical:
		mdlnm_lst_historical.append(fld.aux("source_info").array.tolist())

	mdlnm_lst_scenario=[]
	for fld in fld_scenario:
		mdlnm_lst_scenario.append(fld.aux("source_info").array.tolist())

	'''
	Match historical with Scenario
	'''
	idx_historical_matched=[]
	for i,mdlnm in enumerate(mdlnm_lst_historical):
		if mdlnm in mdlnm_lst_scenario:
			print(i,mdlnm)
			idx_historical_matched.append(i)
		
	var_diff=fld_scenario.copy()
	for i,fld in enumerate(fld_historical[idx_historical_matched]):
		var_diff[i]=fld_scenario[i]-fld
	var_diff_median=var_diff.collapse('median','sample').squeeze()
#	var_scenario_mean=fld_scenario.collapse('mean','sample').squeeze()
#	var_scenario_sd  =fld_scenario.collapse('sd'  ,'sample').squeeze()
#	nsample=fld_scenario.shape[0]
		
#	var_historical_mean=fld_historical[idx_historical_matched].collapse('mean','sample').squeeze()
#	var_historical_sd  =fld_historical[idx_historical_matched].collapse('sd'  ,'sample').squeeze()

	return var_diff_median

'''
aggregate_fieldlist_nD :
	
INPUT :
	flst : cf.FieldList of model mean.
OUTPUT :
	flst_o : cf.Field of model mean, aggregated along coordinate "sample".
'''
def aggregate_fieldlist_nD(flst):

	# Initiate lists for cf.FieldList and auxiliary coordinates
	flst_cp = []; experiment_info = []; institution_info = []; source_info = []; variant_info = []

	for i,x in enumerate(flst):
		xs = x.squeeze()


		# delete degenerated vertical coordiantes
		if xs.get_property("standard_name") == "10m_wind_speed": xs.set_property("variable_id", "uv10")
		if xs.get_property("variable_id") == "uas" or xs.get_property("variable_id") == "vas": xs.del_construct("height")

		xx = cf.Field(properties=xs.properties())

		# Asign cf.DomainAxis
		for size in xs.shape:
			xx.set_construct(cf.DomainAxis(size=size))

		# Insert coordinate
		for coord_name in xs.dimension_coordinates():
			xx.set_construct(flst[0].construct(coord_name))

		# Insert field data
		xx.set_data(xs.data)

		# Set new attribute for aggregation on the next step.
		xx.set_property('sample', i)

		# Correct irregualte psl standard name
		if xx.get_property('standard_name') == 'air_pressure_at_sea_level':
			xx.set_property('standard_name', 'air_pressure_at_mean_sea_level')

		# Gather cf.Field into a list
		flst_cp.append(xx)

		# Cather lists of information for auxiliary coordinates
		institution_info.append(x.get_property('institution_id')); experiment_info.append(x.get_property('experiment_id')); source_info.append(x.get_property('source_id')); variant_info.append(x.get_property('variant_label'))

	# Perform aggregation on property "sample"
	flst_o = cf.aggregate(flst_cp, dimension=('sample',))[0]

	# Convert auxiliary coordinate and attach onto cf.Field
	experiment_info_id = cf.AuxiliaryCoordinate(properties=dict(standard_name='experiment_info'), data=cf.Data(experiment_info)); flst_o.set_construct(experiment_info_id)
	institution_info_id = cf.AuxiliaryCoordinate(properties=dict(standard_name='institution_info'), data=cf.Data(institution_info)); flst_o.set_construct(institution_info_id)
	source_info_id = cf.AuxiliaryCoordinate(properties=dict(standard_name='source_info'), data=cf.Data(source_info)); flst_o.set_construct(source_info_id)
	variant_info_id = cf.AuxiliaryCoordinate(properties=dict(standard_name='variant_info'), data=cf.Data(variant_info)); flst_o.set_construct(variant_info_id)

	# Build coordinate "sample" and attach to cf.Field
	sample_id = cf.DimensionCoordinate(properties={'standard_name':'sample'}, data=cf.Data(np.arange(len(flst)), '1')); flst_o.set_construct(sample_id)

	return flst_o
'''
For 1-D array
'''
def aggregate_fieldlist(flst):

	# Unify time coordinate using CanESM5
	for fld in flst:
		if fld.get_property("source_id") == "CanESM5":
			fld_base = fld

	flst_cp=[]; experiment_info=[]; institution_info=[]; source_info=[]; variant_info=[];
	for i,x in enumerate(flst):
		print(x.get_property("source_id"), x.get_property("variant_label"))

		# Remove degenerated coordinates, if there is any.
		xs = x.squeeze() 

		# Build new field, starting from copying attributes.
		xx = cf.Field(properties=xs.properties())

		#Assign cf.DomainAxis. Check domain size: cf.Field.domain
		for size in xs.shape:
			xx.set_construct(cf.DomainAxis(size=size))

		# Build cf.Coordinates
		for ic,coordinate_size in enumerate(fld_base.shape):
			# Only set the non-degenerated coordiantes
			if coordinate_size > 1:
				xx.set_construct(fld_base.coordinate(identity="dimensioncoordinate{}".format(ic)))

		# Asign data
		xx.set_data(xs.data)

		# Asign more attributes
		xx.set_property('sample', i)
		'''
		Correct irregualte psl standard name
		'''
		if xx.get_property('standard_name') == 'air_pressure_at_sea_level':
			xx.set_property('standard_name','air_pressure_at_mean_sea_level')
		if xx.get_property('standard_name') == 'atmosphere_optical_thickness_due_to_dust_ambient_aerosol':
			xx.set_property('standard_name','atmosphere_optical_thickness_due_to_dust_ambient_aerosol_particles')

		flst_cp.append(xx);institution_info.append(x.get_property('institution_id'));experiment_info.append(x.get_property('experiment_id'));source_info.append(x.get_property('source_id'));variant_info.append(x.get_property('variant_label'))

	# Aggregate the list along property "sample", turn the list into a cf.Field
	flst_o = cf.aggregate(flst_cp,dimension=('sample',))[0]

	# Attach ansiliary properties on the newly built cf.Field
	experiment_info_id=cf.AuxiliaryCoordinate(properties=dict(standard_name='experiment_info'),data=cf.Data(experiment_info))
	flst_o.set_construct(experiment_info_id, axes="domainaxis1")
	institution_info_id=cf.AuxiliaryCoordinate(properties=dict(standard_name='institution_info'),data=cf.Data(institution_info))
	flst_o.set_construct(institution_info_id, axes="domainaxis1")
	source_info_id=cf.AuxiliaryCoordinate(properties=dict(standard_name='source_info'),data=cf.Data(source_info))
	flst_o.set_construct(source_info_id, axes="domainaxis1")
	variant_info_id=cf.AuxiliaryCoordinate(properties=dict(standard_name='variant_info'),data=cf.Data(variant_info))
	flst_o.set_construct(variant_info_id, axes="domainaxis1")

	# Assign cooridinate "sample" and value to the prepended cooridinate "sample"
	sample_id=cf.DimensionCoordinate(properties={'standard_name':'sample'},data=cf.Data(np.arange(len(flst)),'1'))
	flst_o.set_construct(sample_id, axes="domainaxis1")

	return flst_o
'''
Calculate standard deviation of a cf.FieldList
	Each cf.Field needs to be single value;
'''
def stdev_fldlst(fldlst):
	'''
	New list stripted of cf.Field wrapping
	'''
	lst=[]

	for fld in fldlst:
		'''
		Check cf.Field size
		'''
		if fld.size!=1:
			print("Only single value field can be processed ...")
			exit()
		else:
		#	if fld.get_property("source_id")!="CanESM5": lst.append(fld.array[0])
			lst.append(fld.array[0])
	'''
	Calculate standard deviation
		To avoid missing values, list is first converted into numpy array, then the standard deviation is calculated using np.nanstd()
	'''
	lst=np.array(lst)
	lst_std=np.nanstd(lst)

	return lst_std

def quartile_fldlst(fldlst):
	'''
	New list stripted of cf.Field wrapping
	'''
	lst=[]

	for fld in fldlst:
		'''
		Check cf.Field size
		'''
		if fld.size!=1:
			print("Only single value field can be processed ...")
			exit()
		else:
			lst.append(fld.array[0])
	'''
	Calculate the 1st and 3rd quartiles:
	'''
	lst=np.array(lst)
	lst_quartile=np.nanpercentile(lst,[25,75])
	lst_median=np.nanmedian(lst)
	lst_quartile[0]=np.absolute(lst_quartile[0]-lst_median)
	lst_quartile[1]=np.absolute(lst_quartile[1]-lst_median)

	return lst_quartile

def mean_fldlst(fldlst):
	'''
	New list stripted of cf.Field wrapping
	'''
	lst=[]

	for fld in fldlst:
		'''
		Check cf.Field size
		'''
		if fld.size!=1:
			print("Only single value field can be processed ...")
			exit()
		else:
			lst.append(fld.array[0])
	'''
	Calculate standard deviation
		To avoid missing values, list is first converted into numpy array, then the standard deviation is calculated using np.nanstd()
	'''
	lst=np.array(lst)
	lst_mean=np.nanmean(lst)

	return lst_mean

def median_fldlst(fldlst):
	'''
	New list stripted of cf.Field wrapping
	'''
	lst=[]

	for fld in fldlst:
		'''
		Check cf.Field size
		'''
		if fld.size!=1:
			print("Only single value field can be processed ...")
			exit()
		else:
			lst.append(fld.array[0])
	'''
	Calculate standard deviation
		To avoid missing values, list is first converted into numpy array, then the standard deviation is calculated using np.nanstd()
	'''
	lst=np.array(lst)
	lst_median=np.nanmedian(lst)

	return lst_median

def tolist_fldlst(fldlst):
	'''
	New list stripted of cf.Field wrapping
	'''
	lst=[]

	for fld in fldlst:
		'''
		Check cf.Field size
		'''
		if fld.size!=1:
			print("Only single value field can be processed ...")
			exit()
		else:
			lst.append(fld.array[0])
	lst=np.array(lst)

	return lst
'''
Check depth of a dictionary, an elegant solution from https://stackoverflow.com/questions/23499017/know-the-depth-of-a-dictionary
'''
def depth(d):
	if isinstance(d,dict):
		return 1+(max(map(depth,d.values())) if d else 0)
	return 0
