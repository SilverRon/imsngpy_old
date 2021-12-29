#!/usr/bin/env python3
#============================================================
#	IMSNG Pipeline
#	=> Processing
#	Data Monitoring => Processing => Transient Search
#============================================================
#%%
#	Library
#------------------------------------------------------------
import os
import glob
import numpy as np
import warnings
warnings.filterwarnings(action='ignore')
import time
import sys
sys.path.append('/home/paek/imsngpy')
#	IMSNGpy modules
from tableutil import getccdinfo
from preprocess import *
from misc import *
from phot import *
from util import *
#
start_localtime = time.strftime('%Y-%m-%d %H:%M:%S (%Z)', time.localtime())
#	Astropy
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from ccdproc import ImageFileCollection
from astropy.time import Time
from astropy.nddata import CCDData
#	Multiprocess tools
from itertools import repeat
import multiprocessing
#------------------------------------------------------------
#	My library
# from tableutil import *
#============================================================
#	USER SETTING
#============================================================
#	Input
#------------------------------------------------------------
"""
#	[0] Folder to process
try:
	path_raw = sys.argv[1]
except:
	path_raw = input('''# Data folder to process : ''')
#	[1]	Observatory_ccd
try:
	obs = (sys.argv[2]).upper()
except:
	obs = input('''# Observatory(_ccd) to run
--------------------
LOAO
DOAO
SOAO
CBNUO
KCT_ASI1600MM
KCT_STX16803
KHAO
RASA36
LSGT
---------------------
:''').upper()
print('# Observatory : {}'.format(obs.upper()))
#	[3]	The number of cores
try:
	ncore = int(sys.argv[3])
except:
	ncore = 8
"""
#	Test setting
# path_raw = '/data6/obsdata/LOAO/1994_1026'
# path_raw = '/data6/obsdata/LOAO/1994_1003'
# path_raw = '/data6/obsdata/LOAO/1969_0119'
path_raw = '/data6/obsdata/LOAO/test'
obs = 'LOAO'
fast_mode4mframe = True
ncore = 8
#------------------------------------------------------------
#	PATH
#------------------------------------------------------------
path_factory = '/data3/paek/factory'
path_gal = '/data6/IMSNG/IMSNGgalaxies'
# path_config = '/home/paek/config'
path_config = '/home/paek/imsngpy/config'
path_log = '/home/paek/log'
path_bkg = '/data6/bkgdata'
path_table = '/home/paek/imsngpy/table'
path_gppy = '/home/paek/imsngpy'
#------------------------------------------------------------
path_mframe = f'{path_factory}/master_frames'
path_ref = f'{path_factory}/ref_frames/{obs.upper()}'
path_obs = f'{path_factory}/{obs.lower()}'
path_default_gphot = f'{path_config}/gphot.{obs.lower()}.config'
#------------------------------------------------------------
path_save = f'{path_bkg}/{obs.upper()}'
#------------------------------------------------------------
ccddat = f'{path_table}/obs.dat'
#------------------------------------------------------------
#	Codes
path_phot = f'{path_gppy}/imsngpy/gpphot.py'
path_phot_sg = '/home/paek/qsopy/phot/gregoryphot_2021.py'
path_phot_mp = '/home/paek/qsopy/phot/gregoryphot_mp_2021.py'
path_phot_sub = '/home/paek/qsopy/phot/gregoryphot_sub_2021.py'
# path_find = '/home/paek/qsopy/phot/gregoryfind_2021.py'
# path_find = '/home/paek/qsopy/phot/gregoryfind_bulk_2021.py'
path_find = '/home/paek/qsopy/phot/gregoryfind_bulk_mp_2021.py'
#------------------------------------------------------------
#	Table
logtbl = ascii.read(f'{path_log}/{obs.lower()}.log')
hdrtbl = ascii.read(f'{path_table}/changehdr.dat')
alltbl = ascii.read(f'{path_table}/alltarget.dat')
frgtbl = ascii.read(f'{path_table}/fringe.dat')
# ccdtbl = ascii.read(f'{path_table}/ccd.dat') 
ccdtbl = ascii.read(f'{path_table}/ccd.tsv') 
#------------------------------------------------------------
path_data = f'{path_obs}/{os.path.basename(path_raw)}'
print(f"""{'-'*60}\n#\tCOPY DATA\n{'-'*60}""")
#	Remove former data
if os.path.exists(path_data):
	rmcom = f'rm -rf {path_data}'
	print(rmcom)
	os.system(rmcom)
#	Copy to factory director
cpcom = f'cp -r {path_raw} {path_data}'
print(cpcom)
os.system(cpcom)
#%%
#	Identify CCD
ic0 = ImageFileCollection(path_data, keywords='*')

'''
print(f"""{'-'*60}\n#\tIDENTIFY CCD\n{'-'*60}""")
ic0 = ImageFileCollection(path_data, keywords='*')
for key, val, suf, ccd in zip((ccdtbl['key'][ccdtbl['obs']==obs]), (ccdtbl['value'][ccdtbl['obs']==obs]), (ccdtbl['suffix'][ccdtbl['obs']==obs]), (ccdtbl['ccd'][ccdtbl['obs']==obs])):
	if (key.lower() in ic0.keywords) & (val == ic0.summary[key.lower()][0]):
		ccdkey = key
		ccdval = val
		ccdtype = ccd
		if suf.mask == True:
			#	No suffix
			suffix = ''
			obsccd = f'{obs}'
		else:
			suffix = suf
			obsccd = f'{obs}_{suffix}'
		print(f'OBSERVAT : {obs}\nCCD KEYWORD : {key}\nCCD HEADER VALUE : {val}\nCCD NAME : {ccdtype}\nSUFFIX : {suffix}\n==> OBS_CCD : {obsccd}')
'''
ccdkey, ccdval, ccdtype, obsccd = identify_ccdinfo(ic0, obs, ccdtbl)
#	CCD INFO
indx_ccd = np.where(
	(ccdtbl['obs']==obs) &
	(ccdtbl['key']==ccdkey) &
	(ccdtbl['value']==ccdval)
)
print(f"""{'-'*60}\n#\tCCD INFO\n{'-'*60}""")
gain = ccdtbl['gain'][indx_ccd][0]*(u.electron/u.adu)
rdnoise = ccdtbl['readnoise'][indx_ccd][0]*(u.electron)
pixscale = ccdtbl['pixelscale'][indx_ccd][0]*(u.arcsec/u.pixel)
fov = ccdtbl['foveff'][indx_ccd][0]*(u.arcmin)
print(f"""GAIN : {gain}\nREAD NOISE : {rdnoise}\nPIXEL SCALE : {pixscale}\nEffective FoV : {fov}""")
#------------------------------------------------------------
#%%
#	Header correction
#------------------------------------------------------------
comment = f"""{'='*60}\n#\tHEADER CORRECTION\n{'='*60}"""
print(comment)

for i, inim in enumerate(ic0.summary['file']):
	image = f'{path_data}/{inim}'
	#	CCD Type
	'''
	if i == 0:
		for key in set(ccdtbl['key']):
			if key.lower() in ic0.summary.keys():
				ccdtype = ccdtbl[
					(ccdtbl['value'] == ic0.summary[key.lower()][i]) &
					(ccdtbl['obs'] == obs)
				]['ccd'].item()
			else:
				ccdtype = 'UNKNOWN'
	fits.setval(image, 'CCDNAME', value=ccdtype)
	fits.setval(image, 'OBSERVAT', value=obs)'''
	fits.setval(image, 'PATHPRC', value=path_data, comment='Path where data is processed')
	fits.setval(image, 'CCDNAME', value=ccdtype)
	fits.setval(image, 'OBSERVAT', value=obs)
	fits.setval(image, 'OBSCCD', value=obsccd)
	#	Correction with table
	for key, val, nval in zip(hdrtbl['key'], hdrtbl['val'], hdrtbl['newval']):
		if ic0.summary[key.lower()][i] == val:
			print(f'{inim} - {key} : {val} --> {nval}')
			fits.setval(image, key, value=nval)
	#	DATE-OBS, JD, MJD
	if 'T' not in ic0.summary['date-obs'][i]:
		dateobs = f"{ic0.summary['date-obs']}'T'{ic0.summary['time-obs']}"
		fits.setval(image, 'date-obs', value=dateobs)
		del dateobs
	else:
		pass
	t = Time(ic0.summary['date-obs'][i], format='isot')
	fits.setval(image, 'JD', value=t.jd)
	fits.setval(image, 'MJD', value=t.mjd)
	del t
	#	OBJECT name
	if 'ngc' in ic0.summary['object'][i]:
		objectname = ic0.summary['object'][i]
		while len(objectname)<7:
			head = objectname[0:3]
			tail = objectname[3:]
			tail = f'0{tail}'
			objectname = f'{head}{tail}'
		fits.setval(image, 'OBJECT', value=objectname.upper())
		del objectname
		del head
		del tail
	print()
ic1 = ImageFileCollection(path_data, keywords='*')
t_med = np.median(ic1.filter(imagetyp='object').summary['jd'])	#	[JD]
#------------------------------------------------------------
#%%
#	Object Master Table
#
#	Write the status of processing
#	Pointing the original filename and updated filename
#	Each dtype is set to 'strtype' variable
#------------------------------------------------------------
strtype = 'U300'
omtbl = Table()
objtypelist = []
for obj in ic1.filter(imagetyp='object').summary['object']:
	if obj in alltbl['obj']:
		objtypelist.append('IMSNG')
	elif 'GRB' in obj:
		objtypelist.append('GRB')
	elif ('GECKO' in obj) | ('GCK' in obj):
		objtypelist.append('GECKO')
	else:
		objtypelist.append('NON-IMSNG')

omtbl['objtype'] = objtypelist
omtbl['raw'] = [inim for inim in ic1.filter(imagetyp='object').files]
omtbl['now'] = omtbl['raw']
omtbl['preprocess'] = ''
omtbl['defringe'] = ''
omtbl['cosmic_ray_removal'] = ''
omtbl['astrometry'] = ''
omtbl['final'] = [f'{path_data}/{fnamechange(inim)}' for inim in ic1.filter(imagetyp='object').files]
for key in omtbl.keys(): omtbl[key] = omtbl[key].astype(strtype)
#============================================================
#%%
#	Pre-processing
#------------------------------------------------------------
print(f"""{'='*60}\n#\tPREPARE FOR PRE-PROCESS\n{'='*60}""")
mframe = dict()
#------------------------------------------------------------
#	BIAS
#------------------------------------------------------------
print(f"""{'-'*60}\n#\tBIAS\n{'-'*60}""")
if 'bias' in ic1.summary['imagetyp']:
	biaslist = ic1.filter(imagetyp='bias').files
	print(f"""{len(biaslist)} bias frames --> master bias""")
	biasdata, biasim = master_bias(biaslist)
	mframe['bias'] = biasdata
	del biasdata
	path_bias = f'{path_mframe}/{obs}/zero'
	cpcom = f"cp {biasim} {path_bias}"
	print(cpcom)
	# os.system(cpcom)
else:
	print('No bias frame. Borrow from previous data.')
	mftype = 'zero'
	biasim = get_nearest_bias(
		mftype=mftype,
		t_med=t_med,
		ccdkey=ccdkey,
		ccdval=ccdval,
		keyword=f'{path_mframe}/{obs}/{mftype}/????????-{mftype}.fits',
		keywords=[ccdkey.lower(), 'jd'],
		fast_mode4mframe=True,
		)
	print(f'Borrow {os.path.basename(biasim)}')
	mframe['bias'] = CCDData(fits.getdata(biasim), unit="adu", meta=fits.getheader(biasim))
	'''
	mftype = 'zero'
	if fast_mode4mframe == True:
		biaslist = np.array(
			[os.path.basename(inim) for inim in sorted(glob.glob(f'{path_mframe}/{obs}/{mftype}/????????-{mftype}.fits'))]
			)
		deltarr = np.array(
			[np.abs(date2jd(inim.split('-')[0]).jd-t_med) for inim in biaslist]
		)
		indx_bias = np.where(deltarr == np.min(deltarr))
		biasim = f"{path_mframe}/{obs}/{mftype}/{biaslist[indx_bias].item()}"
	else:
		ic_bias = ImageFileCollection(
			location=f'{path_mframe}/{obs}/{mftype}',
			keywords=[
				ccdkey.lower(),
				# 'date-obs',
				# 'imagetyp',
				'jd',
				# 'mjd',
			]
			)
		ic_bias_avail = ic_bias.summary[
			(ic_bias.summary['jd'].mask == False) &
			(ic_bias.summary[ccdkey.lower()]==ccdval)
			]
		biaslist = ic_bias_avail['file']
		deltarr = np.array(np.abs(ic_bias_avail['jd']-t_med))
		indx_bias = np.where(deltarr == np.min(deltarr))
		biasim = f"{path_mframe}/{obs}/{mftype}/{biaslist[indx_bias].item()}"
		del ic_bias_avail
	#	Delete variables
	del mftype
	del deltarr
	del indx_bias
	del biasim'''

#------------------------------------------------------------
#%%
#	DARK
#------------------------------------------------------------
print(f"""{'-'*60}\n#\tDARK\n{'-'*60}""")

darkframe = dict()
if 'dark' in ic1.summary['imagetyp']:
	darkexptime = np.array(list(set(ic1.filter(imagetyp='dark').summary['exptime'])))
	for exptime in darkexptime:
		darkdata, darkim = master_dark(ic1.filter(imagetyp='dark', exptime=exptime).files, mbias=mframe['bias'])
		# darkframe[f'{int(exptime)}'] = master_dark(ic1.filter(imagetyp='dark', exptime=exptime).files, mbias=mframe['bias'])
		darkframe[f'{int(exptime)}'] = darkdata
		del darkdata
		path_dark = f'{path_mframe}/{obs}/dark'
		cpcom = f"cp {darkim} {path_dark}"
		print(cpcom)
		# os.system(cpcom)
else:
	print('No dark frame. Borrow from previous data.')
	mftype = 'dark'
	darkexptime = []
	for exptime in set(ic1.filter(imagetyp='object').summary['exptime']):
		mftype = 'dark'
		darkim, nexptime = get_nearest_dark(
			t_med=t_med,
			keyword=f'{path_mframe}/{obs}/{mftype}/*-????????-{mftype}.fits',
			exptime=exptime,
			ccdkey=ccdkey,
			ccdval=ccdval,
			keywords=[ccdkey.lower(), 'jd', 'exptime'],
			fast_mode4mframe=True,
			)
		print(f'Borrow {os.path.basename(darkim)}')
		# darkframe[f'{int(nexptime)}'] = CCDData(fits.getdata(darkim), unit="adu", meta=fits.getheader(darkim))

		if str(nexptime) not in darkframe.keys():
			darkframe[f'{int(nexptime)}'] = CCDData(fits.getdata(darkim), unit="adu", meta=fits.getheader(darkim))
			darkexptime.append(int(nexptime))
		else:
			print(f'No available dark frame')
			pass
		darkexptime = np.array(darkexptime)

mframe['dark'] = darkframe
del darkframe
#------------------------------------------------------------
#%%
#	FLAT
#------------------------------------------------------------
print(f"""{'-'*60}\n#\tFLAT\n{'-'*60}""")
flatframe = dict()
#	Filter list
if 'object' in ic1.summary['imagetyp']:
	objfilterlist = list(ic1.filter(imagetyp='object').summary['filter'])
	print(f'OBJECT FILTERS : {list(set(objfilterlist))}')
else:
	objfilterlist = []
if 'flat' in ic1.summary['imagetyp']:
	flatfilterlist = list(ic1.filter(imagetyp='flat').summary['filter'])
	print(f'FLAT FILTERS : {list(set(flatfilterlist))}')
else:
	flatfilterlist = []
filterlist = set(objfilterlist+flatfilterlist)
print(f'==> ALL FILTERS : {list(set(filterlist))}')

if 'flat' in ic1.summary['imagetyp']:
	#	Dark exptime should be corrected!
	indx_mindark = np.where(darkexptime == np.min(darkexptime))
	mdark = mframe['dark'][str(int(darkexptime[indx_mindark].item()))]
	for filte in filterlist:
		# print(filte)
		flatdata, flatim = master_flat(ic1.filter(imagetyp='flat', filter=filte).files, mbias=mframe['bias'], mdark=mdark, filte=filte)
		# flatframe[filte] = master_flat(ic1.filter(imagetyp='flat', filter=filte).files, mbias=mframe['bias'], mdark=mdark, filte=filte)
		flatframe[filte] = flatdata
		path_flat = f'{path_mframe}/{obs}/flat'
		cpcom = f"cp {flatim} {path_flat}"
		print(cpcom)
		# os.system(cpcom)
		del flatdata
	del mdark
else:
	print('No Flat frame. Borrow from previous data.')
	mftype = 'flat'
	for filte in set(ic1.filter(imagetyp='object').summary['filter']):

		flatim = get_nearest_flat(
			t_med,
			keyword=f'{path_mframe}/{obs}/{mftype}/????????-n{filte}.fits',
			filte=filte,
			ccdkey=ccdkey,
			ccdval=ccdval,
			keywords=[ccdkey.lower(), 'imagetyp', 'jd', 'filter',],
			fast_mode4mframe=True,
			)

		print(f'Borrow {os.path.basename(flatim)}')
		flatframe[filte] = CCDData(fits.getdata(flatim), unit="adu", meta=fits.getheader(flatim))

mframe['flat'] = flatframe
del flatframe
#------------------------------------------------------------
#%%
#	OBJECT Correction
#------------------------------------------------------------
print(f"""{'-'*60}\n#\tOBJECT CORRECTION ({len(ic1.filter(imagetyp='object').files)})\n{'-'*60}""")
#	Function for multi-process
def obj_process4mp(inim, newim, filte, exptime, darkexptime, rdnoise, mframe,):
	'''
	Routine for multiprocess
	'''
	#	Find the closest exposuretime betw dark and object
	indx_closet = np.where(
		np.abs(exptime-darkexptime) == np.min(np.abs(exptime-darkexptime))
	)
	bestdarkexptime = darkexptime[indx_closet].item()
	print(f"{os.path.basename(inim)} {exptime} sec in {filte}-band <-- (scaled) DARK {int(darkexptime[indx_closet].item())} sec")
	#	Process
	nccd = obj_process(
		inim=inim,
		# gain=ccdinfo['gain'],
		gain=None,
		readnoise=rdnoise,
		mbias=mframe['bias'],
		mdark=mframe['dark'][str(int(bestdarkexptime))],
		mflat=mframe['flat'][filte],
	)
	# nccd.write(f'{os.path.dirname(inim)}/fdz{os.path.basename(inim)}', overwrite=True)
	nccd.write(newim, overwrite=True)
#	Run with multi-process
fdzimlist = add_prefix(ic1.filter(imagetyp='object').files, 'fdz')
if __name__ == '__main__':
	#	Fixed the number of cores (=4)
	with multiprocessing.Pool(processes=4) as pool:
		results = pool.starmap(
			obj_process4mp,
			zip(
				ic1.filter(imagetyp='object').files,
				fdzimlist,
				ic1.filter(imagetyp='object').summary['filter'],
				ic1.filter(imagetyp='object').summary['exptime'],
				repeat(darkexptime),
				repeat(rdnoise),
				repeat(mframe),
				)
			)
#	Image collection for pre-processed image
fdzic = ImageFileCollection(f'{path_data}', glob_include='fdzobj*', keywords='*')
omtbl['preprocess'] = fdzimlist
omtbl['now'] = fdzimlist
del fdzimlist
for key in omtbl.keys(): omtbl[key] = omtbl[key].astype(strtype)
#------------------------------------------------------------
#%%
#	Defringe 
#------------------------------------------------------------
print(f"""{'-'*60}\n#\tFRINGE CORRECTION\n{'-'*60}""")
# for filte in set(frgtbl[(frgtbl['ccd'] == ccdtype) & (frgtbl['obs'] == obs)]['filter']):
for filte in set(fdzic.summary['filter']):
	frgtbl_ = frgtbl[
		(frgtbl['filter']==filte) &
		(frgtbl['ccd'] == ccdtype) &
		(frgtbl['obs'] == obs)
		]
	if len(frgtbl_) > 0:
		if __name__ == '__main__':
			with multiprocessing.Pool(processes=ncore) as pool:
				results = pool.starmap(
					defringe,
					zip(
						fdzic.filter(filter=filte).files,
						repeat(frgtbl_['image'][0]),
						add_prefix(fdzic.filter(filter=filte).files, 'df'),
						repeat(frgtbl_['table'][0]),
						repeat(10)
						)
					)
				for i, inim in enumerate(fdzic.filter(filter=filte).files):
					indx_tmp = np.where(omtbl['now'] == inim)
					omtbl['now'][indx_tmp] = results[i]
					omtbl['defringe'][indx_tmp] = results[i]
				del indx_tmp
				for key in omtbl.keys(): omtbl[key] = omtbl[key].astype(strtype)
	else:
		print(f'{filte} : N/A')
#------------------------------------------------------------
#%%
#	FIX (TMP)
#------------------------------------------------------------


#------------------------------------------------------------
#	COSMIC-RAY REMOVAL
#------------------------------------------------------------
#	Seeing measurement w/ simple SE
prefix = 'simple'
path_conf = f'{path_config}/{prefix}.sex'
path_param = f'{path_config}/{prefix}.param'
path_nnw = f'{path_config}/{prefix}.nnw'
path_conv = f'{path_config}/{prefix}.conv'
#	Single
'''
inim = omtbl['now'][0]
seeing, peeing = get_seeing(
	inim,
	gain, 
	pixscale, 
	fov, 
	path_conf, 
	path_param, 
	path_conv, 
	path_nnw, 
	seeing_assume=3*u.arcsec, 
	frac=0.68, 
	n_min_src=5
	)
'''
if __name__ == '__main__':
	#	Seeing measurement
	print(f"""{'-'*60}\n#\tSEEING MEASUREMENT\n{'-'*60}""")
	with multiprocessing.Pool(processes=ncore) as pool:
		results = pool.starmap(
			get_seeing,
			zip(
					omtbl['now'],
					repeat(gain), 
					repeat(pixscale), 
					repeat(fov),
					repeat(path_conf), 
					repeat(path_param), 
					repeat(path_conv), 
					repeat(path_nnw), 
					repeat(3*u.arcsec),
					repeat(0.68),
					repeat(5),
			)
		)
		print('DONE')
		#	Cosmic-ray removal
		print(f"""{'-'*60}\n#\tCOSMIC-RAY REMOVAL\n{'-'*60}""")
		cleantype = 'medmask'
		if __name__ == '__main__':
			with multiprocessing.Pool(processes=ncore) as pool:
				results = pool.starmap(
					cosmic_ray_removal,
					zip(
							omtbl['now'],
							add_prefix(omtbl['now'], 'cr'),
							add_suffix(omtbl['final'], 'mask'),
							repeat(gain), 
							repeat(rdnoise), 
							[r[0] for r in results], 
							repeat(cleantype)
						)
					)

for i, inim in enumerate(omtbl['now']):
	tmpim = add_prefix(omtbl['now'], 'cr')[i]
	indx_tmp = np.where(omtbl['now'] == inim)
	omtbl['now'][indx_tmp] = tmpim
	omtbl['cosmic_ray_removal'][indx_tmp] = tmpim
del tmpim
del indx_tmp
for key in omtbl.keys(): omtbl[key] = omtbl[key].astype(strtype)
#------------------------------------------------------------
#%%
#	ASTROMETRY
#------------------------------------------------------------
print(f"""{'-'*60}\n#\tASTROMETRY\n{'-'*60}""")
#	Classification IMSNG and non-IMSNG target
frac = 0.10 #	Pixel scale inverval ratio
tralist, tdeclist = [], []
for i, inim in enumerate(omtbl['now']):
	if fits.getheader(inim)['OBJECT'] in alltbl['obj']:
		indx_obj = np.where(fits.getheader(inim)['OBJECT']==alltbl['obj'])
		tra, tdec = alltbl['ra'][indx_obj].item(), alltbl['dec'][indx_obj].item()
	else:
		tra, tdec = None, None
	tralist.append(tra)
	tdeclist.append(tdec)
#	Multi-processing
if __name__ == '__main__':
	with multiprocessing.Pool(processes=ncore) as pool:
		results = pool.starmap(
			astrometry,
			zip(
					omtbl['now'],
					add_prefix(omtbl['now'], 'a'),
					repeat(pixscale), 
					repeat(frac),
					tralist,
					tdeclist,
					repeat(fov),
					repeat(30)
				)
			)
#	Check astrometry results
print(f"""{'-'*60}\n#\tCHECK ASTROMETRY RESULTS\n{'-'*60}""")
c_all = SkyCoord(alltbl['ra'], alltbl['dec'], unit=(u.hourangle, u.deg))
for i, inim in enumerate(add_prefix(omtbl['now'], 'a')):
	hdr = fits.getheader(inim)
	if os.path.exists(inim):
		print(f'[{i}] {os.path.basename(inim)} : Astrometry Success')
		c = SkyCoord(hdr['CRVAL1'], hdr['CRVAL2'], unit=u.deg)
		indx_tmp, sep, _ = c.match_to_catalog_sky(c_all)
		ra_offset, dec_offset = c.spherical_offsets_to(c_all)

		#	Put pointing offset info.
		fits.setval(inim, keyword='CNTSEP', value=round(sep[0].arcmin, 3), comment='Offset between pointing and galaxy position [arcmin]')
		fits.setval(inim, keyword='CNTRAOFF', value=round(ra_offset.arcmin[indx_tmp], 3), comment='RA offset between pointing and galaxy position [arcmin]')
		fits.setval(inim, keyword='CNTDEOFF', value=round(dec_offset.arcmin[indx_tmp], 3), comment='Dec offset between pointing and galaxy position [arcmin]')
		
		print('\tCalculate accuracy')
		astrometry_analysis(
			inim=omtbl['now'][i], 
			incor=f"{os.path.splitext(omtbl['now'][i])[0]}.corr",
			# outpng=f'{os.path.splitext(inim)[0]}.astrm.png',
			# outdat=f'{os.path.splitext(inim)[0]}.astrm.dat'
			outpng=f"{os.path.splitext(omtbl['final'][i])[0]}.astrm.png",
			outdat=f"{os.path.splitext(omtbl['final'][i])[0]}.astrm.dat"
			)
		#	Update
		# omtbl['now'][i] = inim
	else:		
		print(f'{i} {os.path.basename(inim)} : Astrometry Fail')
		#	Suspect of wrong object name
		if omtbl['objtype'][i] == 'IMSNG':
			print('\tThis is IMSNG target. Start re-astronomy.')
			#	Retry astrometry
			astrometry(
				inim=omtbl['now'][i], 
				outim=add_prefix(omtbl['now'], 'a')[i], 
				pixscale=pixscale, 
				frac=frac, 
				cpulimit=60
				)
			if os.path.exists(inim):
				print('\tRe-astrometry SUCCESS!')
				c = SkyCoord(hdr['CRVAL1'], hdr['CRVAL2'], unit=u.deg)
				indx_tmp, sep, _ = c.match_to_catalog_sky(c_all)
				if sep[0] < fov:
					newobj = alltbl['obj'][indx_tmp]
					print(f"\tCorrect OBJECT HEADER, {hdr['OBJECT']} --> {newobj} position.")
					fits.setval(inim, keyword='OBJECT', value=newobj)
					#	Put pointing offset info.
					fits.setval(inim, keyword='CENTSEP', value=round(sep[0].arcmin, 3), comment='Offset between pointing and galaxy position')
					ra_offset, dec_offset = c.spherical_offsets_to(c_all)
					fits.setval(inim, keyword='CNTSEP', value=round(sep[0].arcmin, 3), comment='Offset between pointing and galaxy position [arcmin]')
					fits.setval(inim, keyword='CNTRAOFF', value=round(ra_offset.arcmin, 3), comment='RA offset between pointing and galaxy position [arcmin]')
					fits.setval(inim, keyword='CNTDEOFF', value=round(dec_offset.arcmin, 3), comment='Dec offset between pointing and galaxy position [arcmin]')
		
				astrometry_analysis(
					inim=omtbl['now'][i], 
					incor=f"{os.path.splitext(omtbl['now'][i])[0]}.corr",
					outpng=f"{os.path.splitext(omtbl['final'][i])[0]}.astrm.png",
					outdat=f"{os.path.splitext(omtbl['final'][i])[0]}.astrm.dat"
					)
				# omtbl['now'][i] = inim
				pass
			else:
				print('\tRe-astrometry Fail...')
				pass
		else:
			print('\tNo other actions')
	del hdr
#	
for i, inim in enumerate(omtbl['now']):
	tmpim = add_prefix(omtbl['now'], 'a')[i]
	if os.path.exists(tmpim):
		indx_tmp = np.where(omtbl['now'] == inim)
		omtbl['now'][indx_tmp] = tmpim
		omtbl['astrometry'][indx_tmp] = tmpim
		del indx_tmp
del tmpim
for key in omtbl.keys(): omtbl[key] = omtbl[key].astype(strtype)
#------------------------------------------------------------
#%%
#	File name change
#------------------------------------------------------------
print(f"""{'-'*60}\n#\tFILENAME CHANGE to IMSNG FORMAT\n{'-'*60}""")
for i, inim in enumerate(omtbl['now']):
	newim = f"{omtbl['final'][i]}"
	com = f'cp {inim} {newim}'
	os.system(com)
	print(f'[{i}] {os.path.basename(inim)} --> {os.path.basename(newim)}')
	fits.setval(newim, keyword='IMNAME', value=os.path.basename(newim), comment='Formatted file name by gpPy process')

ic_cal = ImageFileCollection(path_data, glob_include='Calib-*.f*', glob_exclude='*mask*')
#------------------------------------------------------------
#%%
#	Photometry
#------------------------------------------------------------
photcom = f"python {path_phot} '{path_data}/Calib-*0.fits' {ncore}"
print(photcom)
os.system(photcom)
#------------------------------------------------------------
#%%
#	IMAGE COMBINE
#------------------------------------------------------------
def group_images(objtbl, tfrac):
	delt = np.array(objtbl['jd'] - np.min(objtbl['jd']))*(24*60*60) # [days] --> [sec]
	tsub = delt - (np.cumsum(objtbl['exptime']*tfrac) - objtbl['exptime'][0])
	indx = np.where(tsub < 0)
	indx_inv = np.where(tsub > 0)
	return indx, indx_inv
#------------------------------------------------------------
print(f"""{'-'*60}\n#\tGROUP IMAGES\n{'-'*60}""")
tfrac = 1.5 # Time fraction for grouping
comimlist = []
for obj in set(ic_cal.summary['object']):
	for filte in set(ic_cal.filter(object=obj).summary['filter']):
		print(f"{len(ic_cal.filter(object=obj).summary['filter'])} images for {obj} in {filte}")
		print('-'*60)

		ic_obj = ic_cal.filter(object=obj, filter=filte)
		objtbl = ic_obj.summary

		indx_com, indx_inv = group_images(objtbl, tfrac)
		comimlist.append(objtbl[indx_com])

		#	Numbering
		n=0
		for inim in objtbl['file'][indx_com]:
			print(f'[{n}] {os.path.basename(inim)}')
			n+=1
		print('-'*60)

		while len(indx_inv[0]):
			objtbl = objtbl[indx_inv]
			indx_com, indx_inv = group_images(objtbl, tfrac)
			comimlist.append(objtbl[indx_com])
			for inim in objtbl['file'][indx_com]:
				print(f'[{n}] {os.path.basename(inim)}')
				n+=1
			print('-'*60)
#------------------------------------------------------------
# %%
print(f"""{'-'*60}\n#\tALIGN IMAGES AND COMBINE\n{'-'*60}""")
for i in range(len(comimlist)):
	#	Alignment
	##	Target image
	imtbl = comimlist[i]
	indx_ref = np.where(imtbl['seeing']==np.max(imtbl['seeing']))
	tgtim = imtbl['file'][imtbl['file']==imtbl['file'][indx_ref]][0]
	##	Source image
	srcimlist = list(imtbl['file'][imtbl['file']!=tgtim])
	#	Combine
	aligned_imlist = [CCDData(fits.getdata(tgtim), unit='adu', header=fits.getheader(tgtim))]
	for srcim in srcimlist: aligned_imlist.append(align_astroalign(srcim, tgtim, zero=False))
	imcombine_ccddata(aligned_imlist, fluxscale=True, zpkey='ZP', nref=0, imlist=None,)
#------------------------------------------------------------
# %%
#	Photometry
#------------------------------------------------------------
photcom = f"python {path_phot} '{path_data}/Calib-*com.fits' {ncore}"
print(photcom)
os.system(photcom)
#%%
'''
#------------------------------------------------------------
#	IMAGE COMBINE
#------------------------------------------------------------
print(f"""{'='*60}\n#\tIMAGE COMBINE\n{'='*60}""")

def group_images(objtbl, tfrac):
	delt = np.array(objtbl['jd'] - np.min(objtbl['jd']))*(24*60*60) # [days] --> [sec]
	tsub = delt - (np.cumsum(objtbl['exptime']*tfrac) - objtbl['exptime'][0])
	indx = np.where(tsub < 0)
	indx_inv = np.where(tsub > 0)
	return indx, indx_inv

print(f"""{'-'*60}\n#\tGROUP IMAGES\n{'-'*60}""")
tfrac = 1.5 # Time fraction for grouping
comimlist = []
for obj in set(ic_cal.summary['object']):
	for filte in set(ic_cal.filter(object=obj).summary['filter']):
		print(f"{len(ic_cal.filter(object=obj).summary['filter'])} images for {obj} in {filte}")
		print('-'*60)

		ic_obj = ic_cal.filter(object=obj, filter=filte)
		objtbl = ic_obj.summary

		indx_com, indx_inv = group_images(objtbl, tfrac)
		comimlist.append(objtbl[indx_com])

		#	Numbering
		n=0
		for inim in objtbl['file'][indx_com]:
			print(f'[{n}] {os.path.basename(inim)}')
			n+=1
		print('-'*60)

		while len(indx_inv[0]):
			objtbl = objtbl[indx_inv]
			indx_com, indx_inv = group_images(objtbl, tfrac)
			comimlist.append(objtbl[indx_com])
			for inim in objtbl['file'][indx_com]:
				print(f'[{n}] {os.path.basename(inim)}')
				n+=1
			print('-'*60)

#------------------------------------------------------------
# %%
print(f"""{'-'*60}\n#\tALIGN IMAGES AND COMBINE\n{'-'*60}""")
for i in range(len(comimlist)):
	#	Alignment
	##	Target image
	imtbl = comimlist[i]
	indx_ref = np.where(imtbl['seeing']==np.min(imtbl['seeing']))
	tgtim = imtbl['file'][imtbl['file']==imtbl['file'][indx_ref]][0]
	##	Source image
	srcimlist = list(imtbl['file'][imtbl['file']!=tgtim])
	#	Combine
	aligned_imlist = [CCDData(fits.getdata(tgtim), unit='adu', header=fits.getheader(tgtim))]
	for srcim in srcimlist: aligned_imlist.append(align_astroalign(srcim, tgtim,))
	imcombine_ccddata(aligned_imlist, imlist=None)
#------------------------------------------------------------
#	Seeing measurement again
if __name__ == '__main__':
	print(f"""{'-'*60}\n#\tSEEING MEASUREMENT\n{'-'*60}""")
	with multiprocessing.Pool(processes=ncore) as pool:
		results = pool.starmap(
			get_seeing,
			zip(
					omtbl['now'],
					repeat(gain), 
					repeat(pixscale), 
					repeat(fov),
					repeat(path_conf), 
					repeat(path_param), 
					repeat(path_conv), 
					repeat(path_nnw), 
					repeat(3*u.arcsec),
					repeat(0.68),
					repeat(5),
			)
		)
		print('DONE')'''