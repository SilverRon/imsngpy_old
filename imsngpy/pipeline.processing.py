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
import sys
sys.path.append('/home/paek/imsngpy')
#	IMSNGpy modules
from tableutil import getccdinfo
from preprocess import *
from misc import *
from phot import *
import warnings
warnings.filterwarnings(action='ignore')
import time
start_localtime = time.strftime('%Y-%m-%d %H:%M:%S (%Z)', time.localtime())
from astropy.io import ascii
from astropy.io import fits
import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from ccdproc import ImageFileCollection
from astropy.time import Time
#	Multiprocess tools
from itertools import repeat
import multiprocessing
# from __future__ import print_function, division, absolute_import
# from timeit import default_timer as timer
# from numba import jit
# from pyraf import iraf
# import matplotlib.pyplot as plt
# plt.ioff()
# from astropy.nddata import CCDData
# from imsng import calib
# from imsng import tool_tbd
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
path_raw = '/data6/obsdata/LOAO/1969_0119'
obs = 'LOAO'
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
path_table = '/home/paek/table'
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
ccdtbl = ascii.read(f'{path_table}/ccd.dat') 

path_data = f'{path_obs}/{os.path.basename(path_raw)}'

ccdinfo = getccdinfo(obs, ccddat)
gain = ccdinfo['gain']
rdnoise = ccdinfo['rdnoise']
pixscale = ccdinfo['pixelscale']
fov = ccdinfo['fov']


if os.path.exists(path_data):
	rmcom = f'rm -rf {path_data}'
	print(rmcom)
	os.system(rmcom)

print(f"""{'-'*60}\n#\tCOPY DATA\n{'-'*60}""")
cpcom = f'cp -r {path_raw} {path_data}'
print(cpcom)
os.system(cpcom)

ic0 = ImageFileCollection(path_data, keywords='*')
# ic0.summary.write('{}/hdr.raw.dat'.format(path_data), format='ascii.tab', overwrite=True)
#------------------------------------------------------------
#%%
#	Header correction
#------------------------------------------------------------
comment = f"""{'-'*60}\n#\tHEADER CORRECTION\n{'-'*60}"""
print(comment)

for i, inim in enumerate(ic0.summary['file']):
	#	CCD Type
	if i == 0:
		for key in set(ccdtbl['key']):
			if key.lower() in ic0.summary.keys():
				ccdtype = ccdtbl[
					(ccdtbl['value'] == ic0.summary[key.lower()][i]) &
					(ccdtbl['obs'] == obs)
				]['ccd'].item()
			else:
				ccdtype = 'UNKNOWN'
	fits.setval(f'{path_data}/{inim}', 'CCDNAME', value=ccdtype)
	fits.setval(f'{path_data}/{inim}', 'OBSERVAT', value=obs)
	#	Correction with table
	for key, val, nval in zip(hdrtbl['key'], hdrtbl['val'], hdrtbl['newval']):
		if ic0.summary[key.lower()][i] == val:
			print(f'{inim} - {key} : {val} --> {nval}')
			fits.setval(f'{path_data}/{inim}', key, value=nval)
	#	DATE-OBS, JD, MJD
	if 'T' not in ic0.summary['date-obs'][i]:
		dateobs = f"{ic0.summary['date-obs']}'T'{ic0.summary['time-obs']}"
		fits.setval(f'{path_data}/{inim}', 'date-obs', value=dateobs)
		del dateobs
	else:
		pass
	t = Time(ic0.summary['date-obs'][i], format='isot')
	fits.setval(f'{path_data}/{inim}', 'JD', value=t.jd)
	fits.setval(f'{path_data}/{inim}', 'MJD', value=t.mjd)
	del t
	#	OBJECT name
	if 'ngc' in ic0.summary['object'][i]:
		objectname = ic0.summary['object'][i]
		while len(objectname)<7:
			head = objectname[0:3]
			tail = objectname[3:]
			tail = f'0{tail}'
			objectname = f'{head}{tail}'
		fits.setval(f'{path_data}/{inim}', 'OBJECT', value=objectname.upper())
		del objectname
		del head
		del tail
	print()
ic1 = ImageFileCollection(path_data, keywords='*')
#%%
#------------------------------------------------------------
#	Object Master Table
#
#	Write the status of processing
#	Pointing the original filename and updated filename
#	Each dtype is set to 'strtype' variable
strtype = 'U300'
omtbl = Table()
omtbl['raw'] = [inim for inim in ic1.filter(imagetyp='object').files]
omtbl['now'] = omtbl['raw']
omtbl['preprocess'] = ''
omtbl['defringe'] = ''
omtbl['cosmic_ray_removal'] = ''
omtbl['astrometry'] = ''
omtbl['final'] = [fnamechange(inim) for inim in ic1.filter(imagetyp='object').files]
for key in omtbl.keys(): omtbl[key] = omtbl[key].astype(strtype)
#============================================================
#%%
#	Pre-processing
#------------------------------------------------------------
mframe = dict()
#------------------------------------------------------------
#	BIAS
#------------------------------------------------------------
biaslist = ic1.filter(imagetyp='bias').files
if len(biaslist) > 0:
	mframe['bias'] = master_bias(biaslist)
else:
	print('borrow')
#------------------------------------------------------------
#	DARK
#------------------------------------------------------------
darkframe = dict()
darklist = ic1.filter(imagetyp='dark').files
if len(darklist) > 0:
	darkexptime = np.array(list(set(ic1.filter(imagetyp='dark').summary['exptime'])))
	for exptime in darkexptime:
		darkframe[f'{int(exptime)}'] = master_dark(ic1.filter(imagetyp='dark', exptime=exptime).files, mbias=mframe['bias'])
else:
	print('borrow')

mframe['dark'] = darkframe
del darkframe
#------------------------------------------------------------
#	FLAT
#------------------------------------------------------------
flatframe = dict()
flatlist = ic1.filter(imagetyp='flat').files
if len(flatlist) > 0:
	indx_mindark = np.where(darkexptime == np.min(darkexptime))
	mdark = mframe['dark'][str(int(darkexptime[indx_mindark].item()))]
	for filte in set(ic1.filter(imagetyp='flat',).summary['filter']):
		# print(filte)
		flatframe[filte] = master_flat(ic1.filter(imagetyp='flat', filter=filte).files, mbias=mframe['bias'], mdark=mdark, filte=filte)
	del mdark
else:
	print('borrow')
mframe['flat'] = flatframe
del flatframe
#------------------------------------------------------------
#	OBJECT Correction
#------------------------------------------------------------
print(f"""{'-'*60}\n#\tOBJECT CORRECTION ({len(ic1.filter(imagetyp='object').files)})\n{'-'*60}""")
#	Function for multi-process
def obj_process4mp(inim, newim, filte, exptime, darkexptime, ccdinfo, mframe,):
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
		readnoise=ccdinfo['rdnoise'],
		mbias=mframe['bias'],
		mdark=mframe['dark'][str(int(bestdarkexptime))],
		mflat=mframe['flat'][filte],
	)
	# nccd.write(f'{os.path.dirname(inim)}/fdz{os.path.basename(inim)}', overwrite=True)
	nccd.write(newim, overwrite=True)
#	Run with multi-process
fdzimlist = add_prepix(ic1.filter(imagetyp='object').files, 'fdz')
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
				repeat(ccdinfo),
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
						add_prepix(fdzic.filter(filter=filte).files, 'df'),
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
#%%
#------------------------------------------------------------
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
		#	Cosmic-ray removal
		print(f"""{'-'*60}\n#\tCOSMIC-RAY REMOVAL\n{'-'*60}""")
		cleantype = 'medmask'
		if __name__ == '__main__':
			with multiprocessing.Pool(processes=ncore) as pool:
				results = pool.starmap(
					cosmic_ray_removal,
					zip(
							omtbl['now'],
							add_prepix(omtbl['now'], 'cr'),
							repeat(gain), 
							repeat(rdnoise), 
							[r[0] for r in results], 
							repeat(cleantype)
						)
					)

for i, inim in enumerate(omtbl['now']):
	tmpim = add_prepix(omtbl['now'], 'cr')[i]
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
					add_prepix(omtbl['now'], 'a'),
					repeat(pixscale), 
					repeat(frac),
					tralist,
					tdeclist,
					repeat(fov),
					repeat(60)
				)
			)
#	Check astrometry results
c_all = SkyCoord(alltbl['ra'], alltbl['dec'], unit=(u.hourangle, u.deg))
for i, inim in enumerate(add_prepix(omtbl['now'], 'a')):
	if os.path.exists(inim):
		print(f'{i} {os.path.basename(inim)} : Astrometry Success')
		hdr = fits.getheader(inim)
		c = SkyCoord(hdr['CRVAL1'], hdr['CRVAL2'], unit=u.deg)
		indx_tmp, sep, _ = c.match_to_catalog_sky(c_all)
		if sep.item() < fov:
			print(f"Correct {hdr['OBJECT']} position.")
		del hdr
		
	else:
		print(f'{i} {os.path.basename(inim)} : Astrometry Fail')


astrometry_analysis(
	inim=outim,
	incor=f'{os.path.splitext(inim)[0]}.corr',
	outpng=f'{os.path.splitext(outim)[0]}.astrm.png',
	outdat=f'{os.path.splitext(outim)[0]}.astrm.dat'
	)
#%%
