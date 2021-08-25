#============================================================
#	IMSNG Pipeline
#	=> Processing
#	Data Monitoring => Processing => Transient Search
#============================================================
#	Library
#------------------------------------------------------------
from tableutil import getccdinfo
import warnings
warnings.filterwarnings(action='ignore')
import time
start_localtime = time.strftime('%Y-%m-%d %H:%M:%S (%Z)', time.localtime())
import os
import sys
import glob
from astropy.io import ascii
from astropy.io import fits
import numpy as np
from astropy.table import Table
from astropy import units as u
from ccdproc import ImageFileCollection
from astropy.time import Time
# from __future__ import print_function, division, absolute_import
# from timeit import default_timer as timer
# from numba import jit
# from pyraf import iraf
# import numpy as np
# import matplotlib.pyplot as plt
# plt.ioff()
# from astropy.nddata import CCDData
# from imsng import calib
# from imsng import tool_tbd
# from itertools import product
# from itertools import repeat
# import multiprocessing
#------------------------------------------------------------
#	My library
# from tableutil import *
#============================================================
#	USER SETTING
#============================================================
#	Input
#------------------------------------------------------------
#	Observatory_ccd
try:
	path_raw = sys.argv[1]
except:
	path_raw = input('''# Data folder to process : ''')

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
#	The number of cores
try:
	ncores = int(sys.argv[3])
except:
	ncores = 8

#	Test setting
path_raw = '/data6/obsdata/LOAO/1994_1026'
obs = 'LOAO'
ncores = 8

#------------------------------------------------------------
#	PATH
#------------------------------------------------------------
path_factory = '/data3/paek/factory'
path_gal = '/data6/IMSNG/IMSNGgalaxies'
path_config = '/home/paek/config'
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


path_data = f'{path_obs}/{os.path.basename(path_raw)}'





'''
for inim in imlist:
	for key, val in zip(keys, vals):
		fits.setval(inim, key, value=val)
'''














#============================================================
#	MAIN BODY
#============================================================
tdict = dict()
starttime = time.time()
#	Remove old folder and re-copy folder
rmcom = f'rm -rf {path_data}'
print(rmcom)
os.system(rmcom)
#	Copy raw data to factory
cpcom = f'cp -r {path_raw} {path_data}'
print(cpcom)
os.system(cpcom)
#------------------------------------------------------------
#	Process summary status
#------------------------------------------------------------
timetbl = Table()
timetbl['process'] = [
	'master_frame',
	'pre_process',
	'astrometry',
	'cr_removal',
	'defringe',
	'photometry',
	'image_stack',
	'photometry_com',
	'subtraction',
	'photometry_sub',
	'transient_search',
	'total'
	]
timetbl['status'] = False
timetbl['time'] = 0.0 * u.second
#------------------------------------------------------------
#	CCD TYPE
#------------------------------------------------------------
ic0 = ImageFileCollection(path_data, keywords='*')
ic0.summary.write('{}/hdr.raw.dat'.format(path_data), format='ascii.tab', overwrite=True) 

#	Exceptions
#	DOAO CCD
if obs == 'DOAO':
	instrume = ic0.summary['instrume'][0]
	if instrume == 'Apogee USB/Net':
		obs = 'DOAO_APOGEE'
	elif instrume == '':			
		obs = 'DOAO_FLI'
	elif instrume == 'FLI':			
		obs = 'DOAO_FLI'
else:
	pass
#	KHAO Binning
if obs == 'KHAO':
	xbinning = ic0.summary['xbinning'][0]
	if xbinning > 1:
		ccdinfo['pixscale'] = ccdinfo['pixscale']*2
#	RASA Mode
if obs == 'RASA36':
	if 'hdr' in path:
		mode = 'hdr'
		badmode = True
	elif 'high' in path:
		mode = 'high'
		badmode = True
	else:
		biasim = ic0.summary['file'][ic0.summary['imagetyp']=='Bias Frame'][0]
		biaslevel = np.median(fits.getdata(f'{path_data}/{biasim}').flatten())
		if biaslevel > 87:
			mode = 'hdr'
		else:
			mode = 'high'
		badmode = False
	# print(f'RASA36 [{mode} mode (Bad mode : {badmode})]')# : bias level = {biaslevel}')
	print(f'Master Frame Mode:{mode} [Bad Mode:{badmode}]')

# ccdinfo = tableutil.getccdinfo(obs, ccddat)
ccdinfo = getccdinfo(obs, ccddat)



#	Unifying Header
for file in ic0.summary['file']:
	inim = f'{path_data}/{file}'
	with fits.open(inim) as hdul:
		hdr = hdul[0].header
		#	Observatory
		fits.setval(inim, keyword='OBS', value=obs, comment='observation location')
		#	CHECK DATE-OBS
		if 'T' not in hdr['DATE-OBS']:
			dateobs = hdr['DATE-OBS']+'T'+hdr['TIME-OBS']
			fits.setval(inim, keyword='DATE-OBS', value=dateobs,)
		else:
			dateobs = hdr['DATE-OBS']
		t = Time(hdr['DATE-OBS'], format='isot')
		fits.setval(inim, keyword='JD', value=t.jd, comment='Julian Date')
		fits.setval(inim, keyword='MJD', value=t.mjd, comment='Modified Julian Date')
		#	CHECK OBJECT HDR
		obj = hdr['OBJECT']
		if 'ngc' in obj:
			while len(obj)<7:
				head = obj[0:3]
				tail = obj[3: ]
				tail = '0'+tail
				obj = head+tail
		obj = obj.upper()
		fits.setval(inim, keyword='object', value=obj,)
		#	Other modifications
		for key, val, newval in zip(hdrtbl['key'], hdrtbl['val'], hdrtbl['newval']):
			if hdr[key] == val:
				fits.setval(inim, key, value=newval)
				print(inim, key, val, newval)








ic1 = ImageFileCollection(path_data, keywords='*')
ic1.summary.write(f'{path_data}/hdr.cor.dat', format='ascii.tab', overwrite=True)

nobj = len(ic1.filter(imagetyp='OBJECT').summary)
#------------------------------------------------------------
#	Slack messages
#------------------------------------------------------------	
path_keys = '/home/paek/table'
keytbl = ascii.read(f'{path_keys}/keys.dat')
OAuth_Token = keytbl['key'][keytbl['name']=='slack'].item()

channel = '#pipeline'
text = f'[Pipeline/{obs}] Start Processing {os.path.basename(path_data)} Data ({nobj} objects) with {ncores} cores'

param_slack = dict(
	token = OAuth_Token,
	channel = channel,
	text = text,
)

# misc.slack_bot(**param_slack)
slack_bot(**param_slack)


#============================================================
#	Making master frames (BIAS, DARK, FLAT)
#============================================================
st = time.time()
#------------------------------------------------------------
#	BIAS
#------------------------------------------------------------
biastbl = ic1.filter(imagetyp='Bias').summary
biaslist = list(biastbl['file'])
# mbias = preprocess.master_bias(imlist)
mbias = master_bias(imlist)

darktbl = ic1.filter(imagetyp='dark').summary




try:
	biasnumb = len(ic1.filter(imagetyp='Bias').summary)
except:
	biasnumb = 0
# if len(ic1.filter(imagetyp='Bias').summary) != 0:
if biasnumb != 0:
	mzero = calib.master_zero(ic1, fig=False)
	# print(zeroim)
	date = fits.getheader(f'{path_data}/zero.fits')['date-obs'][:10].replace('-', '')
	if obs == 'RASA36':
		zeroim = f'{path_mframe}/{obs}/zero/{date}-zero_{mode}.fits'
	else:
		zeroim = f'{path_mframe}/{obs}/zero/{date}-zero.fits'
	cpcom = f'cp {path_data}/zero.fits {zeroim}'
	print(cpcom)
	os.system(cpcom)
	plt.close('all')
else:
	#	IF THERE IS NO FLAT FRAMES, BORROW FROM CLOSEST OTHER DATE
	print('\nNO BIAS FRAMES\n')
	if obs == 'RASA36':
		pastzero = np.array(glob.glob(f'{path_mframe}/{obs}/zero/*zero_{mode}.fits'))
	else:
		pastzero = np.array(glob.glob(f'{path_mframe}/{obs}/zero/*zero.fits'))
	#	CALCULATE CLOSEST ONE FROM TIME DIFFERENCE
	deltime = []
	for date in pastzero:
		zeromjd = calib.isot_to_mjd((os.path.basename(date)).split('-')[0])
		deltime.append(np.abs(ic1.summary['mjd'][0]-zeromjd))
	indx_closet = np.where(deltime == np.min(deltime))
	tmpzero = path_data+'/'+os.path.basename(np.asscalar(pastzero[indx_closet]))
	cpcom = 'cp {} {}'.format(np.asscalar(pastzero[indx_closet]), tmpzero)
	print(cpcom)
	os.system(cpcom)
	# if obs != 'KCT':
	if 'KCT_ASI1600MM' in obs:
		#KCT Exception
		mzero = CCDData.read(tmpzero, hdu=0, unit='adu')
	elif obs == 'RASA36':
		if (mode == 'high') & (badmode == True):
			mzero = CCDData.read(tmpzero, hdu=0).multiply(20)
			print('[Bad mode] Multiply 20 on the high mode bias.')
		else:
			mzero = CCDData.read(tmpzero, hdu=0)
	else:
		mzero = CCDData.read(tmpzero, hdu=0)#, unit='adu')
	mzero.meta['FILENAME'] = os.path.basename(tmpzero)
'''
#	RASA36 discreminator for high/hdr modes with bias level
if obs == 'RASA36':
	biaslevel = np.median(mzero)
	if biaslevel > 87:
		mode = 'hdr'
	else:
		mode = 'high'
	print(f'RASA36 [{mode} mode] : bias level = {biaslevel}')
'''
#------------------------------------------------------------
#	DARK (ITERATION FOR EACH EXPOSURE TIMES)
#------------------------------------------------------------
try:
	darkexptimelist = sorted(list(set(ic1.filter(imagetyp='dark').summary['exptime'])))
	darknumb = len(darkexptimelist)
except:
	darknumb = 0
darkdict = dict()
# if len(darkexptimelist) != 0:
if darknumb != 0:
	dark_process = True
	for i, exptime in enumerate(darkexptimelist):
		print('PRE PROCESS FOR DARK ({} sec)\t[{}/{}]'.format(exptime, i+1, len(darkexptimelist)))
		mdark = calib.master_dark(ic1, mzero=mzero, exptime=exptime, fig=False)
		darkdict['{}'.format(int(exptime))] = mdark

		date = fits.getheader(f'{path_data}/dark-{int(exptime)}.fits')['date-obs'][:10].replace('-', '')
		if obs == 'RASA36':
			darkim = f'{path_mframe}/{obs}/dark/{int(exptime)}-{date}-dark_{mode}.fits'
		else:
			darkim = f'{path_mframe}/{obs}/dark/{int(exptime)}-{date}-dark.fits'
		# print(zeroim)
		cpcom = 'cp {}/dark-{}.fits {}'.format(path_data, int(exptime), darkim)
		print(cpcom)
		os.system(cpcom)
		plt.close('all')
else:
	#	Borrow
	print('\nNO DARK FRAMES\n')
	objexptimelist = sorted(list(set(ic1.filter(imagetyp='object').summary['exptime'])))
	exptime = objexptimelist[-1]
	# pastdark = np.array(glob.glob('{}/{}/dark/{}*dark*.fits'.format(path_mframe, obs, int(exptime))))
	if obs == 'RASA36':
		pastdark = np.array(glob.glob(f'{path_mframe}/{obs}/dark/{int(exptime)}*dark_{mode}.fits'))
	else:
		pastdark = np.array(glob.glob(f'{path_mframe}/{obs}/dark/{int(exptime)}*dark.fits'))

	if len(pastdark) == 0:
		pastdark = np.array(glob.glob('{}/{}/dark/*dark*.fits'.format(path_mframe, obs)))
	else:
		pass
	#	CALCULATE CLOSEST ONE FROM TIME DIFFERENCE
	deltime = []
	delexptime = []
	darkexptimes = []
	for date in pastdark:
		# darkmjd = calib.isot_to_mjd((os.path.basename(date)).split('-')[0])
		darkmjd = calib.isot_to_mjd((os.path.basename(date)).split('-')[1])
		darkexptime = int( os.path.basename(date).split('-')[0] )
		# darkexptime = delexptime.append(int( os.path.basename(date).split('-')[1] ))
		darkexptimes.append(darkexptime)
		deltime.append(np.abs(ic1.summary['mjd'][0]-darkmjd))
	if 'KCT' in obs:
		indx_closet = np.where(
			(np.abs(np.array(darkexptimes)-exptime) == np.min(np.abs(np.array(darkexptimes)-exptime)))
		)
	else:
		indx_closet = np.where(
			(deltime == np.min(deltime)) &
			(darkexptimes == np.max(darkexptimes))
		)
	if len(indx_closet[0]) == 0:
		indx_closet = np.where(
			(deltime == np.min(deltime))
		)
	else:
		pass
	# tmpdark = path_data+'/'+os.path.basename(pastdark[indx_closet].item())
	# tmpdark = '{}/{}'.format(path_data, os.path.basename(pastdark[indx_closet[0]].item()))
	# tmpdark = pastdark[indx_closet[0]].item()
	tmpdark = pastdark[indx_closet][-1]
	exptime = int(fits.getheader(tmpdark)['exptime'])

	# cpcom = 'cp {} {}/dark-{}.fits'.format(tmpdark, path_data, int(exptime))
	cpcom = 'cp {} {}'.format(tmpdark, path_data, int(exptime))
	print(cpcom)
	os.system(cpcom)
	if 'KCT' in obs:
		#KCT Exception
		mdark = CCDData.read(tmpdark, hdu=0, unit='adu')
	elif obs == 'RASA36':
		if (mode == 'high') & (badmode == True):
			mdark = CCDData.read(tmpdark, hdu=0).multiply(20)
			print('[Bad mode] Multiply 20 on the high mode dark.')
		else:
			mdark = CCDData.read(tmpdark, hdu=0)
	else:
		mdark = CCDData.read(tmpdark, hdu=0)#, unit='adu')
	mdark.meta['FILENAME'] = os.path.basename(tmpdark)
	mdark.meta['EXPTIME'] = exptime
	darkdict['{}'.format(int(exptime))] = mdark
#------------------------------------------------------------
#	FLAT (ITERATION FOR EACH FILTERS)
#------------------------------------------------------------
flatdict = dict()
try:
	flatfilterlist = list(set(ic1.filter(imagetyp='flat').summary['filter']))
	for i, filte in enumerate(flatfilterlist):
		# print(i, filte)
		print('MAKING MASTER FLAT IN {}-BAND'.format(filte))
		mflat = calib.master_flat(ic1, mzero, filte, mdark=mdark, fig=True)
		flatdict[filte] = mflat

		date = fits.getheader(f'{path_data}/dark-{int(exptime)}.fits')['date-obs'][:10].replace('-', '')
		if obs == 'RASA36':
			flatim = f'{path_mframe}/{obs}/flat/{date}-n{file}_{mode}.fits'
		else:
			flatim = f'{path_mframe}/{obs}/flat/{date}-n{file}.fits'
		cpcom = f'cp {path_data}/n{filte}.fits {flatim}'
		print(cpcom)
		os.system(cpcom)
		plt.close('all')
except:
	print('No flat calibration image.')
	# flatdict['None'] = None
	pass
# tdict['masterframe'] = time.time() - st
timetbl['status'][timetbl['process']=='master_frame'] = True
timetbl['time'][timetbl['process']=='master_frame'] = int(time.time() - st)
#------------------------------------------------------------
#	OBJECT CALIBARTION (ZERO, DARK, FLAT)
#------------------------------------------------------------
st_ = time.time()
comment     = '='*60+'\n' \
			+ 'OBJECT CALIBRATION\n' \
			+ '='*60+'\n'
print(comment)
objfilterlist = sorted(list(set(ic1.filter(imagetyp='object').summary['filter'])))
objexptimelist = sorted(list(set(ic1.filter(imagetyp='object').summary['exptime'])))
for i, filte in enumerate(objfilterlist):
	print('PRE PROCESS FOR {} FILTER OBJECT\t[{}/{}]'.format(filte, i+1, len(objfilterlist)))
	if filte in flatdict.keys():
		mflat = flatdict[filte]
	else:
		print('\nNO {} FLAT FRAMES\n'.format(filte))
		#	CALCULATE CLOSEST ONE FROM TIME DIFFERENCE
		deltime = []
		if obs != 'RASA36':
			pastflat = np.array(glob.glob('{}/{}/flat/*n{}*.fits'.format(path_mframe, obs, filte)))
			for date in pastflat:
				flatmjd = calib.isot_to_mjd((os.path.basename(date)).split('-')[0])
				deltime.append(np.abs(ic1.summary['mjd'][0]-flatmjd))
		elif obs == 'RASA36':
			pastflat = np.array(glob.glob(f'{path_mframe}/{obs}/flat/*_{mode}-n{filte}*.fits'))
			for date in pastflat:
				flatmjd = fits.getheader(date)['MJD']
				deltime.append(np.abs(ic1.summary['mjd'][0]-flatmjd))

		indx_closet = np.where(deltime == np.min(deltime))
		tmpflat = '{}/{}'.format(path_data, os.path.basename(pastflat[indx_closet][0].item()))
		# tmpflat = pastflat[indx_closet][0].item()
		cpcom = 'cp {} {}'.format(pastflat[indx_closet][0].item(), tmpflat)
		print(cpcom)
		os.system(cpcom)
		if 'KCT' not in obs:
			mflat = CCDData.read(tmpflat, hdu=0)#, unit='adu')
		else:
			#KCT Exception
			mflat = CCDData.read(tmpflat, hdu=0, unit='adu')
			mflat.meta['FILENAME'] = os.path.basename(tmpflat)
		flatdict[filte] = mflat
	for expt in objexptimelist:
		if str(int(expt)) in darkdict.keys():
			mdark = darkdict[str(int(expt))]
		else:
			mdark = darkdict[list(darkdict.keys())[-1]]
		calib.calibration(ic1, mzero, mflat, filte, mdark=mdark)
# tdict['objectcorrection'] = time.time() - st - tdict[list(tdict.keys())[-1]]
timetbl['status'][timetbl['process']=='pre_process'] = True
timetbl['time'][timetbl['process']=='pre_process'] = int(time.time() - st_)
#------------------------------------------------------------
#	ASTROMETRY
#------------------------------------------------------------
st_ = time.time()
print('ASTROMETRY START')
print('='*60)
astrometryfailist = []
fzimlist = []
for ims in ('{}/fz*.fits'.format(path_data), '{}/fz*.fit'.format(path_data), '{}/fz*.fts'.format(path_data)):
	fzimlist.extend(sorted(glob.glob(ims)))
# fzimlist = sorted(glob.glob(path_data+'/fz*.fits'))


astimlist = []
astotlist = []
astralist = []
astdelist = []

for inim in fzimlist:
	obj = (fits.getheader(inim)['object']).upper()
	if (obj in alltbl['obj']):
		indx_target = np.where(obj == alltbl['obj'])[0][0]
		ra, dec = alltbl['ra'][indx_target].item(), alltbl['dec'][indx_target].item()
		astimlist.append(inim)
		astralist.append(ra)
		astdelist.append(dec)
	else:
		astotlist.append(inim)
#	Astrometry (IMSNG field)
if __name__ == '__main__':
	with multiprocessing.Pool(processes=ncores) as pool:
		results = pool.starmap(calib.astrometry, zip(astimlist, repeat(ccdinfo['pixscale']), astralist, astdelist, repeat(ccdinfo['fov']/60.), repeat(15)))
#	Astrometry (non IMSNG field)
if __name__ == '__main__':
	with multiprocessing.Pool(processes=ncores) as pool:
		results = pool.starmap(calib.astrometry, zip(astotlist, repeat(ccdinfo['pixscale']), repeat(None), repeat(None), repeat(None), repeat(60)))
#	Astrometry (failed IMSNG field)
astfailist = []
for inim in astimlist:
	if (os.path.exists('{}/a{}'.format(path_data, os.path.basename(inim))) == False):
		astfailist.append(inim)
if __name__ == '__main__':
	with multiprocessing.Pool(processes=ncores) as pool:
		results = pool.starmap(calib.astrometry, zip(astfailist, repeat(ccdinfo['pixscale']), repeat(None), repeat(None), repeat(None), repeat(60)))
for inim in astfailist:
	if (os.path.exists('{}/a{}'.format(path_data, os.path.basename(inim))) == False):
		astrometryfailist.append('{}/a{}'.format(path_data, os.path.basename(inim)))

os.system('rm '+path_data+'/*.axy '+path_data+'/*.corr '+path_data+'/*.xyls '+path_data+'/*.match '+path_data+'/*.rdls '+path_data+'/*.solved '+path_data+'/*.wcs ')
print('ASTROMETRY COMPLETE\n'+'='*60)
# tdict['astronomy'] = time.time() - st - tdict[list(tdict.keys())[-1]]
timetbl['status'][timetbl['process']=='astrometry'] = True
timetbl['time'][timetbl['process']=='astrometry'] = int(time.time() - st_)
#------------------------------------------------------------
#	Quick seeing measurement with SE & Cosmic ray removal
#------------------------------------------------------------
st_ = time.time()
print('Quick seeing measurement with SE & Cosmic ray removal')
print('='*60)
gain = ccdinfo['gain'].value
rdnoise = ccdinfo['rdnoise']

# afzimlist = sorted(glob.glob(path_data+'/afz*.fits'))
afzimlist = []
for ims in ('{}/afz*.fits'.format(path_data), '{}/afz*.fit'.format(path_data), '{}/afz*.fts'.format(path_data)):
	afzimlist.extend(sorted(glob.glob(ims)))
outimlist = []
for i, inim in enumerate(afzimlist):
	outim = '{}/cr{}'.format(os.path.dirname(inim), os.path.basename(inim))
	outimlist.append(outim)
if ('KCT' not in obs) & ('RASA36' not in obs):
	#	Seeing measurement
	if __name__ == '__main__':
		with multiprocessing.Pool(processes=ncores) as pool:
			results = pool.starmap(tool_tbd.SE_seeing, zip(afzimlist, repeat(obs), repeat(ccddat), repeat(path_config), repeat(3*u.arcsecond), repeat(0.95), repeat(True)))
	#	Remove cosmic-ray
	if __name__ == '__main__':
		with multiprocessing.Pool(processes=ncores) as pool:
			results = pool.starmap(calib.cr_removal, zip(afzimlist, outimlist, repeat(gain), repeat(rdnoise)))
else:
	print('Skip Seeing measurement & CR remove processes for {}'.format(obs))
	for inim, outim in zip(afzimlist, outimlist):
		cpcom = 'cp {} {}'.format(inim, outim)
		print(cpcom)
		os.system(cpcom)
timetbl['status'][timetbl['process']=='cr_removal'] = True
timetbl['time'][timetbl['process']=='cr_removal'] = int(time.time() - st_)
#------------------------------------------------------------
#	FILE NAME CHANGE
#------------------------------------------------------------
# for inim in sorted(glob.glob(path_data+'/afz*.fits')):
# for inim in sorted(glob.glob(path_data+'/crafz*.fits')):
# 	obj = (fits.getheader(inim)['object']).upper()
# 	#	consider IMSNG galaxy
# 	if len(astrometryfailist) != 0:
# 		print('Astrometry fail list:', astrometryfailist)
# 		if (obj in alltbl['obj']) & (inim in astrometryfailist):
# 			robj, sep = tool_tbd.imsng_name_correction(inim, alltbl, radius=ccdinfo['fov']*u.arcmin)
# 		else:
# 			pass
# 	else:
# 		pass
# 	calib.fnamechange(inim, obs)

fov = ccdinfo['fov']*u.arcmin 
crafzimlist = []
for ims in ('{}/crafz*.fits'.format(path_data), '{}/crafz*.fit'.format(path_data), '{}/crafz*.fts'.format(path_data)):
	crafzimlist.extend(sorted(glob.glob(ims)))
# for inim in sorted(glob.glob('{}/crafz*.fits'.format(path_data))): 
for inim in crafzimlist:
	obj = fits.getheader(inim)['object'] 
	#	Modify incorrect object header
	if (inim.replace('crafz', 'afz') in astrometryfailist) & (obj in alltbl['obj']): 
		robj, sep = tool_tbd.imsng_name_correction(inim, alltbl, radius=fov) 
	else:
		pass
	calib.fnamechange(inim, obs)

caliblist = sorted(glob.glob(path_data+'/Calib*.fits'))
ic_cal = ImageFileCollection(path_data, glob_include='Calib*0.fits', keywords='*')
os.system('chmod 777 {}'.format(path_data))
os.system('chmod 777 {}/*'.format(path_data))
#	Calib-*.fits TO SAVE PATH
f = open(path_data+'/object.txt', 'a')
f.write('obs obj dateobs filter exptime\n')
for inim in caliblist:
	img = os.path.basename(inim)
	part = img.split('-')
	line = '{} {} {} {} {}\n'.format(part[1], part[2], part[3]+'T'+part[4], part[5], part[6])
	print(line)
	f.write(line)
f.close()

#	DATA FOLDER TO SAVE PATH
# os.system('rm {}/afz*.fits {}/fz*.fits'.format(path_data, path_data))
os.system(f'rm {path_data}/*fz*.f*')
os.system(f'rm -rf {path_save}/{os.path.basename(path_data)}')
plt.close('all')
#------------------------------------------------------------
#	Defringe for LOAO I-band images
#------------------------------------------------------------
st_ = time.time()
if (obs == 'LOAO') & ('I' in ic_cal.filter(imagetyp='object').summary['filter']):
	dfim = '/home/paek/qsopy/fringe/LOAO/fringe_i_ori.fits'
	dfdat = '/home/paek/qsopy/fringe/LOAO/fringe_i.dat'
	dfimlist = []
	for inim in ic_cal.filter(imagetyp='object', filter='I').summary['file']:
		# dfimlist.append(calib.defringe(str(inim), dfim, dfdat))
		dfedim = calib.defringe(str(inim), dfim, dfdat)
		mvcom = 'mv {} {}'.format(dfedim, inim)
		print(mvcom)
		os.system(mvcom)
	# tdict['defringe'] = time.time() - st - tdict[list(tdict.keys())[-1]]
else:
	print('No images to defringe')
	pass
timetbl['status'][timetbl['process']=='defringe'] = True
timetbl['time'][timetbl['process']=='defringe'] = int(time.time() - st_)
#------------------------------------------------------------
print('='*60)
print('Calibration IS DONE.\t('+str(int(time.time() - starttime))+' sec)')
print('='*60)
#------------------------------------------------------------
#	Photometry for single images
#------------------------------------------------------------
st_ = time.time()
print('#\tPhotometry')
# path_data = '{}/{}'.format(path_obs, os.path.basename(path))
path_infile = '{}/{}'.format(path_data, os.path.basename(path_default_gphot))
path_new_gphot = '{}/gphot.config'.format(os.path.dirname(path_infile))
#	Read default photometry configuration
os.system('cp {} {}'.format(path_default_gphot, path_new_gphot))
f = open(path_default_gphot, 'r')
lines = f.read().splitlines()
f.close()
#	Write photometry configuration
g = open(path_new_gphot, 'w')
for line in lines:
	if 'imkey' in line:
		line = '{}\t{}/C*0.fits'.format('imkey', path_data)
	else:
		pass
	g.write(line+'\n')
g.close()
#	Execute
if obs == 'DOAO':
	path_phot = path_phot_sg
else:
	path_phot = path_phot_mp
com = 'python {} {} {}'.format(path_phot, path_data, ncores)
print(com)
os.system(com)
# tdict['photometry'] = time.time() - st - tdict[list(tdict.keys())[-1]]
timetbl['status'][timetbl['process']=='photometry'] = True
timetbl['time'][timetbl['process']=='photometry'] = int(time.time() - st_)

#------------------------------------------------------------
#	Image registering & combine
#------------------------------------------------------------
st_ = time.time()
print('IMAGE REGISTERING & COMBINE')
combined_images = []
step = (1/24/60)*60	# 1 hour
ic_cal_phot = ImageFileCollection(path_data, glob_include='Calib*0.fits', keywords='*')
calist = sorted(glob.glob('{}/Calib*.fits'.format(path_data)))


objlist = sorted(list(set(ic_cal_phot.summary['object'])))
filterlist = sorted(list(set(ic_cal_phot.summary['filter'])))
# obj = 'NGC3147'
# filte = 'R'
for obj in objlist:
	for filte in filterlist:
		imlist_tmp = sorted(glob.glob('{}/Calib*-{}-*-{}-*.fits'.format(path_data, obj, filte)))
		if len(imlist_tmp) == 0:
			pass
		elif len(imlist_tmp) == 1:
			inim = imlist_tmp[0]
			comim = inim.replace('.fits', '.com.fits')
			cpcom = f'cp {inim} {comim}'
			print(cpcom)
			os.system(cpcom)
		else:
			print(obj, filte)
			# ic_part = sorted(glob.glob('{}/Calib*{}*{}*.fits'.format(path_data, obj, filte)))

			jds = np.array([fits.getheader(inim)['jd'] for inim in imlist_tmp])
			delts = jds - np.min(jds)

			grouplist = []
			grouplists = []

			i = 0
			for i in range(len(delts)):
				#	Initial setting
				if i == 0:
					t0 = delts[i]
				#	Add last group to grouplists
				elif i == len(delts)-1:
					grouplists.append(grouplist)
				t1 = delts[i]
				# print(t0, t1)
				dif = np.abs(t0-t1)

				if dif < step:
					grouplist.append(imlist_tmp[i])
				#	Generate new group
				else:
					grouplists.append(grouplist)
					grouplist = [imlist_tmp[i]]	
				t0 = t1
				
			for group in grouplists:
				print('-'*60)
				images_to_align = group
				ref_image = images_to_align[0]
				# print(images_to_align, ref_image)
				for inim in images_to_align: print(inim)
				try:
					outim = tool_tbd.imcombine_routine(images_to_align, ref_image)
					combined_images.append(outim)
				except:
					print('Fail to image align & combine routine.')
					print(images_to_align)
					pass


'''
#	old version (~21.04.06)
# obj = 'NGC3147'
# filte = 'R'
for obj in objlist:
	for filte in filterlist:
		imlist_tmp = sorted(glob.glob('{}/Calib*-{}-*-{}-*.fits'.format(path_data, obj, filte)))
		if len(imlist_tmp) == 0:
			pass
		else:
			print(obj, filte)
			# ic_part = sorted(glob.glob('{}/Calib*{}*{}*.fits'.format(path_data, obj, filte)))

			jds = np.array([fits.getheader(inim)['jd'] for inim in imlist_tmp])
			delts = jds - np.min(jds)

			#	Boundaries for grouping
			bds_delt = [0]
			delt0 = delts[0]

			for delt1 in delts:
				delt = np.abs(delt0 - delt1)
				if delt < step:
					pass
				else:
					bds_delt.append(delt1)
				delt0 = delt1
			bds_delt.append(delt1)
			#	Grouping and imcombine routine
			for i in range(len(bds_delt)-1):
				bd_lo = bds_delt[i]
				bd_up = bds_delt[i+1]
				indx_com = np.where((bd_lo <= delts) & (delts <= bd_up))
				# images_to_align = list(ic_part.summary['file'][indx_com])
				images_to_align = tool_tbd.npstr2str(np.array(imlist_tmp)[indx_com])
				ref_image = imlist_tmp[0]
				try:
					outim = tool_tbd.imcombine_routine(images_to_align, ref_image)
					combined_images.append(outim)
				except:
					print('Fail to image align & combine routine.')
					print(images_to_align)
					pass
'''
rmcom = 'rm {}/*Calib*gregister.fits'.format(path_data)
print(rmcom)
os.system(rmcom)
# tdict['imagecombine'] = time.time() - st - tdict[list(tdict.keys())[-1]]
timetbl['status'][timetbl['process']=='image_stack'] = True
timetbl['time'][timetbl['process']=='image_stack'] = int(time.time() - st_)
#------------------------------------------------------------
#	Photometry for combined images
#------------------------------------------------------------
st_ = time.time()
#	Write photometry configuration
h = open(path_new_gphot, 'w')
for line in lines:
	if 'imkey' in line:
		line = '{}\t{}/C*com.fits'.format('imkey', path_data)
	else:
		pass
	h.write(line+'\n')
h.close()
#	Execute
path_phot = path_phot_mp
com = 'python {} {} {}'.format(path_phot, path_data, ncores)
print(com)
os.system(com)
# tdict['photometry_com'] = time.time() - st - tdict[list(tdict.keys())[-1]]
timetbl['status'][timetbl['process']=='photometry_com'] = True
timetbl['time'][timetbl['process']=='photometry_com'] = int(time.time() - st_)

ic_com_phot = ImageFileCollection(path_data, glob_include='Calib*com.fits', keywords='*')	
#	Summary
print('Draw observation summary plots')
# for filte in list(set(ic_cal_phot.summary['filter'])):
for filte in filterlist:
	try:
		tool_tbd.obs_summary(filte, ic_cal_phot, ic_com_phot, path_save=path_data)
	except:
		print('Fail to make summary plots.')
		pass
	plt.close('all')
#------------------------------------------------------------
#	Image subtraction
#------------------------------------------------------------
# print('IMAGE SUBTRACTION')
# images_to_subtract = []
# images_to_ref = []
# ds9comlist = []
# for inim in combined_images:
# 	hdr = fits.getheader(inim)
# 	obj = hdr['object']
# 	filte = hdr['filter']
# 	path_refim = '/data3/paek/factory/ref_frames/{}'.format(obs)
# 	refimlist = glob.glob('{}/Ref*{}*{}*.fits'.format(path_refim, obj, filte))
# 	if len(refimlist) > 0:
# 		# refim = refimlist[0]
# 		# subim, ds9com = tool_tbd.subtraction_routine(inim, refim)
# 		# subtracted_images.append(subim)
# 		# ds9comlist.append(ds9com)
# 		images_to_subtract.append(inim)
# 		images_to_ref.append(refimlist[0])
# 	else:
# 		print('There is no reference image for {}'.format(os.path.basename(inim)))
# 		pass
# #	Image subtraction
# if len(images_to_subtract) > 0:
# 	if __name__ == '__main__':
# 		with multiprocessing.Pool(processes=ncores) as pool:
# 			results = pool.starmap(tool_tbd.subtraction_routine, zip(images_to_subtract, images_to_ref))
# else:
# 	pass
# subtracted_images = sorted(glob.glob('{}/hd*.fits'.format(path_data)))

print('IMAGE SUBTRACTION')
subtracted_images = []
ds9comlist = []
for inim in combined_images:
	hdr = fits.getheader(inim)
	# obs = os.path.basename(inim).split('-')[1]
	# obs = 'LOAO'
	obj = hdr['object']
	filte = hdr['filter']
	path_refim = '/data3/paek/factory/ref_frames/{}'.format(obs)
	refimlist = glob.glob('{}/Ref*{}*{}*.fits'.format(path_refim, obj, filte))
	if len(refimlist) > 0:
		refim = refimlist[0]
		subim, ds9com = tool_tbd.subtraction_routine(inim, refim)
		subtracted_images.append(subim)
		ds9comlist.append(ds9com)
	else:
		print('There is no reference image for {}'.format(os.path.basename(inim)))
		pass
rmcom = 'rm {}/*Ref*gregister.fits'.format(path_data)
print(rmcom)
os.system(rmcom)
# tdict['subtraction'] = time.time() - st - tdict[list(tdict.keys())[-1]]
timetbl['status'][timetbl['process']=='subtraction'] = True
timetbl['time'][timetbl['process']=='subtraction'] = int(time.time() - st_)
#------------------------------------------------------------
#	Photometry for subtracted images
#------------------------------------------------------------
st_ = time.time()
#	Write photometry configuration
s = open(path_new_gphot, 'w')
for line in lines:
	if 'imkey' in line:
		line = '{}\t{}/hd*com.fits'.format('imkey', path_data)
	else:
		pass
	if 'photfraction' in line:
		line = '{}\t{}'.format('photfraction', 1.0)
	else:
		pass
	if 'DETECT_MINAREA' in line:
		line = '{}\t{}'.format('DETECT_MINAREA', 10)
	else:
		pass
	if 'DETECT_THRESH' in line:
		line = '{}\t{}'.format('DETECT_THRESH', 1.25)
	else:
		pass
	s.write(line+'\n')
s.close()
#	Execute
hdimlist = sorted(glob.glob('{}/hd*.fits'.format(path_data)))
if len(hdimlist) > 0:
	com = 'python {} {}'.format(path_phot_sub, path_data)
	print(com)
	os.system(com)
	# tdict['photometry_sub'] = time.time() - st - tdict[list(tdict.keys())[-1]]
else:
	print('No subtracted image.')
	pass
timetbl['status'][timetbl['process']=='photometry_sub'] = True
timetbl['time'][timetbl['process']=='photometry_sub'] = int(time.time() - st_)
#------------------------------------------------------------
#	Transient Search
#------------------------------------------------------------
st_ = time.time()
fovval = fov.value
#	Input table for transient search
tstbl = Table()
hdimlist = sorted(glob.glob(f'{path_data}/hd*com.fits'))
if len(hdimlist) != 0:
	tstbl['hdim'] = hdimlist
	tskeys = ['hdcat', 'hcim', 'inim', 'scicat']
	for key in tskeys:
		tstbl[key] = ' '*100
	tstbl['fovval'] = fovval

	for i, hdim in enumerate(hdimlist):
		hdcat = hdim.replace('.fits','.phot_sub.cat')
		hcim = hdim.replace('hdCalib', 'hcCalib')
		inim = hdim.replace('hdCalib', 'Calib')
		scicat = inim.replace('.fits', '.phot.cat')

		for key, im in zip(tskeys, [hdcat, hcim, inim, scicat]):
			tstbl[key][i] = im

	out_tstbl = f'{path_data}/transient_search.txt'
	tstbl.write(out_tstbl, format='ascii.tab', overwrite=True)

	com = f'python {path_find} {out_tstbl} {ncores}'
	print(com)
	subprocess.call(com, shell=True)		
	'''
	fovval = fov.value
	hdimlist = sorted(glob.glob(f'{path_data}/hd*com.fits'))
	# hdim = hdimlist[0]
	for hdim in hdimlist:
		hdcat = hdim.replace('.fits','.phot_sub.cat')
		hcim = hdim.replace('hdCalib', 'hcCalib')
		inim = hdim.replace('hdCalib', 'Calib')
		scicat = inim.replace('.fits', '.phot.cat')

		try:
			#	Single
			com = f'python {path_find} {hdim} {hdcat} {hcim} {scicat} {fovval} {ncores}'
			print(com)
			# os.system(com)
			subprocess.call(com, shell=True)
		except:
			print(f'Try/Except : Skip transient search process : {hdim}')
	'''
timetbl['status'][timetbl['process']=='transient_search'] = True
timetbl['time'][timetbl['process']=='transient_search'] = int(time.time() - st_)
#------------------------------------------------------------
#	Summary file
#------------------------------------------------------------
timetbl['status'][timetbl['process']=='total'] = True
timetbl['time'][timetbl['process']=='total'] = int(time.time() - st)	
timetbl.write('{}/obs.summary.log'.format(path_data), format='ascii.tab', overwrite=True)
print(timetbl)
#	Write data summary
f = open(path_data+'/obs.summary.log', 'a')
end_localtime = time.strftime('%Y-%m-%d %H:%M:%S (%Z)', time.localtime())
f.write('Pipelne start\t: {}\n'.format(start_localtime))
f.write('Pipelne end\t: {}\n'.format(end_localtime))
try:
	f.write('='*60+'\n')
	f.write('PATH :{}\n'.format(path))
	f.write('OBJECT NUMBER # :{}\n'.format(len(ic_cal.summary)))
	objkind = sorted(set(ic_cal.summary['object']))
	f.write('OBJECTS # : {}\n'.format(objkind))
	for obj in objkind:
		f.write('-'*60+'\n')
		for filte in list(set(ic_cal.summary['filter'])):
			indx_tmp = ic_cal.files_filtered(filter=filte, object=obj)
			if len(indx_tmp) > 0:
				f.write('{}\t{}\n'.format(obj, filte))
except:
	pass
f.close()
'''
#------------------------------------------------------------
#	Send E-mail
#------------------------------------------------------------
param_email = dict(
						path_machine = '/home/paek/table/mail.machine.dat',
						path_address = '/home/paek/table/mail.recivers.txt',
						subject = '[IMSNG] Pipeline is done for {}.'.format(obs),
						path_contents = '{}/obs.summary.log'.format(path_data),
						path_attach_png = '{}/obs.summary*.png'.format(path_data),							
						path_attach_txt = '{}/phot.*.dat'.format(path_data),
					)
tool_tbd.sendmail(**param_email)
'''
#------------------------------------------------------------
#	Slack message
#------------------------------------------------------------	
total_time = round(timetbl['time'][timetbl['process']=='total'].item()/60., 1)

channel = '#pipeline'
text = f'[Pipeline/{obs}] End Processing {os.path.basename(path)} Data ({nobj} objects) with {ncores} cores taking {total_time} mins'

param_slack = dict(
	token = OAuth_Token,
	channel = channel,
	text = text,
)

tool_tbd.slack_bot(**param_slack)
#------------------------------------------------------------
# tdict['total'] = time.time() - st
# tool_tbd.dict2table(tdict, '{}/{}'.format(path_data, 'times.log'))
#------------------------------------------------------------
#	Finish the process
#------------------------------------------------------------
# rmcom = 'rm {}/hc*gregister.fits {}/inv*.fits'.format(path_data, path_data)
rmcom = 'rm {}/inv*.*'.format(path_data, path_data)
print(rmcom)
os.system(rmcom)
tails = ['.transients.', '.new.', '.ref.', '.sub.', '']
for obj in objlist:
	for filte in filterlist:
		for tail in tails:
			# tail = 'transients'
			# obj = 'NGC3147'
			# filte = 'B'

			pathto = f'{path_gal}/{obj}/{obs}/{filte}'
			files = f'{path_data}/*Calib*-{obj}-*-{filte}-*{tail}*'
			nfiles = len(glob.glob(files))
			# print(files, nfiles)
			# if nfiles >0:
			# 	print(obj, filte, pathto, files, glob.glob(files)[-1])
			if nfiles !=0:
				#	Save path
				if tail == '':
					pathto = f'{path_gal}/{obj}/{obs}/{filte}'
				else:
					pathto = f'{path_gal}/{obj}/{obs}/{filte}/transients'
				#	Make path
				if (not os.path.exists(pathto)):
					os.makedirs(pathto)
				mvcom = f'mv {files} {pathto}'
				print(mvcom)
				os.system(mvcom)
'''
caliblist = sorted(glob.glob('{}/*Calib*'.format(path_data)))
for inim in caliblist:
	if ('.new.' in inim) | ('.ref.' in inim) | ('.sub.' in inim) | ('.transients.' in inim):
		calib.movecalib(inim, path_gal, 'transients')
	else:
		calib.movecalib(inim, path_gal)
'''
mvcom = f'mv {path_data} {path_save}'
os.system(mvcom)
#	WRITE LOG
f = open(path_log, 'a')
# f.write(path_raw+'/'+os.path.basename(path_data)+'\n')
# f.write('{}/{}\n'.format(path_raw, os.path.basename(path_data)))
f.write(f'{path_raw}/{os.path.basename(path_data)}\n')
f.close()
#============================================================
#	Time
#============================================================
delt = time.time() - starttime
dimen = 'seconds'
if delt > 60.:
	delt = delt/60.
	dimen = 'mins'
if delt > 60.:
	delt = delt/60.
	dimen = 'hours'
print('ALL PROCESS IS DONE.\t({} {})'.format(round(delt, 3), dimen))
print(obs, newlist)