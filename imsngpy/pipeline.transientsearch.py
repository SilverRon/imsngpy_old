#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#============================================================
#%%
#	Library
#------------------------------------------------------------
import time
import os
import sys
from datetime import date
import numpy as np
import subprocess
#	IMSNGpy modules
sys.path.append('/home/paek/imsngpy')
from tableutil import *
from preprocess import *
from misc import *
from phot import *
from util import *
from query import *
#	Astropy
from astropy.io import ascii
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.coordinates import SkyCoord
from astropy.table import hstack
from astropy.wcs import WCS
from astropy import units as u
#	Bottle Neck
import bottleneck as bn
#	Multiprocess tools
from itertools import repeat
import multiprocessing
st_ = time.time()
#============================================================
#%%
#	Path
#------------------------------------------------------------
path_base = '/home/paek/imsngpy'
path_table = f'{path_base}/table'
path_config = f'{path_base}/config'
path_gphot = f'{path_config}/gphot.config'
#	Table
ccdtbl = ascii.read(f'{path_table}/ccd.tsv') 
#	gphot Configuration
if os.path.exists(path_gphot) == True:
	gphot_dict = file2dict(path_gphot)
else:
	#	Default gpphot.config
	print('[NOTICE] There is no gphot.config. Use default configuration.')
	gphot_dict = {
		'imkey': './Calib*com.fits',
		'photfraction': '0.75',
		'refcatname': 'PS1',
		'refqueryradius': '1',
		'refmaglower': '14',
		'refmagupper': '18',
		'refmagerupper': '0.05',
		'inmagerupper': '0.05',
		'flagcut': '0',
		'sepfrac': '0.5',
		'DETECT_MINAREA': '5',
		'DETECT_THRESH': '3.0',
		'DEBLEND_NTHRESH': '64',
		'DEBLEND_MINCONT': '0.0001',
		'BACK_SIZE': '64',
		'BACK_FILTERSIZE': '3',
		'BACKPHOTO_TYPE': 'LOCAL',
		'check': 'False'
		}
#------------------------------------------------------------
#	Configuration files
#------------------------------------------------------------
#	gpphot
prefix_gp = 'gpphot'
path_param = f'{path_config}/{prefix_gp}.param'
path_conv = f'{path_config}/{prefix_gp}.conv'
path_nnw = f'{path_config}/{prefix_gp}.nnw'
path_conf = f'{path_config}/{prefix_gp}.sex'
#	invert
prefix_iv = 'invert'
path_param_iv = f'{path_config}/{prefix_iv}.param'
path_conv_iv = f'{path_config}/{prefix_iv}.conv'
path_nnw_iv = f'{path_config}/{prefix_iv}.nnw'
path_conf_iv = f'{path_config}/{prefix_iv}.sex'
#============================================================
#	Input
#------------------------------------------------------------
# tstablename = sys.argv[1]
tstablename = '/data3/paek/factory/loao/test_fast/transient_search.ecsv'
print(f'Transient Search Table  : {tstablename}')
tstbl = ascii.read(tstablename)
try:
	ncore = int(sys.argv[2])
except:
	ncore = 1
print(f'NCORE  : {ncore}')
print(f"{'='*60}\n")
#------------------------------------------------------------
#	Function
#------------------------------------------------------------
def routine_phot(inim,):
	print(f'PHOTOMETRY START for {os.path.basename(inim)}')
	#
	global path_param
	global path_conv
	global path_nnw
	global path_conf
	#------------------------------------------------------------
	#	Image info
	#------------------------------------------------------------
	hdr = fits.getheader(inim)
	#	IMAGE PHYSICAL CENTER
	a = hdr['naxis1']/2.
	b = hdr['naxis2']/2.
	#	Small squre based on frac
	frac = float(gphot_dict['photfraction'])
	a_ = a*np.sqrt(frac)
	b_ = b*np.sqrt(frac)
	obj = hdr['OBJECT']
	filte = hdr['FILTER']
	#------------------------------------------------------------
	#	CCD INFO
	#------------------------------------------------------------
	if ('OBSERVAT' in hdr.keys()) & ('CCDNAME' in hdr.keys()):
		obs = hdr['OBSERVAT']
		ccd = hdr['CCDNAME']
	else:
		#	observatory from filename
		obs = os.path.basename(inim).split('-')[1]
		if '_' not in obs:
			ccd = ''
		else:
			obs = obs.split('_')[0]
			ccd = obs.split('_')[1]
	indx_ccd = np.where(
		(ccdtbl['obs']==obs) &
		(ccdtbl['ccd']==ccd)
	)
	
	# print(f"""{'-'*60}\n#\tCCD INFO\n{'-'*60}""")
	gain = ccdtbl['gain'][indx_ccd][0]*(u.electron/u.adu)
	rdnoise = ccdtbl['readnoise'][indx_ccd][0]*(u.electron)
	pixscale = ccdtbl['pixelscale'][indx_ccd][0]*(u.arcsec/u.pixel)
	fov = ccdtbl['foveff'][indx_ccd][0]*(u.arcmin)
	# print(f"""GAIN : {gain}\nREAD NOISE : {rdnoise}\nPIXEL SCALE : {pixscale}\nEffective FoV : {fov}""")
	#------------------------------------------------------------
	#	Single
	if ('SEEING' in hdr.keys()) & ('PEEING' in hdr.keys()):
		seeing = hdr['SEEING']*u.arcsec
		peeing = hdr['PEEING']*u.pix
	else:
		print('No seeing information on the header. Calculate ')
		#	If no seeing info, get it with simple conf
		prefix_sp = 'simple'
		path_conf = f'{path_config}/{prefix_sp}.sex'
		path_param = f'{path_config}/{prefix_sp}.param'
		path_nnw = f'{path_config}/{prefix_sp}.nnw'
		path_conv = f'{path_config}/{prefix_sp}.conv'
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

	#------------------------------------------------------------
	#%%
	#	OPTIMIZED APERTURE FROM SNR CURVE
	#------------------------------------------------------------
	aper_opt = hdr['APER']*pixscale.value
	#------------------------------------------------------------
	#%%
	#	APERTURE DICTIONARY
	#------------------------------------------------------------
	aper_dict = dict(
		MAG_AUTO = dict(
			size=0.0,
			errkey='MAGERR_AUTO',
			suffix='AUTO',
			comment='',
		),
		MAG_APER = dict(
			size=aper_opt/pixscale.value,
			errkey='MAGERR_APER',
			suffix='APER',
			comment='Best aperture diameter based on SNR curve [pix]',
		),
		MAG_APER_1 = dict(
			size=2*0.6731*peeing.value,
			errkey='MAGERR_APER_1',
			suffix='APER_1',
			comment='2*0.6731*seeing (gaussian profile) [pix]',
		),
		MAG_APER_2 = dict(
			size=2*peeing.value,
			errkey='MAGERR_APER_2',
			suffix='APER_2',
			comment='2*seeing aperture diameter [pix]',
		),
		MAG_APER_3 = dict(
			size=3*peeing.value,
			errkey='MAGERR_APER_3',
			suffix='APER_3',
			comment='3*seeing aperture diameter [pix]',
		),
		MAG_APER_4 = dict(
			size=(3*u.arcsec/pixscale).value,
			errkey='MAGERR_APER_4',
			suffix='APER_4',
			comment='Fixed 3\" aperture diameter [pix]',
		),
		MAG_APER_5 = dict(
			size=(5*u.arcsec/pixscale).value,
			errkey='MAGERR_APER_5',
			suffix='APER_5',
			comment='Fixed 5\" aperture diameter [pix]',
		),
	)
	#------------------------------------------------------------
	#%%
	#	SET APERTURES
	#	MAG_AUTO, MAG_APER, MAG_APER_1, MAG_APER_2, MAG_APER_3, MAG_APER_4, MAG_APER_5
	apertures_phot = []
	for magkey in list(aper_dict.keys()):
		if magkey != 'MAG_AUTO':
			apertures_phot.append(aper_dict[magkey]['size'])
	apertures_phot = np.array(apertures_phot)
	apertures_input = ','.join([str(d) for d in apertures_phot])
	#	Run source extractor for photometry
	outcat = f'{os.path.splitext(inim)[0]}.cat'
	param_insex = dict(
		#------------------------------
		#	CATALOG
		#------------------------------
		CATALOG_NAME = outcat,
		#------------------------------
		#	CONFIG FILES
		#------------------------------
		CONF_NAME = path_conf,
		PARAMETERS_NAME = path_param,
		FILTER_NAME = path_conv,    
		STARNNW_NAME = path_nnw,
		#------------------------------
		#	PHOTOMETRY
		#------------------------------
		PHOT_APERTURES = apertures_input,
		GAIN = str(gain.value),
		PIXEL_SCALE = str(pixscale.value),
		#------------------------------
		#	STAR/GALAXY SEPARATION
		#------------------------------
		SEEING_FWHM = str(seeing.value),
		#------------------------------
		#	EXTRACTION
		#------------------------------
		DETECT_MINAREA = gphot_dict['DETECT_MINAREA'],
		DETECT_THRESH = gphot_dict['DETECT_THRESH'],
		DEBLEND_NTHRESH = gphot_dict['DEBLEND_NTHRESH'],
		DEBLEND_MINCONT = gphot_dict['DEBLEND_MINCONT'],
		#------------------------------
		#	BACKGROUND
		#------------------------------
		BACK_SIZE = '128',
		BACK_FILTERSIZE = '10',
		BACKPHOTO_TYPE = 'LOCAL',
		#------------------------------
		#	CHECK IMAGE
		#------------------------------
		CHECKIMAGE_TYPE = 'NONE',
		CHECKIMAGE_NAME = 'check.fits',
	)
	# os.system(sexcom(inim, param_insex))
	# rawtbl = ascii.read(outcat)
	
	sexout = subprocess.getoutput(sexcom(inim, param_insex))
	line = [s for s in sexout.split('\n') if 'RMS' in s]
	skyval, skysig = float(line[0].split('Background:')[1].split('RMS:')[0]), float(line[0].split('RMS:')[1].split('/')[0])
	rawtbl = ascii.read(outcat)
	
	#------------------------------------------------------------
	#	GENERAL INFO. REGARDLESS OF APERTURE
	#------------------------------------------------------------
	#	VALUES
	rsky = 1e2*skysig/skyval
	#	PUT HEADER ON THE IMAGE
	fits.setval(inim, 'AUTHOR', value='Gregory S.H. Paek', comment='PHOTOMETRY AUTHOR')
	fits.setval(inim, 'PHTDATE', value=date.today().isoformat(), comment='PHOTOMETRY DATE [KOR]')
	fits.setval(inim, 'SKYSIG', value=round(skysig, 3), comment='STD of BACKGROUND [pixel]')
	fits.setval(inim, 'SKYVAL', value=round(skyval, 3), comment='MEDIAN of BACKGROUND [pixel]')
	fits.setval(inim, 'RSKY', value=round(rsky, 3), comment='RATIO of SKYSIG/SKYVAL [%]')

	#============================================================
	#%%
	#	ZP CALCULATION
	#------------------------------------------------------------
	# n = 0
	# inmagkey = list(aper_dict.keys())[n]
	for n, inmagkey in enumerate(list(aper_dict.keys())):
		inerrkey = aper_dict[inmagkey]['errkey']
		inzpkey = f"ZP_{aper_dict[inmagkey]['suffix']}"
		inul3key = f"UL3{aper_dict[inmagkey]['suffix'].replace('APER', '')}"
		inul5key = f"UL5{aper_dict[inmagkey]['suffix'].replace('APER', '')}"
		hdrzpkey = inzpkey.replace('_APER', '')
		hdrzperkey = hdrzpkey.replace('ZP', 'ZPER')
		hdrnstdkey = hdrzpkey.replace('ZP', 'NZP')

		#	MERGED TABLE
		rawtbl[inmagkey.lower()] = rawtbl[inmagkey]+hdr[hdrzpkey]
		rawtbl[aper_dict[inmagkey]['errkey'].lower()] = sqsum(rawtbl[inerrkey], hdr[hdrzperkey])

	rawtbl.write(outcat, format='ascii.ecsv', overwrite=True)
	delt = time.time()-st_
	print(f"DONE ({round(delt, 1)} sec)\n")
#------------------------------------------------------------
def mp_phot(imlist, ncore=4):
	if __name__ == '__main__':
		#	Fixed the number of cores (=4)
		with multiprocessing.Pool(processes=ncore) as pool:
			results = pool.starmap(
				routine_phot,
				zip(
					imlist,
					)
				)
#------------------------------------------------------------
def phot(imlist):
	for i, inim in enumerate(imlist):
		print(f"[{i+1}/{len(imlist)}] {os.path.basename(inim)}")
		routine_phot(inim)
#------------------------------------------------------------
def invert_data(inim, outim):
	data, hdr = fits.getdata(inim, header=True)
	fits.writeto(outim, data*-1, header=hdr, overwrite=True)
#------------------------------------------------------------
#	Photometry
mp_phot(tstbl['sub'], ncore=ncore)
#------------------------------------------------------------
n=0
sciim = tstbl['sci'][n]
refim = tstbl['ref'][n]
subim = tstbl['sub'][n]

hdr = fits.getheader(sciim)
c_cent = SkyCoord(hdr['CRVAL1'], hdr['CRVAL2'], unit=u.deg)
dateobs = hdr['DATE-OBS']
epoch = Time(dateobs, format='isot')
seeing = hdr['SEEING']*u.arcsec
ccdinfo = get_ccdinfo(sciim, ccdtbl)

isciim = f"{os.path.dirname(sciim)}/inv{os.path.basename(sciim)}"
isubim = f"{os.path.dirname(subim)}/inv{os.path.basename(subim)}"
irefim = f"{os.path.dirname(refim)}/inv{os.path.basename(refim)}"

invert_data(sciim, isciim)
invert_data(refim, irefim)
invert_data(subim, isubim)
#------------------------------------------------------------
scicat = f'{os.path.splitext(sciim)[0]}.cat'
#------------------------------------------------------------
subcat = f'{os.path.splitext(subim)[0]}.cat'
subtbl = ascii.read(subcat)
c_sub = SkyCoord(subtbl['ALPHA_J2000'], subtbl['DELTA_J2000'])
#------------------------------------------------------------
isubcat = f'{os.path.splitext(isubim)[0]}.cat'
ccdinfo = get_ccdinfo(sciim, ccdtbl)
param_insex4isub = dict(
	#------------------------------
	#	CATALOG
	#------------------------------
	CATALOG_NAME = isubcat,
	#------------------------------
	#	CONFIG FILES
	#------------------------------
	CONF_NAME = path_conf_iv,
	PARAMETERS_NAME = path_param_iv,
	FILTER_NAME = path_conv_iv,    
	STARNNW_NAME = path_nnw_iv,
	#------------------------------
	#	PHOTOMETRY
	#------------------------------
	GAIN = str(ccdinfo['gain'].value),
	PIXEL_SCALE = str(ccdinfo['pixscale'].value),
	#------------------------------
	#	STAR/GALAXY SEPARATION
	#------------------------------
	SEEING_FWHM = str(seeing.value),
	#------------------------------
	#	EXTRACTION
	#------------------------------
	DETECT_MINAREA = gphot_dict['DETECT_MINAREA'],
	DETECT_THRESH = gphot_dict['DETECT_THRESH'],
	DEBLEND_NTHRESH = gphot_dict['DEBLEND_NTHRESH'],
	DEBLEND_MINCONT = gphot_dict['DEBLEND_MINCONT'],
	#------------------------------
	#	BACKGROUND
	#------------------------------
	BACK_SIZE = '128',
	BACK_FILTERSIZE = '10',
	BACKPHOTO_TYPE = 'LOCAL',
	#------------------------------
	#	CHECK IMAGE
	#------------------------------
	CHECKIMAGE_TYPE = 'NONE',
	CHECKIMAGE_NAME = 'check.fits',
)
os.system(sexcom(isubim, param_insex4isub))
isubtbl = ascii.read(isubcat)
c_isub = SkyCoord(isubtbl['ALPHA_J2000'], isubtbl['DELTA_J2000'])
#------------------------------------------------------------
irefcat = f'{os.path.splitext(irefim)[0]}.cat'
param_insex4iref = dict(
	#------------------------------
	#	CATALOG
	#------------------------------
	CATALOG_NAME = irefcat,
	#------------------------------
	#	CONFIG FILES
	#------------------------------
	CONF_NAME = path_conf_iv,
	PARAMETERS_NAME = path_param_iv,
	FILTER_NAME = path_conv_iv,    
	STARNNW_NAME = path_nnw_iv,
	#------------------------------
	#	PHOTOMETRY
	#------------------------------
	GAIN = str(ccdinfo['gain'].value),
	PIXEL_SCALE = str(ccdinfo['pixscale'].value),
	#------------------------------
	#	STAR/GALAXY SEPARATION
	#------------------------------
	SEEING_FWHM = str(seeing.value),
	#------------------------------
	#	EXTRACTION
	#------------------------------
	DETECT_MINAREA = gphot_dict['DETECT_MINAREA'],
	DETECT_THRESH = gphot_dict['DETECT_THRESH'],
	DEBLEND_NTHRESH = gphot_dict['DEBLEND_NTHRESH'],
	DEBLEND_MINCONT = gphot_dict['DEBLEND_MINCONT'],
	#------------------------------
	#	BACKGROUND
	#------------------------------
	BACK_SIZE = '128',
	BACK_FILTERSIZE = '10',
	BACKPHOTO_TYPE = 'LOCAL',
	#------------------------------
	#	CHECK IMAGE
	#------------------------------
	CHECKIMAGE_TYPE = 'NONE',
	CHECKIMAGE_NAME = 'check.fits',
)
os.system(sexcom(irefim, param_insex4iref))
ireftbl = ascii.read(irefcat)
c_iref = SkyCoord(ireftbl['ALPHA_J2000'], ireftbl['DELTA_J2000'])
#------------------------------------------------------------
indx_isub, sep_isub, _ = c_sub.match_to_catalog_sky(c_isub)
indx_iref, sep_iref, _ = c_sub.match_to_catalog_sky(c_iref)

subtbl['sep_isub'] = sep_isub
subtbl['sep_iref'] = sep_iref

plt.close('all')
plt.hist(sep_isub.arcsec, bins=np.arange(0, 10+0.5, 0.5), alpha=0.5, label='sub')
plt.hist(sep_iref.arcsec, bins=np.arange(0, 10+0.5, 0.5), alpha=0.5, label='ref')
plt.axvline(x=seeing.value, ls='--', color='grey', alpha=0.5)
plt.legend(fontsize=20)
plt.show()

print(len(subtbl[sep_isub<seeing*0.5]))
print(len(subtbl[sep_iref<seeing*0.5]))

#
w_ref = WCS(refim)
x, y = w_ref.world_to_pixel(c_sub)
plt.plot(x, y, marker='o', ls='none', c='k', alpha=0.25)
#%%


numbers = np.arange(0, 11+1, 1)	#	flag 0-11
for num in numbers: subtbl[f'flag_{num}'] = False
subtbl['flag'] = False
#------------------------------------------------------------
#	flag 0
#------------------------------------------------------------
#	Skybot query
from astroquery.imcce import Skybot
frac_sb = 5.0
try:
	sbtbl = Skybot.cone_search(c_cent, ccdinfo['fov'], epoch)
	c_sb = SkyCoord(sbtbl['RA'], sbtbl['DEC'])
	sbtbl['sep'] = c_cent.separation(c_sb).to(u.arcmin)
	
	#	Skybot matching
	indx_sb, sep_sb, _ = c_sub.match_to_catalog_sky(c_sb)
	subtbl['sep_skybot'] = sep_sb
	subtbl['flag_0'][
		(sep_sb<seeing*frac_sb)
		] = True
except:
	#	RuntimeError: No solar system object was found in the requested FOV
	print(f"No solar system object was found in the requested FOV ({ccdinfo['fov']})")
	pass
#------------------------------------------------------------
#	flag 1
#------------------------------------------------------------
frac_inv = 1.0
if len(isubtbl)>0:
	#	Matching with inverted images
	indx_isub, sep_isub, _ = c_sub.match_to_catalog_sky(c_isub)
	subtbl['flag_1'][
		(sep_isub<seeing*frac_inv)
		] = True
else:
	print('Inverted subtraction image has no source. ==> pass flag1')
	pass
#------------------------------------------------------------
#	flag 2
#------------------------------------------------------------
if len(ireftbl)>0:
	#	Coordinate
	c_iref = SkyCoord(ireftbl['ALPHA_J2000'], ireftbl['DELTA_J2000'], unit=u.deg)
	#	Matching with inverted images
	indx_iref, sep_iref, _ = c_sub.match_to_catalog_sky(c_iref)
	subtbl['flag_2'][
		(sep_iref<seeing*frac_inv)
		] = True
else:
	print('Inverted reference image has no source. ==> pass flag2')
	pass
#------------------------------------------------------------
#	SEtractor criterion
#------------------------------------------------------------
#	flag 3
#------------------------------------------------------------
#	Sources @edge
frac = 0.9
#	IMAGE PHYSICAL CENTER
a = hdr['naxis1']/2.
b = hdr['naxis2']/2.
#	Small squre based on frac
a_ = a*np.sqrt(frac)
b_ = b*np.sqrt(frac)
subtbl['flag_3'][
	(
		(subtbl['X_IMAGE']>a-a_) |
		(subtbl['X_IMAGE']<a+a_) |
		(subtbl['Y_IMAGE']>b-b_) |
		(subtbl['Y_IMAGE']<b+b_)
		)
	] = True
#------------------------------------------------------------
#	flag 4
#------------------------------------------------------------
#	More than 5 sigma signal
subtbl['flag_4'][
	(subtbl['mag_aper']>hdr['ul5'])
	] = True
#	Empirical criterion
#------------------------------------------------------------
#	flag 5
#------------------------------------------------------------
subtbl['ratio_ellipticity'] = subtbl['ELLIPTICITY']/hdr['ELLIP']
subtbl['ratio_elongation'] = subtbl['ELONGATION']/hdr['ELONG']

subtbl['flag_5'][
	(subtbl['ratio_ellipticity'] > 5)
	] = True
#------------------------------------------------------------
#	flag 6
#------------------------------------------------------------
subtbl['flag_6'][
	(subtbl['FLAGS'] > 4.0)
	] = True
#------------------------------------------------------------
#	flag 7
#------------------------------------------------------------
seeing_up = 3.0
seeing_lo = 0.5

subtbl['ratio_seeing'] = subtbl['FWHM_WORLD'].to(u.arcsec)/seeing
subtbl['flag_7'][
	(subtbl['ratio_seeing']>seeing_up) |
	(subtbl['ratio_seeing']<seeing_lo)
	] = True
#------------------------------------------------------------
#	flag 8
#------------------------------------------------------------
subtbl['flag_8'][
	(subtbl['BACKGROUND']<-50) |
	(subtbl['BACKGROUND']>+50)
	] = True
#------------------------------------------------------------
#	flag 9
#------------------------------------------------------------
scitbl = ascii.read(scicat)
scitbl = scitbl[
	(scitbl['FLAGS']==0) &
	(scitbl['CLASS_STAR']>0.5)
]

aperdict = {
	'mag_aper':'SNR_curve',
	'mag_aper_1':'Best_Aperture',
	'mag_aper_2':'2seeing',
	'mag_aper_3':'3seeing',
	'mag_aper_4':'3arcsec',
	'mag_aper_5':'5arcsec',	
}

key0 = 'mag_aper_1'
key1 = 'mag_aper_3'
#	Sci. sources magnitude diff.
scidelm = scitbl[key0] - scitbl[key1]
#	Subt. sources magnitude diff.
subdelm = subtbl[key0] - subtbl[key1]
subtbl['del_mag'] = subdelm
#	MED & MAD
scidelm_med = np.median(scidelm)
scidelm_mad = get_mad(scidelm)
subtbl['del_mag_med'] = scidelm_med
subtbl['del_mag_mad'] = scidelm_mad
subtbl['N_del_mag_mad'] = np.abs((subtbl['del_mag']-subtbl['del_mag_med'])/subtbl['del_mag_mad'])
#	out
n = 10
indx_out = np.where(
	(subdelm<scidelm_med-scidelm_mad*n) |
	(subdelm>scidelm_med+scidelm_mad*n)
	)
subtbl['flag_9'][indx_out] = True
#------------------------------------------------------------
#	flag 10+11
#------------------------------------------------------------
peeing = hdr['PEEING']
skysig = hdr['SKYSIG']

nbadlist = []
ratiobadlist = []
nnulllist = []
if 'KCT' in inim:
	f = 0.05	# Set tighter criterion 
else:
	f = 0.3
for i in range(len(subtbl)):
	tx, ty = subtbl['X_IMAGE'][i], subtbl['Y_IMAGE'][i]
	bkg = subtbl['BACKGROUND'][i]
	#	Snapshot
	tsize = peeing
	y0, y1 = int(ty-tsize), int(ty+tsize)
	x0, x1 = int(tx-tsize), int(tx+tsize)
	cdata = data[y0:y1, x0:x1]
	# plt.close()
	# plt.imshow(cdata)
	crt = bkg - skysig
	cutline = cdata.size*f
	nbad = len(cdata[cdata<crt])
	try:
		ratiobad = nbad/cdata.size
	except:
		ratiobad = -99.0
	nnull = len(np.where(cdata == 1e-30)[0])
	#	Dipole
	if nbad > cutline:
		subtbl['flag_10'][i] = True
	#	HOTPANTS Null value
	if nnull != 0:
		subtbl['flag_11'][i] = True
	nbadlist.append(nbad)
	ratiobadlist.append(ratiobad)
	nnulllist.append(nnull)

subtbl['n_bad'] = nbadlist
subtbl['ratio_bad'] = ratiobadlist
subtbl['n_null'] = nnulllist
#------------------------------------------------------------
#	Final flag
#------------------------------------------------------------
flag = subtbl['flag']
n_all = len(subtbl)
for n in numbers:
	tmptbl = subtbl[subtbl[f'flag_{n}']==True] 
	print(f'flag=={n} : {len(tmptbl)} {int(100*len(tmptbl)/n_all)}%')
	flag = flag + subtbl[f'flag_{n}']
subtbl['flag'] = flag
indx_sb = np.where(subtbl['flag_0']==True)
subtbl['flag'][indx_sb] = False

outcat = hdcat.replace('phot_sub.cat', 'transients.cat')
subtbl.write(outcat, format='ascii.tab', overwrite=True)
print('-'*60)
bgstbl = subtbl[subtbl[f'flag']==True]
tctbl = subtbl[subtbl[f'flag']==False]
print(f'Filtered sources\t: {len(bgstbl)} {int(100*len(bgstbl)/n_all)}%')
print(f'Transient Candidates\t: {len(tctbl)} {int(100*len(tctbl)/n_all)}%')