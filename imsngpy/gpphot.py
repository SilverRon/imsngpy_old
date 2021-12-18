#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#============================================================
#%%
#	Library
#------------------------------------------------------------
import time
st = time.time()
import os
import sys
sys.path.append('/home/paek/imsngpy')
#	IMSNGpy modules
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
#	Bottle Neck
import bottleneck as bn
delt = time.time()-st
print(delt, 'sec')
#============================================================
#%%
#	Path
#------------------------------------------------------------
path_base = '/home/paek/imsngpy'
path_table = f'{path_base}/table'
path_config = f'{path_base}/config'
#	Table
ccdtbl = ascii.read(f'{path_table}/ccd.tsv') 

path_gphot = f'{path_config}/gphot.config'

if os.path.exists(path_gphot) == True:
	gphot_dict = file2dict(path_gphot)
else:
	#	Default gpphot.config
	print('[NOTICE] There is no gphot.config. Use default configuration.')
	gphot_dict = {
		'imkey': 'Calib*com.fits',
		'photfraction': '0.75',
		'refcatname': 'PS1',
		'refqueryradius': '1',
		'refmaglower': '14',
		'refmagupper': '18',
		'refmagerupper': '0.05',
		'inmagerupper': '0.05',
		'flagcut': '0',
		'DETECT_MINAREA': '5',
		'DETECT_THRESH': '3.0',
		'DEBLEND_NTHRESH': '64',
		'DEBLEND_MINCONT': '0.0001',
		'BACK_SIZE': '64',
		'BACK_FILTERSIZE': '3',
		'BACKPHOTO_TYPE': 'LOCAL',
		'check': 'False'
		}

#	Test image
# inim = '/data3/paek/factory/test/phot/Calib-LOAO-M99-20210421-063118-R-180.com.fits'
# inim = '/data3/paek/factory/test/phot/Calib-LOAO-M99-20210421-063001-R-60.fits'
# inim = '/data3/paek/factory/test/phot/Calib-LOAO-M99-20210421-063118-R-180.com.fits'
inim = '/data3/paek/factory/test/phot/Calib-LOAO-NGC6946-20201112-022807-R-180.com.fits'
#	Image info.
hdr = fits.getheader(inim)
filte = hdr['FILTER']
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
#------------------------------------------------------------
#	CCD INFO
import numpy as np
indx_ccd = np.where(
	(ccdtbl['obs']==obs) &
	(ccdtbl['ccd']==ccd)
)
print(f"""{'-'*60}\n#\tCCD INFO\n{'-'*60}""")
from astropy import units as u
gain = ccdtbl['gain'][indx_ccd][0]*(u.electron/u.adu)
rdnoise = ccdtbl['readnoise'][indx_ccd][0]*(u.electron)
pixscale = ccdtbl['pixelscale'][indx_ccd][0]*(u.arcsec/u.pixel)
fov = ccdtbl['foveff'][indx_ccd][0]*(u.arcmin)
print(f"""GAIN : {gain}\nREAD NOISE : {rdnoise}\nPIXEL SCALE : {pixscale}\nEffective FoV : {fov}""")
#------------------------------------------------------------
#	Single
if ('SEEING' in hdr.keys()) & ('PEEING' in hdr.keys()):
	seeing = hdr['SEEING']*u.arcsec
	peeing = hdr['PEEING']*u.pix
else:
	print('No seeing information on the header. Calculate ')
	prefix = 'simple'
	path_conf = f'{path_config}/{prefix}.sex'
	path_param = f'{path_config}/{prefix}.param'
	path_nnw = f'{path_config}/{prefix}.nnw'
	path_conv = f'{path_config}/{prefix}.conv'
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
#	Aperture steps for SNR curve
apertures = peeing.value*np.linspace(0.25, 5, 32)
apertures_input = ','.join([str(d) for d in apertures])
delt = time.time()-st
print(delt, 'sec')
#------------------------------------------------------------
#%%
#	GROWTH CURVE
#------------------------------------------------------------
prefix = 'growthcurve'
path_param_gc = f'{path_config}/{prefix}.param'
path_conv_gc = f'{path_config}/{prefix}.conv'
path_nnw_gc = f'{path_config}/{prefix}.nnw'
path_conf_gc = f'{path_config}/{prefix}.sex'
#------------------------------------------------------------
gccat = f'{os.path.splitext(inim)[0]}.gc.cat'
#	SE parameters
param_insex = dict(
	#------------------------------
	#	CATALOG
	#------------------------------
	CATALOG_NAME = gccat,
	#------------------------------
	#	CONFIG FILES
	#------------------------------
	CONF_NAME = path_conf_gc,
	PARAMETERS_NAME = path_param_gc,
	FILTER_NAME = path_conv_gc,    
	STARNNW_NAME = path_nnw_gc,
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
os.system(sexcom(inim, param_insex))
gcrawtbl = ascii.read(gccat)
#	IMAGE PHYSICAL CENTER
a = hdr['naxis1']/2.
b = hdr['naxis2']/2.
#	Small squre based on frac
frac = 0.95
a_ = a*np.sqrt(frac)
b_ = b*np.sqrt(frac)
#	Conditions for growth-curve
#	Point sources and nearby center
gctbl = gcrawtbl[
	(gcrawtbl['FLAGS']==0) &
	(gcrawtbl['CLASS_STAR']>0.9) &
	(gcrawtbl['X_IMAGE']>a-a_) & (gcrawtbl['X_IMAGE']<a+a_) &
	(gcrawtbl['Y_IMAGE']>b-b_) & (gcrawtbl['Y_IMAGE']<b+b_)
]
gctbl['APER_OPT'] = 0.0
for n in range(len(apertures)):
	if n==0:
		gctbl[f'SNR'] = gctbl[f'FLUX_APER']/gctbl[f'FLUXERR_APER']
	else:
		gctbl[f'SNR_{n}'] = gctbl[f'FLUX_APER_{n}']/gctbl[f'FLUXERR_APER_{n}']
delt = time.time()-st
print(delt, 'sec')
#------------------------------------------------------------
#%%
#	OPTIMIZED APERTURE FROM SNR CURVE
#------------------------------------------------------------
indx_col = np.where('SNR'==np.array(gctbl.keys()))
x=apertures*pixscale.value
for raw in range(len(gctbl)):
	y = np.array(list(gctbl[raw])[indx_col[0].item():])
	y[np.isnan(y)] = 0.0
	indx_peak = np.where(y==np.max(y))
	if len(y)-1 in indx_peak:
		x_opt=None
	else:
		x_opt=x[indx_peak].item()
		plt.plot(x, y, color='silver', alpha=0.125)
		plt.axvline(x=x_opt, ls='-', linewidth=0.5, color='dodgerblue', alpha=0.125)
		gctbl['APER_OPT'][raw] = x_opt
aper_opt = np.median(gctbl['APER_OPT'])	#	[arcsec]
plt.axvline(x=aper_opt, ls='-', linewidth=2.0, color='tomato', alpha=0.5, label=f'OPT.APERTURE : {round(aper_opt, 3)}\"\n(SEEING*{round(aper_opt/seeing.value, 3)})')
plt.axvline(x=seeing.value, ls='-', linewidth=2.0, color='gold', alpha=0.5, label=f'SEEING : {round(seeing.value, 3)} \"')
plt.title(os.path.basename(inim), fontsize=14)
plt.grid('both', ls='--', color='silver', alpha=0.5)
plt.xlabel('Aperture Diameter [arcsec]', fontsize=14)
plt.ylabel('SNR', fontsize=14)
plt.legend(fontsize=14, framealpha=0.0, loc='upper right')
# plt.yscale('log')
gcoutpng = f'{os.path.splitext(inim)[0]}.gc.png'
# plt.savefig(gcoutpng, dpi=500, overwrite=True)
delt = time.time()-st
print(delt, 'sec')
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
		comment='Best aperture diameter assuming the gaussian profile [pix]',
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
#	Configurations for photometry
#------------------------------------------------------------
prefix = 'gpphot'
path_param = f'{path_config}/{prefix}.param'
path_conv = f'{path_config}/{prefix}.conv'
path_nnw = f'{path_config}/{prefix}.nnw'
path_conf = f'{path_config}/{prefix}.sex'
#	SET APERTURES
#	MAG_AUTO, MAG_APER, MAG_APER_1, MAG_APER_2, MAG_APER_3, MAG_APER_4, MAG_APER_5
apertures = []
for magkey in list(aper_dict.keys()):
	if magkey != 'MAG_AUTO':
		apertures.append(aper_dict[magkey]['size'])
apertures_input = ','.join([str(d) for d in apertures])

outcat = f'{os.path.splitext(inim)[0]}.cat'
#	SE parameters
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
os.system(sexcom(inim, param_insex))
rawtbl = ascii.read(outcat)

n=1
inmagkey = list(aper_dict.keys())[n]

from astropy.wcs import WCS
w = WCS(inim)
imcent = w.all_pix2world(a, b, 1)
raim = imcent[0].item()
decim = imcent[1].item()

c_cent = SkyCoord(raim, decim, unit=u.deg)

from astropy.coordinates import SkyCoord
# c_cent = SkyCoord(hdr['CRVAL1'], hdr['CRVAL2'], unit=u.deg)
c_raw = SkyCoord(rawtbl['ALPHA_J2000'], rawtbl['DELTA_J2000'], unit=u.deg)

sep = c_raw.separation(c_cent)
#	Select 
indx_sel = np.where(
	(rawtbl['CLASS_STAR']>0.9) &
	(rawtbl['FLAGS']<=float(gphot_dict['flagcut'])) &
	(rawtbl[aper_dict[magkey]['errkey']]<=float(gphot_dict['inmagerupper'])) &
	(sep<fov/2*0.9)
	# (sep<(fov/2)*float(gphot_dict['photfraction']))
)
#	QUERY REFERENCE CATALOG
seltbl = rawtbl[indx_sel]
refcatname = gphot_dict['refcatname']
reftbl = querybox(
	refcatname=refcatname,
	racent=c_cent.ra.value,
	decent=c_cent.dec.value,
	outcat=f'{os.path.dirname(inim)}/ref.ngc6946.{refcatname}.cat',
	radius=0.5,
	refmagkey=''
	)
delt = time.time()-st
print(delt, 'sec')
#%%
#	Space distribution
'''
plt.close('all')
plt.figure(figsize=(10,10))
plt.plot(rawtbl['ALPHA_J2000'], rawtbl['DELTA_J2000'], ls='none', marker='o', alpha=0.125, label=f'IMAGE ({len(rawtbl)})')
plt.plot(seltbl['ALPHA_J2000'], seltbl['DELTA_J2000'], ls='none', marker='+', alpha=0.75, label=f'SELECTED ({len(seltbl)})')
plt.plot(reftbl['ra'], reftbl['dec'], ls='none', marker='.', mfc='none', alpha=0.25, label=f'REFERENCE ({len(reftbl)})')

plt.plot(c_cent.ra.value, c_cent.dec.value, ls='none', marker='x', ms=10, c='tomato', mfc='none', alpha=1)
plt.legend(fontsize=14, framealpha=1.0, loc='upper right')
plt.show()'''
#%%
#	Matching catalog
c_ref = SkyCoord(reftbl['ra'], reftbl['dec'], unit=u.deg)
c = SkyCoord(seltbl['ALPHA_J2000'], seltbl['DELTA_J2000'], unit=u.deg)
#	REFERENCE CATALOG FILTERING
#	Magnitude cut and error cut
indx, sep, _ = c.match_to_catalog_sky(c_ref)

mtbl = hstack([seltbl, reftbl[indx]])
mtbl['sep'] = sep
mtbl = mtbl[mtbl['sep']<seeing]
mmtbl = mtbl[
	(
		(mtbl[filte]>float(gphot_dict['refmaglower'])) &
		(mtbl[filte]<float(gphot_dict['refmagupper'])) &
		(mtbl[f'{filte}err']<float(gphot_dict['refmagerupper']))
	)
	]

#	Roop for magnitude kinds
#	e.g. ZP_APER
inzpkey = f"ZP_{aper_dict[inmagkey]['suffix']}"
mmtbl[inzpkey] = mmtbl[filte]-mmtbl[inmagkey] 

# plt.plot(mmtbl[filte], mmtbl[inzpkey], marker='+', ls='none')

#	Valid elements
indx_zp = np.where(
	(~np.isnan(mmtbl[filte].mask)) &
	(~np.isnan(mmtbl[inmagkey].mask))	
)

zptbl = mmtbl[indx_zp]

'''
plt.close()
# sigma=1

for iteration in [1, 2, 3, 4, 5]:
	zplist, zperlist = [], []
	for sigma in np.arange(0.1, 5.0+0.1, 0.1):
		czparr = sigma_clip(
			np.copy(zptbl['ZP_APER']),
			sigma=sigma,
			# maxiters=None,
			# maxiters=3,
			maxiters=iteration,
			# cenfunc=np.median,
			cenfunc=bn.nanmedian,
			copy=False
			)

		zp = np.median(czparr.data[~czparr.mask])
		zper = np.std(czparr.data[~czparr.mask])
		zplist.append(zp)
		zperlist.append(zper)
		print(sigma, round(zp, 3), round(zper, 3), len(czparr[~czparr.mask]))
	plt.plot(np.arange(0.1, 5.0+0.1, 0.1), zperlist, ls='-', label=f'iter:{iteration}')

plt.legend(fontsize=14)
plt.xlabel('SIGMA', fontsize=20)
plt.ylabel('ZPERR', fontsize=20)
plt.title('Sigma clip test', fontsize=14)
plt.grid('both', ls='--', color='silver', alpha=0.5)
plt.show()
'''
#	Sigma-Clipping
sigma=1.5
iteration=5

czparr = sigma_clip(
	np.copy(zptbl[inzpkey]),
	sigma=sigma,
	# maxiters=None,
	# maxiters=3,
	maxiters=iteration,
	# cenfunc=np.median,
	cenfunc=bn.nanmedian,
	copy=False
	)

zp = np.median(czparr.data[~czparr.mask])
zper = np.std(czparr.data[~czparr.mask])
delt = time.time()-st
print(delt, 'sec')
#%%
'''#	ZP
plt.close()
plt.figure(figsize=(9,6))
#
plt.errorbar(zptbl[filte][~czparr.mask], zptbl['ZP_APER'][~czparr.mask], yerr=zptbl[f'{filte}err'][~czparr.mask], mfc='none', marker='s', c='dodgerblue', ls='none', alpha=0.25, label=f'{len(np.where(czparr.mask==False)[0])} ({round(1e2*len(np.where(czparr.mask==False)[0])/len(zptbl), 1)}%)')
plt.plot(zptbl[filte][czparr.mask], zptbl['ZP_APER'][czparr.mask], marker='x', c='tomato', ls='none', alpha=0.5, label=f'{len(np.where(czparr.mask==True)[0])} ({round(1e2*len(np.where(czparr.mask==True)[0])/len(zptbl), 1)}%)')
plt.axhspan(zp-zper, zp+zper, alpha=0.125, color='dodgerblue', label=f'{round(zp, 3)}+/-{round(zper, 3)}')
plt.axhline(y=zp, ls='--', c='b', alpha=0.5)

plt.axvline(x=float(gphot_dict['refmaglower']), ls='--', c='k', alpha=0.5)
plt.axvline(x=float(gphot_dict['refmagupper']), ls='--', c='k', alpha=0.5)
# plt.xlim()
plt.legend(loc='upper center', framealpha=0, fontsize=20)
plt.ylim([zp-0.5, zp+0.5])
plt.xlabel('Ref. mag [AB]', fontsize=20)
plt.ylabel('ZP', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid('both', ls='--', color='silver', alpha=0.5)
plt.tight_layout()
'''
#%%
from astropy.visualization import (ZScaleInterval, MinMaxInterval, SqrtStretch,ImageNormalize)
'''
# Create interval object
interval = ZScaleInterval()
vmin, vmax = interval.get_limits(fits.getdata(inim))

# Create an ImageNormalize object using a SqrtStretch object
norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=SqrtStretch())

plt.close()
plt.figure(figsize=(10,10))
# transform=ax.get_transform('world')
wcs = WCS(inim)
ax = plt.subplot(projection=wcs)
ax.imshow(fits.getdata(inim), origin='lower', cmap='gray', norm=norm)
# plt.scatter(rawtbl['ALPHA_J2000'], rawtbl['DELTA_J2000'], marker='.', color='grey', alpha=0.25)

vmin, vmax = MinMaxInterval().get_limits(zptbl[inzpkey][~czparr.mask])

plt.scatter(zptbl['ALPHA_J2000'][~czparr.mask], zptbl['DELTA_J2000'][~czparr.mask], c=zptbl[inzpkey][~czparr.mask], transform=ax.get_transform('world'), alpha=0.25, vmin=vmin, vmax=vmax)
# plt.scatter(mmtbl['X_IMAGE'], mmtbl['Y_IMAGE'], c=mmtbl['ZP_APER'])
plt.scatter(zptbl['ALPHA_J2000'][czparr.mask], zptbl['DELTA_J2000'][czparr.mask], c='tomato', marker='x', alpha=0.5, transform=ax.get_transform('world'),)
plt.plot(c_cent.ra.value, c_cent.dec.value, ls='none', marker='+', ms=10, c='tomato', mfc='none', alpha=1, transform=ax.get_transform('world'))

plt.colorbar()
# plt.show()
xl, xr = plt.xlim()
plt.xlim([xr, xl])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('RA [deg]', fontsize=20)
plt.ylabel('Dec [deg]', fontsize=20)
plt.tight_layout()'''
# plt.savefig(f'{os.path.splitext(inim)[0]}.star.png', dpi=500, overwrite=True)
delt = time.time()-st
print(delt, 'sec')
# %%
