#%%
#	Path
path_table = '/home/paek/imsngpy/table'
path_config = '/home/paek/imsngpy/config'

prefix = 'gpphot'
path_param = f'{path_config}/{prefix}.param'
path_conv = f'{path_config}/{prefix}.conv'
path_nnw = f'{path_config}/{prefix}.nnw'
path_conf = f'{path_config}/{prefix}.sex'

#	Table
from astropy.io import ascii
ccdtbl = ascii.read(f'{path_table}/ccd.tsv') 

def file2dict(path_infile):
	out_dict = dict()
	f = open(path_infile)
	for line in f:
		key, val = line.split()
		out_dict[key] = val
	return out_dict

path_gphot = f'{path_config}/gphot.config'

import os
if os.path.exists(path_gphot) == True:
	gphot_dict = file2dict(path_gphot)
else:
	print('[NOTICE] There is no gphot.config. Use default configuration.')
	gphot_dict = {
		'imkey': 'Calib*com.fits',
		'photfraction': '0.75',
		'refcatname': 'PS1',
		'refqueryradius': '1',
		'refmaglower': '12',
		'refmagupper': '20',
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
inim = '/data3/paek/factory/test/phot/Calib-LOAO-M99-20210421-063118-R-180.com.fits'
from astropy.io import fits
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
import sys
sys.path.append('/home/paek/imsngpy')
from misc import *


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

aperture_dict = dict(
	MAG_AUTO=dict(
		errkey='MAGERR_AUTO',
		aperture=0.,
		comment='',
	),
	#	BEST SNR ASSUMING GAUSSIAN PROFILE
	MAG_APER_1=dict(
		errkey='MAGERR_APER_1',
		aperture=2*0.6731*peeing.value,
		comment='',
	),
)

# apertures = peeing.value*np.arange(0.25, 8+0.25, 0.25)
apertures = peeing.value*np.linspace(0.25, 5, 32)
apertures_input = ','.join([str(d) for d in apertures])
#%%
prefix = 'growthcurve'
path_param = f'{path_config}/{prefix}.param'
path_conv = f'{path_config}/{prefix}.conv'
path_nnw = f'{path_config}/{prefix}.nnw'
path_conf = f'{path_config}/{prefix}.sex'

outcat = f'{os.path.splitext(inim)[0]}.gc.cat'
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
#	
a = hdr['naxis1']/2.
b = hdr['naxis2']/2.
#	Small squre based on frac
frac = 0.95
a_ = a*np.sqrt(frac)
b_ = b*np.sqrt(frac)
gctbl = rawtbl[
	(rawtbl['FLAGS']==0) &
	(rawtbl['CLASS_STAR']>0.9) &
	(rawtbl['X_IMAGE']>a-a_) & (rawtbl['X_IMAGE']<a+a_) &
	(rawtbl['Y_IMAGE']>b-b_) & (rawtbl['Y_IMAGE']<b+b_)
]
gctbl['APER_OPT'] = 0.0
for n in range(len(apertures)):
	if n==0:
		gctbl[f'SNR'] = gctbl[f'FLUX_APER']/gctbl[f'FLUXERR_APER']
	else:
		gctbl[f'SNR_{n}'] = gctbl[f'FLUX_APER_{n}']/gctbl[f'FLUXERR_APER_{n}']
#%%
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
plt.savefig(gcoutpng, dpi=500, overwrite=True)
#%%

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
		suffix='AUTO',
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

#%%
prefix = 'gpphot'
path_param = f'{path_config}/{prefix}.param'
path_conv = f'{path_config}/{prefix}.conv'
path_nnw = f'{path_config}/{prefix}.nnw'
path_conf = f'{path_config}/{prefix}.sex'

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
magkey = list(aper_dict.keys())[n]


c_cent = SkyCoord(hdr['CRVAL1'], hdr['CRVAL2'], unit=u.deg)
c_raw = SkyCoord(rawtbl['ALPHA_J2000'], rawtbl['DELTA_J2000'], unit=u.deg)

sep = c_raw.separation(c_cent)



indx_sel = np.where(
	(rawtbl['CLASS_STAR']>0.9) &
	(rawtbl['FLAGS']<=float(gphot_dict['flagcut'])) &
	(rawtbl[aper_dict[magkey]['errkey']]<=float(gphot_dict['inmagerupper'])) &
	(sep<fov/2*float(gphot_dict['photfraction']))
)

seltbl = rawtbl[indx_sel]

from query import *
refcatname = 'PS1'
reftbl = querybox(
	refcatname=refcatname,
	racent=184.683,
	decent=14.401,
	outcat=f'{os.path.dirname(inim)}/ref.test.{refcatname}.cat',
	radius=0.5,
	refmagkey=''
	)

from astropy.coordinates import SkyCoord
c_ref = SkyCoord(reftbl['ra'], reftbl['dec'], unit=u.deg)
c = SkyCoord(seltbl['ALPHA_J2000'], seltbl['DELTA_J2000'], unit=u.deg)

indx, sep, _ = c.match_to_catalog_sky(c_ref)
from astropy.table import hstack
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

inmagkey = 'MAG_APER'

zparr = mmtbl[filte]-mmtbl[inmagkey]
plt.plot(mmtbl[filte], zparr, marker='+', ls='none')


indx_zp = np.where(
	(~np.isnan(mmtbl[filte].mask)) &
	(~np.isnan(mmtbl[inmagkey].mask))	
)
zparr = mmtbl[filte]-mmtbl[inmagkey]
zparrm = zparr[indx_zp]


zparrmc = sigma_clip(np.copy(zparrm))

zp = np.median(zparrmc.data[~zparrmc.mask])
zper = np.std(zparrmc.data[~zparrmc.mask])

print(zp, zper)

plt.scatter(rawtbl['X_IMAGE'], rawtbl['Y_IMAGE'], marker='.', color='grey', alpha=0.25)
plt.scatter(mmtbl['X_IMAGE'][indx_zp], mmtbl['Y_IMAGE'][indx_zp], c=zparrmc)
plt.colorbar()

'''
from astropy.stats import sigma_clip
#	REMOVE BLANK ROW (=99)	
indx_avail      = np.where( (mmtbl[filte] != 99) & (mmtbl[filte] != 99) )
intbl           = intbl[indx_avail]
zplist          = np.copy(intbl[refmagkey] - intbl[filte])
intbl['zp']		= zplist
#	SIGMA CLIPPING
zplist_clip     = sigma_clip(zplist, sigma=sigma, maxiters=None, cenfunc=np.median, copy=False)
indx_alive      = np.where( zplist_clip.mask == False )
indx_exile      = np.where( zplist_clip.mask == True )
#	RE-DEF. ZP LIST AND INDEXING CLIPPED & NON-CLIPPED
intbl_alive     = intbl[indx_alive]
intbl_exile     = intbl[indx_exile]
#	ZP & ZP ERR. CALC.
if method == 'default':
	zp              = np.median(np.copy(intbl_alive['zp']))
	zper			= np.std(np.copy(intbl_alive['zp']))
elif method == 'weightedmean':
	print(method)
	mager = sqsum(intbl_alive[inmagerkey], intbl_alive[refmagerkey])
	w0 = 1/mager
	w = w0/np.sum(w0)
	zp = np.sum(w*intbl_alive['zp'])/np.sum(w)
	zper = 1/np.sqrt(np.sum(w))'''