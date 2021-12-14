
#	Path
path_table = '/home/paek/imsngpy/table'
path_config = '/home/paek/imsngpy/config'

path_param = f'{path_config}/gpphot.param'
path_conv = f'{path_config}/gpphot.conv'
path_nnw = f'{path_config}/gpphot.nnw'
path_conf = f'{path_config}/gpphot.sex'

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
photdict = file2dict(path_gphot)

#	Test image
inim = '/data3/paek/factory/test/phot/Calib-LOAO-M99-20210421-063118-R-180.com.fits'
from astropy.io import fits
hdr = fits.getheader(inim)
if ('OBSERVAT' in hdr.keys()) & ('CCDNAME' in hdr.keys()):
	obs = hdr['OBSERVAT']
	ccd = hdr['CCDNAME']
else:
	#	observatory from filename
	obs = os.path.basename(inim).split('-')[1]
	ccd = hdr['CCDNAME']
print(obs)
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

prefix = 'simple'
path_conf = f'{path_config}/{prefix}.sex'
path_param = f'{path_config}/{prefix}.param'
path_nnw = f'{path_config}/{prefix}.nnw'
path_conv = f'{path_config}/{prefix}.conv'

#	Single
if ('SEEING' in hdr.keys()) & ('PEEING' in hdr.keys()):
	seeing = hdr['SEEING']*u.arcsec
	peeing = hdr['PEEING']*u.pix
else:
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
					GAIN = str(gain.value),
					PIXEL_SCALE = str(pixscale.value),
					#------------------------------
					#	STAR/GALAXY SEPARATION
					#------------------------------
					SEEING_FWHM = str(seeing.value),
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
		print('DONE')'''