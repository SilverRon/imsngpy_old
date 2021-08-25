from astropy.io import ascii
from astropy import units as u

def getccdinfo(obs, path_ccd):
	'''GET CCD INFORMATION (GAIN, PIXEL SCALE, FOV)
	
	ccddict = getccdinfo(obs, path_ccd)
	
	INPUT:
	path_ccd = '/home/sonic/Research/table/obs.txt'

	OUTPUT:
	Type : Dictionary
	obs, gain, pixelscale, fov, rdnoise
	'''
	ccdtbl_ = ascii.read(path_ccd)
	ccdtbl = ccdtbl_[ccdtbl_['obs']==obs]

	outdict = dict()
	outdict['obs'] = obs
	outdict['gain'] = ccdtbl['gain'].item() * u.electron / u.second
	outdict['pixelscale'] = ccdtbl['pixelscale'].item() *  u.arcsecond / u.pixel
	outdict['fov'] = ccdtbl['fov'].item() * u.deg
	outdict['rdnoise'] = ccdtbl['RDnoise'].item()

	return outdict