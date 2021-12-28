#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from astropy.io import ascii
from astropy import units as u

def getccdinfo(obs, path_ccd):
	'''GET CCD INFORMATION (GAIN, PIXEL SCALE, FOV)
	
	ccddict = getccdinfo(obs, path_ccd)
	
	INPUT:
	path_ccd = '/home/sonic/Research/table/obs.dat'

	OUTPUT:
	Type : Dictionary
	obs, gain, pixelscale, fov, rdnoise
	'''
	ccdtbl_ = ascii.read(path_ccd)
	ccdtbl = ccdtbl_[ccdtbl_['obs']==obs]

	outdict = dict()
	outdict['obs'] = obs
	# outdict['gain'] = ccdtbl['gain'].item() * u.electron / u.second
	outdict['gain'] = ccdtbl['gain'].item() * u.electron / u.adu
	outdict['pixelscale'] = ccdtbl['pixelscale'].item() *  u.arcsecond / u.pixel
	outdict['fov'] = ccdtbl['fov'].item() * u.arcmin
	outdict['rdnoise'] = ccdtbl['RDnoise'].item() * u.electron

	return outdict

def file2dict(path_infile):
	out_dict = dict()
	f = open(path_infile)
	for line in f:
		key, val = line.split()
		out_dict[key] = val
	return out_dict