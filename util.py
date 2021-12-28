import os
from astropy.time import Time
from astropy.io import fits
import numpy as np
import astroalign as aa
from astropy.nddata import CCDData
#------------------------------------------------------------
def scaling_func(arr): return 1/np.ma.median(arr)
#------------------------------------------------------------
def bn_median(masked_array, axis=None):
    """
    Perform fast median on masked array

    Parameters

    masked_array : `numpy.ma.masked_array`
        Array of which to find the median.

    axis : int, optional
        Axis along which to perform the median. Default is to find the median of
        the flattened array.
    """
    import numpy as np
    import bottleneck as bn
    data = masked_array.filled(fill_value=np.NaN)
    med = bn.nanmedian(data, axis=axis)
    # construct a masked array result, setting the mask from any NaN entries
    return np.ma.array(med, mask=np.isnan(med))
#------------------------------------------------------------
def align_astroalign(srcim, tgtim,):
	'''
	http://quatrope.github.io/astroalign/
	'''
	print(f'Align {os.path.basename(srcim)} to {os.path.basename(tgtim)}')
	tdata, thdr = fits.getdata(tgtim, header=True)
	tzero = np.median(tdata[~np.isnan(tdata)].flatten())

	#	Registered image
	rdata, footprint = aa.register(
		fits.getdata(srcim),
		fits.getdata(tgtim),
		fill_value=np.NaN,
		)
	rzero = np.median(rdata[~np.isnan(rdata)].flatten())
	zero_offset = tzero-rzero
	print(f'\t--> (Aligned image) = {os.path.basename(srcim)} - ({round(zero_offset, 3)})')
	rdata = rdata+zero_offset
	return CCDData(rdata, unit='adu', header=fits.getheader(srcim))
#------------------------------------------------------------
from astropy.io import fits
from astropy.nddata import CCDData
import ccdproc

def imcombine_ccddata(aligned_imlist, imlist=None):
	"""
	"""
	comment = f"""{'-'*60}\n#\tCOMBINE {len(aligned_imlist)} OBJECT FRAMES\n{'-'*60}"""
	print(comment)
	#	Info. for combined image
	if imlist==None:
		imlist = [f"{indata.header['PATHPRC']}/{indata.header['IMNAME']}" for indata in aligned_imlist]
	jd = calc_jd_median(imlist)
	comim = combine_name(imlist)
	#	Print
	for n, inim in enumerate(imlist):
		print(f'[{n}] {os.path.basename(inim)}')
	print(f'==> {os.path.basename(comim)}')

	#	Combine
	combiner = ccdproc.Combiner(aligned_imlist, dtype=np.float32)
	comdata = combiner.median_combine(median_func=bn_median)	

	#	Header
	comdata.header = aligned_imlist[0].header
	comdata.header['IMNAME'] = os.path.basename(comim)
	comdata.header['NCOMBINE'] = len(aligned_imlist)
	for n, inim in enumerate(imlist):
		comdata.header[f'COMB{n}'] = os.path.basename(inim)
	comdata.header['DATE-OBS'] = (jd.isot, 'YYYY-MM-DDThh:mm:ss observation start, UT')
	comdata.header['JD'] = (jd.value, 'Julian Date at start of exposure')
	comdata.header['MJD'] = (jd.mjd, 'Modified Julian Date at start of exposure')
	comdata.header['EXPTIME'] = (calc_total_exptime(imlist), 'Total Exposure time [sec]')
	#	Save
	# comdata.write(comim, overwrite=True)
	fits.writeto(comim, comdata.data, header=comdata.meta, overwrite=True)
	# return comdata
#------------------------------------------------------------
#	Alignment
def gregistering(imlist, refim):
	import alipy
	identifications = alipy.ident.run(
		refim,
		imlist,
		visu=False
		)
	# for inim, id in zip(imlist, identifications):
	for id in identifications:
		if id.ok == True: # i.e., if it worked
			# print("%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio))
			print(f"{id.ukn.name} : {id.trans}, flux ratio {round(id.medfluxratio, 2)}" % (id.ukn.name, id.trans, id.medfluxratio))
			params_align = dict(
				filepath = id.ukn.filepath,
				uknstarlist = id.uknmatchstars,
				refstarlist = id.refmatchstars,
				shape = alipy.align.shape(refim),
				outdir = os.path.dirname(refim),
				makepng = False
				)
			alipy.align.irafalign(**params_align)
		else:
			print(f"{id.ukn.name} : No transformation found!")
	outputshape = alipy.align.shape(refim)
#------------------------------------------------------------
def imcombine_iraf(imlist):
	'''
	path_data = '/data3/paek/factory/doao/20201209-1m-IMSNG'
	images_to_align = sorted(glob.glob('/data3/paek/factory/doao/20201209-1m-IMSNG/Calib-DOAO*-R-60.fits'))
	ref_image = '/data3/paek/factory/doao/20201209-1m-IMSNG/Calib-DOAO-NGC6946-20201209-094720-R-60.fits'
	'''
	from pyraf import iraf
	import glob, os
	from astropy.io import fits
	from astropy.time import Time
	import numpy as np

	path_data = os.path.dirname(imlist[0])
	make_list_file(
		imlist,
		outname=f'{path_data}/imcombine.tmp',
		)
	jd = calc_jd_median(imlist)

	comim = combine_name(imlist)

	print(f'#\t{len(imlist)} IMAGE IMCOMBINE --> {os.path.basename(comim)}')
	param_imcomb = dict(
						input=f"@{path_data}/imcombine.tmp",
						output=comim,
						combine="median",
						project="no",
						reject="none",
						scale="none",
						zero="mode",
						)
	iraf.imcombine(**param_imcomb)

	for i, inim in enumerate(imlist):
		fits.setval(comim, keyword=f'IMCOMB{i}', value=os.path.basename(inim), comment=f'Combined image {i}')
	fits.setval(comim, keyword='DATE-OBS', value=jd.isot, comment='YYYY-MM-DDThh:mm:ss observation start, UT')
	fits.setval(comim, keyword='JD', value=jd.value, comment='Julian Date at start of exposure')
	fits.setval(comim, keyword='MJD', value=jd.mjd, comment='Modified Julian Date at start of exposure')
	fits.setval(comim, keyword='EXPTIME', value=exptime, comment='Total Exposure time [sec]')

	return outim
#------------------------------------------------------------
def calc_jd_median(imlist):
	return Time(np.median([fits.getheader(inim)['JD'] for inim in imlist]), format='jd')
#------------------------------------------------------------
def calc_total_exptime(imlist):
	return int(np.sum([fits.getheader(inim)['EXPTIME'] for inim in imlist]))
#------------------------------------------------------------
def make_list_file(imlist, outname,):
	f = open(f'{outname}', 'w')
	for i, inim in enumerate(imlist):
		f.write(f'{inim}\n')
	f.close()
#------------------------------------------------------------
def split_dateobs(dateobs):
	utdate = dateobs.split('T')[0].replace('-', '')
	uttime = dateobs.split('T')[1].replace(':', '')[:6]
	return utdate, uttime
#------------------------------------------------------------
def combine_name(imlist):
	"""
	"""
	#	The first image
	hdr = fits.getheader(imlist[0])
	path_data = os.path.dirname(imlist[0])
	#	Component
	jdmed = calc_jd_median(imlist)
	dateobsmed = jdmed.isot
	utdate, uttime = split_dateobs(dateobsmed)
	exptime = calc_total_exptime(imlist)
	filte = hdr['FILTER']
	obj = hdr['OBJECT']
	obs = hdr['OBSERVAT']

	outim = f'{path_data}/Calib-{obs}-{obj}-{utdate}-{uttime}-{filte}-{exptime}.com.fits'
	return outim
#------------------------------------------------------------
def subtraction_routine(inim, refim):
	'''
	obs = 'LOAO'
	path_refim = '/data3/paek/factory/ref_frames/{}'.format(obs)
	inim = '/data3/paek/factory/test/Calib-LOAO-NGC6946-20201213-014607-R-180-com.fits'

	obj = 'NGC6946'
	filte = 'R'
	'''
	# inseeing = fits.getheader(inim)['seeing']
	# refseeing = fits.getheader(refim)['seeing']

	# if inseeing > refseeing:
	# 	images_to_align = [inim]
	# 	ref_image = refim
	# else:
	# 	images_to_align = [refim]
	# 	ref_image = inim
	gregistering([refim], inim)
	#	Registered reference image
	grefim = '{}/{}'.format(os.path.dirname(inim), os.path.basename(outim_gregistering(refim)))
	subim = hotpants(inim, grefim, iu=60000, tu=6000000000, tl=-100000)
	ds9com = 'ds9 {} {} {}&'.format(inim, grefim, subim)
	# os.system(ds9com)
	return subim, ds9com
#------------------------------------------------------------
def scale_flux_zp(zp, zp_ref):
	fref_f = 10.**(-0.4*(zp-zp_ref))
	f_fref = 1./fref_f
	return f_fref