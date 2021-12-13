import os


def scaling_func(arr): return 1/np.ma.median(arr)

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
from astropy.io import fits
from astropy.nddata import CCDData
import ccdproc

def imcombine_ccddata(ccddatalist):
	"""
	"""
	comment = f"""{'-'*60}\n#\t{filte} FLAT MASTER FRAME (<--{len(imlist)} frames)\n{'-'*60}"""
	print(comment)

	#	Combine
	combiner = ccdproc.Combiner(ccddatalist)
	# combiner.minmax_clipping()
	# combiner.scaling = scaling_func
	com = combiner.median_combine(median_func=bn_median)	
	#	Header
	nmflat.header = fits.getheader(imlist[0])
	nmflat.header['NCOMBINE'] = len(imlist)
	nmflat.header['SUBBIAS'] = mbias.header['FILENAME']
	nmflat.header['SUBDARK'] = mdark.header['FILENAME']
	for n, inim in enumerate(imlist): nmflat.header[f'COMB{n}'] = inim
	#	Save
	dateobs = nmflat.header['DATE-OBS'].split('T')[0].replace('-', '')
	filename = f'{dateobs}-n{filte}.fits'
	nmflat.header['FILENAME'] = filename
	nmflat.write(f'{os.path.dirname(inim)}/{filename}', overwrite=True)
	return nmflat
#------------------------------------------------------------
#	Alignment
def gregistering(images_to_align, refim):
	import alipy

	identifications = alipy.ident.run(
		refim,
		images_to_align,
		visu=False
		)
	# for inim, id in zip(images_to_align, identifications):
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
			# print("%20s : no transformation found !" % (id.ukn.name))
			print(f"{id.ukn.name} : No transformation found!")
	outputshape = alipy.align.shape(refim)
#------------------------------------------------------------
def imcombine_routine(imlist):
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

	for i, inim in imlist: fits.setval(inim, keyword=f'IMCOMB{i}', value=os.path.basename(inim), comment=f'Combined image {i}')
	fits.setval(inim, keyword='DATE-OBS', value=jd.isot, comment='YYYY-MM-DDThh:mm:ss observation start, UT')
	fits.setval(inim, keyword='JD', value=jd.value, comment='Julian Date at start of exposure')
	fits.setval(inim, keyword='MJD', value=jd.mjd, comment='Modified Julian Date at start of exposure')
	fits.setval(inim, keyword='EXPTIME', value=exptime, comment='Total Exposure time [sec]')

	return outim

def calc_jd_median(imlist):
	return 	Time(np.median([fits.getheader(inim)['JD'] for inim in imlist]), format='jd')

def make_list_file(imlist, outname,):
	f = open(f'{outname}', 'w')
	for i, inim in enumerate(imlist):
		f.write(f'{inim}\n')
	f.close()

def split_dateobs(dateobs):
	utdate = dateobs.split('T')[0].replace('-', '')
	uttime = dateobs.split('T')[1].replace(':', '')[:6]
	return utdate, uttime

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
	exptime = int(np.sum([fits.getheader(inim)['EXPTIME'] for inim in imlist]))
	filte = hdr['FILTER']
	obj = hdr['OBJECT']

	outim = f'{path_data}/Calib-{obs}-{obj}-{utdate}-{uttime}-{filte}-{exptime}.com.fits'
	return outim