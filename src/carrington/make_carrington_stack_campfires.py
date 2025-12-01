import glob
from rectify import rectify
from astropy.io import fits
from astropy.io.fits.hdu.compressed import CompImageHDU
import astropy.constants
import astropy.time
import os

# Select instrument
#instrument = r'HRIEUV'
#instrument = r'AIA'

# Change input data path as required
#path = r'/home/francisv/CLOSE-UP/new_base_code_20250911/test_run_20200530/data/hrieuv/'
#path = r'/home/francisv/CLOSE-UP/new_base_code_20250911/test_run_20200530/data/aia/171'

# Change output (Carrington) folder as required
#outbase = r'/home/francisv/CLOSE-UP/new_base_code_20250911/test_run_20200530/carrington/hrieuv/'
#outbase = r'/home/francisv/CLOSE-UP/new_base_code_20250911/test_run_20200530/carrington/aia/171'

# Check lonlims and latlims in jHelioviewer, the shape remains mostly the same 
#shape = (2400, 2400)
#lonlims = (248, 289)
#latlims = (-12, 26)


def make_carrington_stack_campfires(instrument, path, outbase, shape, lonlims, latlims):
    
    if instrument == 'HRIEUV':
        files = glob.glob(os.path.join(path, 'solo_L2_eui-hrieuv*-image*.fits'))
    else:
        files = glob.glob(os.path.join(path, 'aia.lev1.171A*lev1.fits'))  # TEMPORARY SOLUTION!
       #files = glob.glob(os.path.join(path, 'aia.lev1.171A*lev1p5.fits'))  # We want to take the level 1.5 files as input but we need a CROTA fix in rectify for that!

    first_file = files[0]
   
    with fits.open(first_file) as hdul:
        if (any(isinstance(hdu, CompImageHDU) for hdu in hdul)):
            first_header = fits.getheader(first_file, 1)  # for compressed fits files
        else:
            first_header = fits.getheader(first_file)  # for uncompressed fits files

    reference_date = first_header['DATE-OBS']
    rate_wave = '171'  # or None for default Carrington rate 
    
    solar_r = 1.004

    for i, f in enumerate(files):

        instrumentmap = fits.getdata(f)
      
        with fits.open(f) as hdul:
            if (any(isinstance(hdu, CompImageHDU) for hdu in hdul)):
                header = fits.getheader(f, 1)  # for compressed fits files
            else:
                header = fits.getheader(f)  # for uncompressed fits files

        
       # the following line throws an error because there are no CROTA or CROTA2 in the header of the AIA Level 1.5 files. They do contain PCi_j, and the Level 1 files contain CROTA2
        spherical = rectify.CarringtonTransform(header, radius_correction = solar_r, reference_date = reference_date, rate_wave = rate_wave)  # with differential rotation
       #spherical = rectify.CarringtonTransform(header, radius_correction=solar_r)  # without differential rotation
        spherizer = rectify.Rectifier(spherical)
        carrington = spherizer(instrumentmap, shape, lonlims, latlims, opencv=False, order=2)

        out = os.path.join(outbase, os.path.basename(f))
    
        if instrument == 'HRIEUV':
            out = out.replace('_L2_', '_L3_')
        else:
            out = out.replace('_lev1', '_lev3')  # TEMPORARY SOLUTION!
           #out = out.replace('_lev1p5', '_lev3')  # switch to this line of code when the CROTA issue is fixed in rectify
            
        out = out.replace('.fits', '_carrington.fits')
        
        header['MAPPINGR'] = solar_r*astropy.constants.R_sun.value
        header['CACRPIX1'] = (shape[0] + 1)/2
        header['CACRPIX2'] = (shape[1] + 1)/2
        header['CACRVAL1'] = (lonlims[1] + lonlims[0])/2
        header['CACRVAL2'] = (latlims[1] + latlims[0])/2
        header['CACDELT1'] = (lonlims[1] - lonlims[0])/(shape[0]-1)
        header['CACDELT2'] = (latlims[1] - latlims[0])/(shape[1]-1)

       #header['WCSNAME'] = 'Carrington heliographic'
       #header['CTYPE1']  = 'CRLN'
       #header['CTYPE2']  = 'CRLT'
       #header['CUNIT1']  = 'deg     '
       #header['CRPIX1']  = (shape[0] + 1)/2
       #header['CRPIX2']  = (shape[1] + 1)/2
       #header['CRVAL1']  = (lonlims[1] + lonlims[0])/2
       #header['CRVAL2']  = (latlims[1] + latlims[0])/2
       #header['CDELT1']  = (lonlims[1] - lonlims[0])/(shape[0]-1)
       #header['CDELT2']  = (latlims[1] - latlims[0])/(shape[1]-1)

        if ("OSCNMEAN" in header):
            del header["OSCNMEAN"]
        if ("OSCNRMS" in header):
            del header["OSCNRMS"]

        fits.writeto(out, carrington, header=header, overwrite=True)