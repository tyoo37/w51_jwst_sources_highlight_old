from spectral_cube import SpectralCube
from astropy.io import fits
import astropy.units as u
import matplotlib.pyplot as plt

from spectral_cube import SpectralCube
from astropy.io import fits
import astropy.units as u
from astropy.wcs import WCS
import os
from reproject import reproject_interp
import matplotlib.pyplot as plt
from regions import Regions
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import sys
sys.path.append('/home/t.yoo/Paths')
import Paths.Paths as paths

Path = paths.filepaths()
image_filenames ={
    "f140m": "/orange/adamginsburg/jwst/w51/F140M/pipeline/jw06151-o001_t001_nircam_clear-f140m-merged_i2d.fits",
    "f150w": "/orange/adamginsburg/jwst/w51/F150W/pipeline/jw06151-o001_t001_nircam_clear-f150w-merged_i2d.fits",
    "f162m": "/orange/adamginsburg/jwst/w51/F162M/pipeline/jw06151-o001_t001_nircam_clear-f162m-merged_i2d.fits",
    "f182m": "/orange/adamginsburg/jwst/w51/F182M/pipeline/jw06151-o001_t001_nircam_clear-f182m-merged_i2d.fits",
    "f187n": "/orange/adamginsburg/jwst/w51/F187N/pipeline/jw06151-o001_t001_nircam_clear-f187n-merged_i2d.fits",
    "f210m": "/orange/adamginsburg/jwst/w51/F210M/pipeline/jw06151-o001_t001_nircam_clear-f210m-merged_i2d.fits",
    "f335m": "/orange/adamginsburg/jwst/w51/F335M/pipeline/jw06151-o001_t001_nircam_clear-f335m-merged_i2d.fits",
    "f360m": "/orange/adamginsburg/jwst/w51/F360M/pipeline/jw06151-o001_t001_nircam_clear-f360m-merged_i2d.fits",
    "f405n": "/orange/adamginsburg/jwst/w51/F405N/pipeline/jw06151-o001_t001_nircam_clear-f405n-merged_i2d.fits",
    "f410m": "/orange/adamginsburg/jwst/w51/F410M/pipeline/jw06151-o001_t001_nircam_clear-f410m-merged_i2d.fits", # weird, the filename is different from what is downloaded with the STScI pipeline...
    "f480m": "/orange/adamginsburg/jwst/w51/F480M/pipeline/jw06151-o001_t001_nircam_clear-f480m-merged_i2d.fits",
    "f560w": "/orange/adamginsburg/jwst/w51/F560W/pipeline/jw06151-o002_t001_miri_f560w_i2d.fits",
    "f770w": "/orange/adamginsburg/jwst/w51/F770W/pipeline/jw06151-o002_t001_miri_f770w_i2d.fits",
    "f1000w": "/orange/adamginsburg/jwst/w51/F1000W/pipeline/jw06151-o002_t001_miri_f1000w_i2d.fits",
    "f1280w": "/orange/adamginsburg/jwst/w51/F1280W/pipeline/jw06151-o002_t001_miri_f1280w_i2d.fits",
    "f1500w": "/orange/adamginsburg/jwst/w51/F1500W/pipeline/jw06151-o002_t001_miri_f1500w_i2d.fits",
    "f2100w": "/orange/adamginsburg/jwst/w51/F2100W/pipeline/jw06151-o002_t001_miri_f2100w_i2d.fits",
    
}


"""
HCN_moment0_fits = fits.open('/orange/adamginsburg/w51/TaehwaYoo/2017.1.00293.S_W51_B3_LB/moments/w51n_HCN1-0_v=0_red.fits')
HCN_moment0_wcs = WCS(HCN_moment0_fits[0].header, naxis=2)
HCN_moment0_image = HCN_moment0_fits[0].data
reprojected_image, footprint = reproject_interp(
    (HCN_moment0_image, HCN_moment0_wcs),
    alma_b3_wcs,
    shape_out=alma_b3_image.shape
)
hdu = fits.PrimaryHDU(data=reprojected_image, header=alma_b3_wcs.to_header())
hdu.writeto(f'/orange/adamginsburg/jwst/w51/data_cutout_for_knot/HCN_moment0_reprojected_to_alma_b3.fits', overwrite=True)
reprojected_image, footprint = reproject_interp(
    (HCN_moment0_image, HCN_moment0_wcs),
    alma_b6_wcs,
    shape_out=alma_b6_image.shape
)
hdu = fits.PrimaryHDU(data=reprojected_image, header=alma_b6_wcs.to_header())
hdu.writeto(f'/orange/adamginsburg/jwst/w51/data_cutout_for_knot/HCN_moment0_reprojected_to_alma_b6.fits', overwrite=True)



h2k_gtc_image = fits.getdata('/orange/adamginsburg/w51/TaehwaYoo/gtc/adendawson/real_reduction/reduced_images/H2_minus_K.fits')
h2k_gtc_wcs = WCS(fits.getheader('/orange/adamginsburg/w51/TaehwaYoo/gtc/adendawson/real_reduction/reduced_images/H2_minus_K.fits'))

reprojected_image, footprint = reproject_interp(
    (h2k_gtc_image, h2k_gtc_wcs),
    alma_b3_wcs,
    shape_out=alma_b3_image.shape
)
hdu = fits.PrimaryHDU(data=reprojected_image, header=alma_b3_wcs.to_header())
hdu.writeto(f'/orange/adamginsburg/jwst/w51/data_cutout_for_knot/h2k_gtc_reprojected_to_alma_b3.fits', overwrite=True)
reprojected_image, footprint = reproject_interp(
    (h2k_gtc_image, h2k_gtc_wcs),
    alma_b6_wcs,
    shape_out=alma_b6_image.shape
)
hdu = fits.PrimaryHDU(data=reprojected_image, header=alma_b6_wcs.to_header())
hdu.writeto(f'/orange/adamginsburg/jwst/w51/data_cutout_for_knot/h2k_gtc_reprojected_to_alma_b6.fits', overwrite=True)


"""
"""
alma_b3_fits = fits.open(Path.w51n_b3_tt0)
alma_b3_image = alma_b3_fits[0].data
if len(alma_b3_image.shape)!=2:
    alma_b3_image = alma_b3_image[0, 0, :, :]  # Adjust if the data is 3D (e.g., spectral cube)
alma_b3_wcs = WCS(alma_b3_fits[0].header, naxis=2)
print()
alma_b6_fits = fits.open(Path.w51n_b6_tt0)
alma_b6_image = alma_b6_fits[0].data
if len(alma_b6_image.shape)!=2:
    alma_b6_image = alma_b6_image[0, 0, :, :]  # Adjust if the data is 3D (e.g., spectral cube)
alma_b6_wcs = WCS(alma_b6_fits[0].header, naxis=2)

filts = ['f140m', 'f162m', 'f210m', 'f187n', 'f405n']
filts = image_filenames.keys()
for filt in filts:
    if os.path.exists(image_filenames[filt]):
        print(filt, image_filenames[filt])
        # Use naxis=2 to ensure celestial WCS
        input_wcs = WCS(fits.getheader(image_filenames[filt], ext=('SCI', 1)), naxis=2)
        reprojected_image, footprint = reproject_interp(
            (fits.getdata(image_filenames[filt]), input_wcs),
            alma_b3_wcs,
            shape_out=alma_b3_image.shape
        )
        hdu = fits.PrimaryHDU(data=reprojected_image, header=alma_b3_wcs.to_header())
        hdu.writeto(f'/orange/adamginsburg/jwst/w51/reproject_to_alma/{filt}_reprojected_to_alma_w51n_b3.fits', overwrite=True)

        input_wcs = WCS(fits.getheader(image_filenames[filt], ext=('SCI', 1)), naxis=2)
        reprojected_image, footprint = reproject_interp(
            (fits.getdata(image_filenames[filt]), input_wcs),
            alma_b6_wcs,
            shape_out=alma_b6_image.shape
        )
        hdu = fits.PrimaryHDU(data=reprojected_image, header=alma_b6_wcs.to_header())
        hdu.writeto(f'/orange/adamginsburg/jwst/w51/reproject_to_alma/{filt}_reprojected_to_alma_w51n_b6.fits', overwrite=True)
"""

alma_b3_fits = fits.open(Path.w51e_b3_tt0)
alma_b3_image = alma_b3_fits[0].data
if len(alma_b3_image.shape)!=2:
    alma_b3_image = alma_b3_image[0, 0, :, :]  # Adjust if the data is 3D (e.g., spectral cube)
alma_b3_wcs = WCS(alma_b3_fits[0].header, naxis=2)
print()
alma_b6_fits = fits.open(Path.w51e_b6_tt0)
alma_b6_image = alma_b6_fits[0].data
if len(alma_b6_image.shape)!=2:
    alma_b6_image = alma_b6_image[0, 0, :, :]  # Adjust if the data is 3D (e.g., spectral cube)
alma_b6_wcs = WCS(alma_b6_fits[0].header, naxis=2)

filts = ['f140m', 'f162m', 'f210m', 'f187n', 'f405n']
filts = image_filenames.keys()
for filt in filts:
    if os.path.exists(image_filenames[filt]):
        print(filt, image_filenames[filt])
        # Use naxis=2 to ensure celestial WCS
        input_wcs = WCS(fits.getheader(image_filenames[filt], ext=('SCI', 1)), naxis=2)
        reprojected_image, footprint = reproject_interp(
            (fits.getdata(image_filenames[filt]), input_wcs),
            alma_b3_wcs,
            shape_out=alma_b3_image.shape
        )
        hdu = fits.PrimaryHDU(data=reprojected_image, header=alma_b3_wcs.to_header())
        hdu.writeto(f'/orange/adamginsburg/jwst/w51/reproject_to_alma/{filt}_reprojected_to_alma_w51e_b3.fits', overwrite=True)

        input_wcs = WCS(fits.getheader(image_filenames[filt], ext=('SCI', 1)), naxis=2)
        reprojected_image, footprint = reproject_interp(
            (fits.getdata(image_filenames[filt]), input_wcs),
            alma_b6_wcs,
            shape_out=alma_b6_image.shape
        )
        hdu = fits.PrimaryHDU(data=reprojected_image, header=alma_b6_wcs.to_header())
        hdu.writeto(f'/orange/adamginsburg/jwst/w51/reproject_to_alma/{filt}_reprojected_to_alma_w51e_b6.fits', overwrite=True)
"""

for filt in filts:
    if os.path.exists(image_filenames[filt]):
        print(filt, image_filenames[filt])
        # Use naxis=2 to ensure celestial WCS
        input_wcs = WCS(fits.getheader(image_filenames[filt], ext=('SCI', 1)), naxis=2)
        reprojected_image, footprint = reproject_interp(
            (fits.getdata(image_filenames[filt]), input_wcs),
            alma_b3_wcs,
            shape_out=alma_b3_image.shape
        )
        hdu = fits.PrimaryHDU(data=reprojected_image, header=alma_b3_wcs.to_header())
        hdu.writeto(f'/orange/adamginsburg/jwst/w51/data_cutout_for_knot/{filt}_reprojected_to_alma_b3.fits', overwrite=True)

        input_wcs = WCS(fits.getheader(image_filenames[filt], ext=('SCI', 1)), naxis=2)
        reprojected_image, footprint = reproject_interp(
            (fits.getdata(image_filenames[filt]), input_wcs),
            alma_b6_wcs,
            shape_out=alma_b6_image.shape
        )
        hdu = fits.PrimaryHDU(data=reprojected_image, header=alma_b6_wcs.to_header())
        hdu.writeto(f'/orange/adamginsburg/jwst/w51/data_cutout_for_knot/{filt}_reprojected_to_alma_b6.fits', overwrite=True)
"""
"""

CO_moment0_fits = fits.open('/orange/adamginsburg/w51/TaehwaYoo/fitsfile_jwst/moment0.fits')
CO_moment0_wcs = WCS(CO_moment0_fits[0].header, naxis=2)
CO_moment0_image = CO_moment0_fits[0].data

vla_kuband_fits = fits.open('/orange/adamginsburg/w51/TaehwaYoo/vla/2016paper/W51Ku_C_Aarray_continuum_2048_high_uniform.clean.image.fits')
vla_kuband_wcs = WCS(vla_kuband_fits[0].header, naxis=2)
vla_kuband_image = vla_kuband_fits[0].data

reprojected_image, footprint = reproject_interp(
    (CO_moment0_image, CO_moment0_wcs),
    alma_b3_wcs,
    shape_out=alma_b3_image.shape
)
hdu = fits.PrimaryHDU(data=reprojected_image, header=alma_b3_wcs.to_header())
hdu.writeto(f'/orange/adamginsburg/jwst/w51/data_cutout_for_knot/CO_moment0_reprojected_to_alma_b3.fits', overwrite=True)

reprojected_image, footprint = reproject_interp(
    (vla_kuband_image, vla_kuband_wcs),
    alma_b3_wcs,
    shape_out=alma_b3_image.shape
)
hdu = fits.PrimaryHDU(data=reprojected_image, header=alma_b3_wcs.to_header())
hdu.writeto(f'/orange/adamginsburg/jwst/w51/data_cutout_for_knot/vla_kuband_reprojected_to_alma_b3.fits', overwrite=True)

reprojected_image, footprint = reproject_interp(
    (CO_moment0_image, CO_moment0_wcs),
    alma_b6_wcs,
    shape_out=alma_b6_image.shape
)
hdu = fits.PrimaryHDU(data=reprojected_image, header=alma_b6_wcs.to_header())
hdu.writeto(f'/orange/adamginsburg/jwst/w51/data_cutout_for_knot/CO_moment0_reprojected_to_alma_b6.fits', overwrite=True)   
reprojected_image, footprint = reproject_interp(
    (vla_kuband_image, vla_kuband_wcs),
    alma_b6_wcs,
    shape_out=alma_b6_image.shape
)
hdu = fits.PrimaryHDU(data=reprojected_image, header=alma_b6_wcs.to_header())
hdu.writeto(f'/orange/adamginsburg/jwst/w51/data_cutout_for_knot/vla_kuband_reprojected_to_alma_b6.fits', overwrite=True)
"""