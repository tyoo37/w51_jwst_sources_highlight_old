from reproject import reproject_interp
from astropy.io import fits
from astropy.wcs import WCS
# Make a new WCS covering B's footprint, using A's pixel scale
from astropy.wcs.utils import skycoord_to_pixel
from astropy.coordinates import SkyCoord
from astropy.nddata.utils import NoOverlapError
from astropy.wcs.utils import proj_plane_pixel_scales

from astropy.wcs.utils import fit_wcs_from_points
from astropy import units as u
import numpy as np


import sys
sys.path.append('/home/t.yoo/Paths')
import Paths.Paths as paths
image_filenames ={
    "f140m": "/orange/adamginsburg/jwst/w51/F140M/pipeline/jw06151-o001_t001_nircam_clear-f140m-merged_i2d.fits",
    "f150w": "/orange/adamginsburg/jwst/w51/F150W/pipeline/jw06151-o001_t001_nircam_clear-f150w-merged_i2d.fits",
    "f162m": "/orange/adamginsburg/jwst/w51/F162M/pipeline/jw06151-o001_t001_nircam_clear-f162m-merged_i2d.fits",
    "f182m": "/orange/adamginsburg/jwst/w51/F182M/pipeline/jw06151-o001_t001_nircam_clear-f182m-merged_i2d.fits",
    "f187n": "/orange/adamginsburg/jwst/w51/F187N/pipeline/jw06151-o001_t001_nircam_clear-f187n-merged_i2d.fits",
    "f210m": "/orange/adamginsburg/jwst/w51/F210M/pipeline/jw06151-o001_t001_nircam_clear-f210m-merged_i2d.fits",
    "f335m": "/orange/adamginsburg/jwst/w51/F300M/pipeline/jw06151-o001_t001_nircam_clear-f335m-merged_i2d.fits",
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

image_sub_filenames = {
    "f405n-f410m": "/orange/adamginsburg/jwst/w51/F405_minus_F410cont_pipeline_v0.1.fits",
    "f410m-f405n": "/orange/adamginsburg/jwst/w51/F410_minus_F405_fractional_bandwidth_pipeline_v0.1.fits",
    "f187n-f182m": "/orange/adamginsburg/jwst/w51/F187_minus_F182cont_pipeline_v0.1.fits",

}
Path = paths.filepaths()
alma_b3_fits = fits.open(Path.w51n_b3_tt0)

alma_b3_image = alma_b3_fits[0].data.squeeze()
  # Adjust if the data is 3D (e.g., spectral cube)

header_b3 = alma_b3_fits[0].header.copy()
alma_b3_wcs = WCS(header_b3, naxis=2)

#alma_b3_wcs = WCS(alma_b3_fits[0].header, naxis=2)

alma_b6_fits = fits.open(Path.w51n_b6_tt0)

alma_b6_image = alma_b6_fits[0].data.squeeze()  # Adjust if the data is 3D (e.g., spectral cube)
header_b6 = alma_b6_fits[0].header.copy()
from reproject import reproject_interp
from astropy.io import fits
# Make a new WCS covering B's footprint, using A's pixel scale
from astropy.wcs.utils import skycoord_to_pixel
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.wcs.utils import WCS

for header in [header_b3, header_b6]:
    keys_to_remove = []
    for key in header.keys():
        # Remove any keyword ending with 3-9 (e.g., CTYPE3, CRVAL4, NAXIS5, etc.)
        if key.endswith(('3','4','5','6','7','8','9')):
            keys_to_remove.append(key)
        # Remove PCi_j where i or j > 2 (e.g., PC3_1, PC1_3, PC4_4, etc.)
        if key.startswith('PC'):
            parts = key.split('_')
            if len(parts) == 2 and (parts[0][2:] not in ['1','2'] or parts[1] not in ['1','2']):
                keys_to_remove.append(key)
        # Remove WCSAXES, CNAME3, CROTA3, etc.
        if key.startswith(('WCSAXES', 'CNAME', 'CROTA')) and (key[-1] in '3456789'):
            keys_to_remove.append(key)
    for key in set(keys_to_remove):
        if key in header:
            del header[key]
header_b3['CTYPE1'] = 'RA---TAN'
header_b3['CTYPE2'] = 'DEC--TAN'
header_b6['CTYPE1'] = 'RA---TAN'
header_b6['CTYPE2'] = 'DEC--TAN'



alma_b6_wcs = WCS(header_b6, naxis=2)    
#alma_b6_wcs = WCS(alma_b6_fits[0].header, naxis=2)


ref_fits = fits.open(image_filenames['f480m'])['SCI']
tgt_header = fits.getheader(image_filenames['f480m'], ext=('SCI', 1))

wcs_header = WCS(tgt_header, naxis=2)




pixscale_b3 = proj_plane_pixel_scales(alma_b3_wcs)  # in deg/pixel
pixscale_b6 = proj_plane_pixel_scales(alma_b6_wcs)  # in deg/pixel

# Image B's footprint in WCS
shape_ref = ref_fits.data.shape
ny, nx = shape_ref
wcs_header_corners = wcs_header.calc_footprint()



# Calculate bounding box of image B in world coordinates
ra_min, dec_min = wcs_header_corners.min(axis=0)
ra_max, dec_max = wcs_header_corners.max(axis=0)

# Define new output WCS using image A's pixel scale
cdelt_b3 = pixscale_b3  # [dy, dx] in deg/pixel from proj_plane_pixel_scales
cdelt_b6 = pixscale_b6

# Assign correctly:
naxis1_b3 = int(np.ceil((ra_max - ra_min) / cdelt_b3[1]))  # dx
naxis2_b3 = int(np.ceil((dec_max - dec_min) / cdelt_b3[0]))  # dy

naxis1_b6 = int(np.ceil((ra_max - ra_min) / cdelt_b6[1]))
naxis2_b6 = int(np.ceil((dec_max - dec_min) / cdelt_b6[0]))

# Build new header
new_wcs_b3 = WCS(naxis=2)
new_wcs_b3.wcs.crpix = [naxis1_b3 / 2, naxis2_b3 / 2]
new_wcs_b3.wcs.cdelt = [-cdelt_b3[1], cdelt_b3[0]]  # [dx, dy], negative for RA
new_wcs_b3.wcs.crval = [(ra_max + ra_min) / 2, (dec_max + dec_min) / 2]
new_wcs_b3.wcs.ctype = ["RA---TAN", "DEC--TAN"]

new_wcs_b6 = WCS(naxis=2)
new_wcs_b6.wcs.crpix = [naxis1_b6 / 2, naxis2_b6 / 2]
new_wcs_b6.wcs.cdelt = [-cdelt_b6[1], cdelt_b6[0]]
new_wcs_b6.wcs.crval = [(ra_max + ra_min) / 2, (dec_max + dec_min) / 2]
new_wcs_b6.wcs.ctype = ["RA---TAN", "DEC--TAN"]

# Reproject A onto the new WCS with A's resolution and B's footprint
shape_out_b3 = (naxis2_b3, naxis1_b3)  # (ny, nx)

shape_out_b3 = (naxis2_b3, naxis1_b3)
shape_out_b6 = (naxis2_b6, naxis1_b6)
print(alma_b3_image.shape)  # should be (ny, nx)
print(alma_b3_wcs.wcs.naxis)  # should be 2
print(alma_b3_wcs.wcs.ctype) 
array_b3, footprint_b3 = reproject_interp(
    (alma_b3_image, alma_b3_wcs), output_projection=new_wcs_b3, shape_out=shape_out_b3
)
shape_out_b6 = (naxis2_b6, naxis1_b6)
array_b6, footprint_b6 = reproject_interp(
    (alma_b6_image, alma_b6_wcs), output_projection=new_wcs_b6, shape_out=shape_out_b6
)



# Save the output if needed
hdu_out_b3 = fits.PrimaryHDU(data=array_b3, header=new_wcs_b3.to_header())
hdu_out_b3.writeto('fitsfiles/alma_b3_reproj.fits', overwrite=True)
hdu_out_b6 = fits.PrimaryHDU(data=array_b6, header=new_wcs_b6.to_header())
hdu_out_b6.writeto('fitsfiles/alma_b6_reproj.fits', overwrite=True)





