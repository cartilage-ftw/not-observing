import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astroquery.xmatch import XMatch
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astroquery.hips2fits import hips2fits
from matplotlib.patches import Rectangle
from astropy.table import Table
from dustmaps.sfd import SFDQuery
#from dustmaps.config import config

#config.reset() # to create ~/.dustmapsrc config file.


def load_targets(path='./targets.txt'):
    return pd.read_csv(path, header=1, sep='\t')

def load_list_terese(path):
    return pd.read_csv(path, sep='\t')

def fetch_photometry(ra_list, dec_list):
    coords_df = Table({"ra": ra_list, "dec":dec_list})
    result = XMatch.query(
        cat1=coords_df,
        cat2="vizier:II/336/apass9",
        max_distance=2.0 * u.arcsec,
        colRA1="ra",
        colDec1="dec",
    )
    return result.to_pandas()

def get_dust_extinction(phot_data):
    sfd = SFDQuery()

    coords = SkyCoord(
    ra=phot_data["ra"].values * u.deg,
    dec=phot_data["dec"].values * u.deg,
    frame="icrs",
    )
    ebv_sfd = sfd(coords)

    # Extinction estimates
    A_V = 2.742 * ebv_sfd
    A_B = 3.626 * ebv_sfd

    # Dereddened magnitudes
    B_0 = phot_data["Bmag"] - A_B
    V_0 = phot_data["Vmag"] - A_V

def make_fast_not_chart(
    coord=None,
    target_name='M51',
    fov_arcmin=8.0,
    hips_survey='CDS/P/DSS2/red',  # other bands could be 'CDS/P/PanSTARRS/DR1/g', 'CDS/P/SDSS/g', etc.
    show_alfosc_fov=True,
    output_file=None,
    comments=None,
):
    """
    Generates a high-contrast finding chart in ~1 second via CDS HiPS2FITS.
    """
    if coord is None:
        coord = SkyCoord.from_name(target_name)
        title_name = target_name
    else:
        title_name = f"RA {coord.ra.deg:.4f}, Dec {coord.dec.deg:.4f}"
    if target_name == None:
        target_name = title_name

    print("Fetching image cutout from CDS...")
    hdu_list = hips2fits.query(
        hips=hips_survey,
        ra=coord.ra,
        dec=coord.dec,
        fov=fov_arcmin * u.arcmin,
        width=800,
        height=800,
        projection="TAN"
    )
    
    hdu = hdu_list[0]
    wcs = WCS(hdu.header)
    data = hdu.data

    fig = plt.figure(figsize=(8, 8), dpi=300)
    ax = fig.add_subplot(1, 1, 1, projection=wcs)

    # Inverted grayscale with percentile contrast stretch
    vmin, vmax = np.percentile(data[~np.isnan(data)], [1, 99.5])
    ax.imshow(data, origin='lower', cmap='gray_r', vmin=vmin, vmax=vmax, aspect='equal')

    # Target crosshair with central gap
    target_x, target_y = wcs.world_to_pixel(coord)
    gap = 12        # Gap around target in pixels
    arm_length = 40 # Arm length in pixels
    crosshair_kwargs = dict(color='#e74c3c', lw=1.2, zorder=5)
    
    # Horizontal arms
    ax.plot([target_x - arm_length, target_x - gap], [target_y, target_y], **crosshair_kwargs)
    ax.plot([target_x + gap, target_x + arm_length], [target_y, target_y], **crosshair_kwargs)
    # Vertical arms
    ax.plot([target_x, target_x], [target_y - arm_length, target_y - gap], **crosshair_kwargs)
    ax.plot([target_x, target_x], [target_y + gap, target_y + arm_length], **crosshair_kwargs)

    # FoV Overlay (for ALFOSC it's 6.4' x 6.4', for StanCam it's 3' x 3')
    if show_alfosc_fov:
        alf_size = 6.4 * u.arcmin
        bl_coord = SkyCoord(coord.ra - alf_size/2, coord.dec - alf_size/2)
        tr_coord = SkyCoord(coord.ra + alf_size/2, coord.dec + alf_size/2)
        
        bl_x, bl_y = wcs.world_to_pixel(bl_coord)
        tr_x, tr_y = wcs.world_to_pixel(tr_coord)
        
        alf_rect = Rectangle(
            (bl_x, bl_y), tr_x - bl_x, tr_y - bl_y,
            edgecolor='#27ae60', facecolor='none', lw=1.5, ls='--',
            label="ALFOSC FoV (6.4' × 6.4')"
        )
        ax.add_patch(alf_rect)
        ax.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.9, fontsize=9)

    # 6. Formatting & Grid
    survey_name = hips_survey.split('/')[-1].upper()
    ax.set_title(
        f"NOT Finding Chart: {target_name} ({coords})\n"
        f"Survey: {survey_name} | FoV: {fov_arcmin:.1f}′ × {fov_arcmin:.1f}′\n"
        f"Comments: {comments}",
        fontsize=11, fontweight='bold', pad=12
    )
    ax.set_xlabel('Right Ascension (J2000)', fontsize=10)
    ax.set_ylabel('Declination (J2000)', fontsize=10)
    ax.coords.grid(color='gray', alpha=0.4, linestyle=':')

    plt.tight_layout()
    #if output_file:
    if target_name is not None:
        output_file = target_name + '.png'
    else:
        output_file = title_name + '.png'
    plt.savefig('finding_charts/' + output_file, dpi=300, bbox_inches='tight')
    print(f"Saved chart to {output_file}")

if __name__ == "__main__":

    # the weird tab-separated format I use
    targets = load_targets('targets.txt')
    for i in range(len(targets)):
        coords, target_name, comments = targets.iloc[i].values
        make_fast_not_chart(SkyCoord(coords, unit=(u.hourangle, u.deg)), target_name, comments=comments)
        coord = SkyCoord(coords, unit=(u.hourangle, u.deg))
    # the target lists that Terese gave us.
    for rpa_list in ['RPA_targets.tab', 'RPA_Bright_targets.tab']:
        targets = load_list_terese(rpa_list)

        # try to fetch photometry
        coords = SkyCoord(targets['RA'], targets['DEC'], unit=(u.hourangle, u.deg))
        phot_data = fetch_photometry(coords.ra.deg, coords.dec.deg)
        print(phot_data.head())

        #get_dust_extinction(phot_data)

        for i in range(len(targets)):
            target_id, ra, dec, v_mag, fe_h = targets.iloc[i].values
            coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
            comments = f'[Fe/H]={fe_h}, V={v_mag}'
            make_fast_not_chart(SkyCoord(" ".join([ra, dec]),unit=(u.hourangle, u.deg)),
                                         target_id, comments=comments)