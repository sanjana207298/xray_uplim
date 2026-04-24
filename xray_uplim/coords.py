"""
xray_uplim.coords
-----------------
Shared coordinate utilities used by all observatory modules.

parse_coord()
    Parse an (RA, Dec) string or float pair into an astropy SkyCoord.

sky_to_evt_pixel()
    Convert RA/Dec to sky pixel (X, Y) in an event file using per-column
    WCS keywords (TCRPXn, TCRVLn, TCDLTn).  Works for any observatory
    whose event files use this FITS convention — confirmed for both
    NuSTAR and XMM-Newton EPIC.

sky_to_img_pixel()
    Convert RA/Dec to pixel position in a standard FITS image (e.g. an
    exposure map or vignetting map) using the primary-array WCS.
    Works for any observatory whose image products carry standard
    CRPIX/CRVAL/CDELT or CD-matrix WCS keywords.
"""

import numpy as np
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u


def parse_coord(ra_str, dec_str):
    """
    Parse an (RA, Dec) pair into an astropy SkyCoord.

    Accepts
    -------
    Decimal degrees  :  304.297   58.202
    Sexagesimal      : "20:17:11.360"   "+58:12:08.10"
    astropy strings  : "20h17m11.36s"  "+58d12m08.1s"

    Returns
    -------
    astropy.coordinates.SkyCoord (ICRS)
    """
    try:
        return SkyCoord(ra=float(ra_str)*u.deg,
                        dec=float(dec_str)*u.deg, frame='icrs')
    except (ValueError, TypeError):
        pass
    ra_fmt  = str(ra_str).replace(':', 'h', 1).replace(':', 'm', 1) + 's'
    dec_fmt = str(dec_str).replace(':', 'd', 1).replace(':', 'm', 1) + 's'
    return SkyCoord(ra_fmt, dec_fmt, frame='icrs')


def sky_to_evt_pixel(ra_deg, dec_deg, evt_hdr):
    """
    Convert RA/Dec (degrees) to event-file sky pixel (X, Y).

    Reads per-column WCS keywords from the EVENTS binary table header:
    TCTYPn, TCRPXn, TCRVLn, TCDLTn for the columns named 'X' and 'Y'.
    This convention is used by NuSTAR and XMM-Newton EPIC event files
    (and likely other missions).

    WHY NOT astropy WCS(evt_hdr, naxis=2)?
    ---------------------------------------
    The primary-array WCS keywords (CRPIX/CRVAL/CDELT) in event files
    are absent or refer to something unrelated.  The correct sky
    projection is stored in the per-column keywords above.

    Parameters
    ----------
    ra_deg, dec_deg : float
        Source sky position in decimal degrees (ICRS).
    evt_hdr : astropy.io.fits.Header
        Header of the EVENTS binary table extension.

    Returns
    -------
    cx, cy  : float  — sky pixel coordinates matching the X/Y column values
    pscale  : float  — pixel scale in arcsec/pix (from the Dec/Y axis)
    """
    # Find column indices labelled 'X' and 'Y'
    x_col = y_col = None
    for i in range(1, 300):
        if f'TTYPE{i}' not in evt_hdr:
            break
        name = evt_hdr[f'TTYPE{i}'].strip().upper()
        if name == 'X':
            x_col = i
        elif name == 'Y':
            y_col = i

    if x_col is None or y_col is None:
        raise RuntimeError(
            f"Could not find X or Y columns in EVENTS header "
            f"(X_col={x_col}, Y_col={y_col}). "
            "Is this a valid cleaned event file?")

    try:
        crpx_x = float(evt_hdr[f'TCRPX{x_col}'])   # reference pixel (RA)
        crvl_x = float(evt_hdr[f'TCRVL{x_col}'])   # reference RA  (deg)
        cdlt_x = float(evt_hdr[f'TCDLT{x_col}'])   # deg/pix — negative for RA
        crpx_y = float(evt_hdr[f'TCRPX{y_col}'])   # reference pixel (Dec)
        crvl_y = float(evt_hdr[f'TCRVL{y_col}'])   # reference Dec (deg)
        cdlt_y = float(evt_hdr[f'TCDLT{y_col}'])   # deg/pix — positive
    except KeyError as exc:
        raise RuntimeError(
            f"Missing column WCS keyword {exc} in EVENTS header. "
            "File may be non-standard or from a different pipeline.") from exc

    # Linear TAN projection — dRA scaled by cos(Dec_ref) for foreshortening
    cos_dec = np.cos(np.radians(crvl_y))
    cx = crpx_x + (ra_deg  - crvl_x) * cos_dec / cdlt_x
    cy = crpx_y + (dec_deg - crvl_y)            / cdlt_y

    pscale = abs(cdlt_y) * 3600.0   # arcsec/pixel from Dec axis

    return cx, cy, pscale


def sky_to_img_pixel(ra_deg, dec_deg, img_hdr):
    """
    Convert RA/Dec to pixel position in a standard FITS image.

    Uses the primary-array WCS (CRPIX/CRVAL/CDELT or CD matrix), which
    works correctly for exposure maps, vignetting maps, and similar image
    products from any observatory.

    Parameters
    ----------
    ra_deg, dec_deg : float
    img_hdr         : astropy.io.fits.Header  (primary image HDU)

    Returns
    -------
    cx, cy  : float
    pscale  : float  — arcsec/pix
    """
    import warnings

    wcs    = WCS(img_hdr, naxis=2)
    cx, cy = wcs.all_world2pix([[ra_deg, dec_deg]], 0)[0]

    pscale = None
    for key in ('CDELT2', 'CD2_2'):
        if key in img_hdr:
            pscale = abs(img_hdr[key]) * 3600.0
            break
    if pscale is None:
        warnings.warn(
            "Pixel scale not found in image header (tried CDELT2, CD2_2). "
            "Falling back to NuSTAR default 2.459 \"/pix — "
            "set the correct value if using a different observatory.",
            RuntimeWarning, stacklevel=2)
        pscale = 2.459

    return cx, cy, pscale
