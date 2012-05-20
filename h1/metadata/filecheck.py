"""
Code to return metadata for specified mime types.
"""

import os
import collections
from scipy.io import netcdf as scipy_netcdf


class UnknownMetadataError(Exception):
    """Cannot make sense of file metadata.
    
    The file smells like the correct mimetype, but cannot 
    understand what the metadata is. 
    """
    pass

def png(filename):
    
    try:
        return tryMagWellPNG(filename)
    except UnknownMetadataError:
        pass

    try:
        return tryPoincarePNG(filename)
    except UnknownMetadataError:
        pass

    try:
        return tryIotaBarPNG(filename)
    except UnknownMetadataError:
        pass


    raise NotImplementedError("File metadata not supported.")
    
def netcdf(filename):
    
    try:
        return tryBOOZER(filename)
    except UnknownMetadataError:
        pass

    try:
        return tryVMEC(filename)
    except UnknownMetadataError:
        pass

    raise NotImplementedError("File metadata not supported.")

def svg(filename):

    try:
        return tryMagWellSVG(filename)
    except UnknownMetadataError:
        pass

    try:
        return tryPoincareSVG(filename)
    except UnknownMetadataError:
        pass

    try:
        return tryIotaBarSVG(filename)
    except UnknownMetadataError:
        pass
    
    raise NotImplementedError("File metadata not supported.")


def tryPoincareSVG(filename):
    basename = os.path.basename(filename)
    if not '-phi' in basename:
        raise UnknownMetadataError

    filetype = "Poincare plot (SVG)"
    kh_str = basename[2:6]
    beta_av = 0
    phi = float(basename.split('.')[0][10:])
    ret = collections.OrderedDict([
            ("filename", basename),
            ("filetype", filetype),
            ("beta_av", beta_av),
            ("k_h", float(kh_str)),
            ("phi", phi),
            ])
    return ret

def tryMagWellSVG(filename):
    basename = os.path.basename(filename)
    if not basename.startswith("magwell_"):
        raise UnknownMetadataError

    filetype = "Magnetic Well (SVG)"
    kh_str = basename.split("_")[1][2:-4]
    beta_av = 0
    
    ret = collections.OrderedDict([
            ("filename", basename),
            ("filetype", filetype),
            ("beta_av", beta_av),
            ("k_h", float(kh_str)),
            ])
    return ret

def tryIotaBarSVG(filename):
    basename = os.path.basename(filename)
    if not basename.startswith("iotabar_"):
        raise UnknownMetadataError

    filetype = "Rotational Transform (SVG)"
    kh_str = basename.split("_")[1][2:-4]
    beta_av = 0
    
    ret = collections.OrderedDict([
            ("filename", basename),
            ("filetype", filetype),
            ("beta_av", beta_av),
            ("k_h", float(kh_str)),
            ])
    return ret


def tryMagWellPNG(filename):
    basename = os.path.basename(filename)
    if not basename.startswith("magwell_"):
        raise UnknownMetadataError
    if not filename.endswith("_large.png"):
        raise UnknownMetadataError

    filetype = "Magnetic Well (PNG)"
    kh_str = basename.split("_")[1][2:]
    beta_av = 0
    
    ret = collections.OrderedDict([
            ("filename", basename),
            ("filetype", filetype),
            ("beta_av", beta_av),
            ("k_h", float(kh_str)),
            ])
    return ret

def tryPoincarePNG(filename):
    basename = os.path.basename(filename)
    if not '-phi' in basename:
        raise UnknownMetadataError
    if not filename.endswith("_large.png"):
        raise UnknownMetadataError

    filetype = "Poincare plot (PNG)"
    kh_str = basename[2:6]
    beta_av = 0
    phi = float(basename.split('_')[0][10:])
    ret = collections.OrderedDict([
            ("filename", basename),
            ("filetype", filetype),
            ("beta_av", beta_av),
            ("k_h", float(kh_str)),
            ("phi", phi),
            ])
    return ret

def tryIotaBarPNG(filename):
    basename = os.path.basename(filename)
    if not basename.startswith("iotabar_"):
        raise UnknownMetadataError
    if not filename.endswith("_large.png"):
        raise UnknownMetadataError

    filetype = "Rotational Transform (PNG)"
    kh_str = basename.split("_")[1][2:]
    beta_av = 0
    
    ret = collections.OrderedDict([
            ("filename", basename),
            ("filetype", filetype),
            ("beta_av", beta_av),
            ("k_h", float(kh_str)),
            ])
    return ret



def tryBOOZER(ncfile):
    """
    Check if it is a Boozer file and extract some data.


    INPUT:
      ncfile .... name of nc file

    RETURN:
      dictionary with some extracted data
    """
    
    ncd = scipy_netcdf.netcdf_file(ncfile, 'r')

    if (not (ncd.variables.has_key("rmnc_b") and
             ncd.variables.has_key("zmns_b") and
             ncd.variables.has_key("bmnc_b") )
        ):
        ncd.close()
        raise UnknownMetadataError

    ret = {}
    filename = os.path.basename(ncfile)
    filetype = "BOOZER"

    try:
        mpol = ncd.variables['mbooz_b'].getValue()
    except KeyError:
        mpol = None

    try:
        ntor = ncd.variables['nbooz_b'].getValue()
    except KeyError:
        ntor = None

    # beta_b is a half mesh quantity
    try:
        beta_av = sum(ncd.variables['beta_b'].data) / (ncd.variables['beta_b'].shape[0]-1)
    except:
        beta_av = None

    # iota_b is a half mesh quantity
    try:
        iota_av = sum(ncd.variables['iota_b'].data) / (ncd.variables['iota_b'].shape[0]-1)
    except:
        iota_av = None

    ncd.close()

    ret = collections.OrderedDict([
            ("filename", filename),
            ("filetype", filetype),
            ("beta_av", beta_av),
            ("iota_av", iota_av),
            ("mpol", mpol),
            ("ntor", ntor)
            ])

    return(ret)


def tryVMEC(ncfile):
    """
    Check if it is a VMEC file and extract some data.


    INPUT:
      ncfile .... name of nc file

    RETURN:
      dictionary with some extracted data
    """

    ncd = scipy_netcdf.netcdf_file(ncfile, 'r')

    if (not (ncd.variables.has_key("rmnc") and
             ncd.variables.has_key("zmns") and
             ncd.variables.has_key("bmnc") )
        ):
        ncd.close()
        raise UnknownMetadataError

    ret = {}
    filename = os.path.basename(ncfile)
    filetype = "VMEC"
    
    try:
        phiedge = ncd.variables['phi'][-1]
    except KeyError:
        phiedge = None

    try:
        mpol = ncd.variables['mpol'].getValue()
    except KeyError:
        mpol = None

    try:
        ntor = ncd.variables['ntor'].getValue()
    except KeyError:
        ntor = None

    # beta_vol is a half mesh quantity
    try:
        beta_av = sum(ncd.variables['beta_vol'].data) / (ncd.variables['beta_vol'].shape[0]-1)
    except:
        beta_av = None

    # iotaf is a fill mesh quantity
    try:
        iota_av = sum(ncd.variables['iotaf'].data) / ncd.variables['iotaf'].shape[0]
    except:
        iota_av = None

    if(ncd.variables['nextcur'].getValue() == 0):
        extcur = []
    else:
        extcur = ncd.variables['extcur'].data.tolist()

    mgrid_filename = "".join(ncd.variables['mgrid_file'].data).strip()
    ncd.close()
    # sometimes in VMEC netCDF output we may find mgrid_file = "none    "
    # This we translate to mgrid_file = "None"
    if(mgrid_filename == 'none'):
        mgrid_filename = 'None'
        mgrid_file_exists = False
    else:
        mgrid_file_exists = os.path.exists(mgrid_filename)

    raw_coil_cur = []
    coil_group = []

    if(mgrid_file_exists):
        try:
            nc_mgrid = scipy_netcdf.netcdf_file(mgrid_filename, 'r')
            raw_coil_cur = nc_mgrid.variables['raw_coil_cur'].data.tolist()
            coil_group = nc_mgrid.variables['coil_group'].data
            nc_mgrid.close()
            print "  data from mgrid file {0:s} loaded".format(mgrid_filename)
        except:
            raw_coil_cur = []
    else:
        raw_coil_cur = []

    coil_cur_currents = collections.OrderedDict()
    for i,coiln in enumerate(coil_group):
        coil_cur_currents.update({"".join(coiln).strip():raw_coil_cur[i]})

    mgrid_file = collections.OrderedDict([
            ("mgrid_filename", mgrid_filename),
            ("coil_currents", coil_cur_currents)
            ])

    ret = collections.OrderedDict([
            ("filename", filename),
            ("filetype", filetype),
            ("beta_av", beta_av),
            ("iota_av", iota_av),
            ("mpol", mpol),
            ("ntor", ntor),
            ("phiedge", phiedge),
            ("extcur", extcur),
            ("mgrid_file", mgrid_file)
            ])

    return(ret)


def tryBline(ncfile):
    """

    to be implemented!!!

    Check if it is a Bline model file and extract some data.


    INPUT:
      ncfile .... name of nc file

    RETURN:
      dictionary with some extracted data
    """

    try:
        ncd = netCDF4.Dataset(ncfile, 'r')
    except:
        raise Exception, 'no netCDF3/4 file'

    ncd.close()

    return(ret)

