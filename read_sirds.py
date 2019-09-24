#========================================================================#
# Program to read monthly SIRDS file from Met Office
#------------------------------------------------------------------------#
# R. Maidment
#========================================================================#

# Import modules
import netCDF4
import pandas as pd
import numpy as np
from datetime import datetime as dt
from datetime import timedelta as td

# Read file
def read_insitu(filename):
    """Read insitu file
    
    """
    #filename = '/gws/nopw/j04/cds_c3s_sst/input/refdata/raw/sirds/SSTCCI2_refdata_drifter_cmems_201803.nc'
    print "Reading file: %s" % filename
    ncin         = netCDF4.Dataset(filename)
    dim_name     = 'N_OBS'
    dim_size     = ncin.dimensions[dim_name].size
    ob_id        = ncin.variables['OB_ID'][:]
    collection   = ncin.variables['COLLECTION'][:]
    subcol1      = ncin.variables['SUBCOL1'][:]
    subcol2      = ncin.variables['SUBCOL2'][:]
    prof         = ncin.variables['PROF_ID'][:]
    callsign     = ncin.variables['PLAT_ID'][:]; callsign = [int("".join(map(str, i))) for i in callsign]
    lon          = ncin.variables['LONGITUDE'][:]
    lat          = ncin.variables['LATITUDE'][:]
    year         = ncin.variables['YEAR'][:]
    month        = ncin.variables['MONTH'][:]
    day          = ncin.variables['DAY'][:]
    hour         = ncin.variables['HOUR'][:]
    minute       = ncin.variables['MINUTE'][:]
    second       = ncin.variables['SECOND'][:]
    depth        = ncin.variables['DEPTH'][:]
    depth_corr   = ncin.variables['DEPTH_CORR'][:]
    sst          = ncin.variables['SST'][:]
    sst_type     = ncin.variables['SST_TYPE_CORR'][:]
    sst_type_unc = ncin.variables['SST_TYPE_CORR_UNC'][:]
    sst_plat     = ncin.variables['SST_PLAT_CORR'][:]
    sst_plat_unc = ncin.variables['SST_PLAT_CORR_UNC'][:]
    ran          = ncin.variables['SST_RAND_UNC'][:]
    sst_unc      = ncin.variables['SST_COMB_UNC'][:]
    qc1          = ncin.variables['QC1'][:]
    qc2          = ncin.variables['QC1'][:]
    
    # Get time in seconds since 1st Jan 1981
    ep = dt(1981,1,1,0,0,0)
    time_seconds = []
    time = []
    for idx, item in enumerate(year):
        date_tmp = dt.strptime(str(year[idx]) + "{:02d}".format(month[idx]) + "{:02d}".format(day[idx]) + ' ' + "{:02d}".format(hour[idx]) +':'+ "{:02d}".format(minute[idx]) +':'+ "{:02d}".format(second[idx]), "%Y%m%d %H:%M:%S")
        time_seconds.append(int((date_tmp - ep).total_seconds()))
        time.append((np.datetime64(ep + td(seconds=time_seconds[idx]))))
    
    # Quality control data - not yet implemented
    
    
    # Create output structure (pandas dataframe)
    data = {'ob_id': ob_id,
            'time_seconds': np.array(time_seconds),
            'time': np.array(time),
            'lat':lat,
            'lon':lon,
            'sst':sst,
            'sst_unc': sst_unc,
            'depth': depth,
            'callsign': callsign,
            'collection': collection,
            'subcol1': subcol1,
            'subcol2': subcol2,
            'prof': prof,
            'depth_corr': depth_corr,
            'sst_type': sst_type,
            'sst_type_unc': sst_type_unc,
            'sst_plat': sst_plat,
            'sst_plat_unc': sst_plat_unc,
            'ran': ran,
            'qc1': qc1,
            'qc2': qc2}
    
    data = pd.DataFrame(data)
    return(data)
