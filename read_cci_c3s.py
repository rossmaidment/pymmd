#========================================================================#
# Program to read fields from CCI/C3S files for generating match-ups
#------------------------------------------------------------------------#
# R. Maidment (July 2019)
#========================================================================#

# Import modules
import netCDF4
import numpy as np
import copy
import xarray as xr

# Read file
def read_satellite(filename, ftype):
    """Read l3c or l4 satellite file
    
    """
    #ftype = 'l3c'
    #filename = '/gws/nopw/j04/cds_c3s_sst/output/v2.6.0/l3c/AVHRR19_G/2018/03/01/20180301120000-C3S-L3C_GHRSST-SSTskin-AVHRR19_G-ICDR2.0_day-v02.0-fv01.0.nc'
    #ftype = 'l4'
    #filename = '/gws/nopw/j04/cds_c3s_sst/public/data/ICDR_v2/Analysis/L4/v2.0/2018/01/01/20180101120000-C3S-L4_GHRSST-SSTdepth-OSTIA-GLOB_ICDR2.0-v02.0-fv01.0.nc'
    print "Reading %s file: %s" % (ftype, filename)
    
    # Read data - L4 or L3C (note L4 mask and L3C quality level have same array name)
    ncin = netCDF4.Dataset(filename)
    if ftype == 'l4':
        lon          = ncin.variables['lon'][:]
        lat          = ncin.variables['lat'][:]
        time_read    = ncin.variables['time'][:]
        sst          = ncin.variables['analysed_sst'][:]
        unc          = ncin.variables['analysis_uncertainty'][:]
        sea_ice_frac = ncin.variables['sea_ice_fraction'][:]
        ql           = ncin.variables['mask'][:]
        sstfill      = ncin.variables['analysed_sst']._FillValue
        sstao        = ncin.variables['analysed_sst'].add_offset
        sstsf        = ncin.variables['analysed_sst'].scale_factor
    elif ftype == 'l3c':
        lon                 = ncin.variables['lon'][:]
        lat                 = ncin.variables['lat'][:]
        time_read           = ncin.variables['time'][:]
        time_bnds           = ncin.variables['time_bnds'][:]
        sst                 = ncin.variables['sea_surface_temperature'][:]
        sst_depth           = ncin.variables['sea_surface_temperature_depth'][:]
        sst_dtime           = ncin.variables['sst_dtime'][:]
        sst_depth_dtime     = ncin.variables['sst_depth_dtime'][:]
        sses_bias           = ncin.variables['sses_bias'][:]
        sses_sd             = ncin.variables['sses_standard_deviation'][:]
        sst_depth_total_unc = ncin.variables['sst_depth_total_uncertainty'][:]
        l2p_flags           = ncin.variables['l2p_flags'][:]
        ql                  = ncin.variables['quality_level'][:]
        wind_speed          = ncin.variables['wind_speed'][:]
        large_scale_cor_unc = ncin.variables['large_scale_correlated_uncertainty'][:]
        synop_cor_unc       = ncin.variables['synoptically_correlated_uncertainty'][:]
        uncor_unc           = ncin.variables['uncorrelated_uncertainty'][:]
        adj_unc             = ncin.variables['adjustment_uncertainty'][:]
        aerosol_dyn_ind     = ncin.variables['aerosol_dynamic_indicator'][:]
        sens                = ncin.variables['sensitivity'][:]
        tfill               = ncin.variables['sst_dtime']._FillValue
        sstfill             = ncin.variables['sea_surface_temperature']._FillValue
        sstao               = ncin.variables['sea_surface_temperature'].add_offset
        sstsf               = ncin.variables['sea_surface_temperature'].scale_factor
    else:
        print 'ftype not recognised or supported'
    
    # Create time field
    # -> If L4 then create a time field set to time in L4 file
    # -> Also add a time fill value to keep coding simple later on
    if ftype == 'l4':
        time = np.empty((7200,3600))
        time[:,:] = time_read
        tfill = -2147483648
    else:
        time = copy.deepcopy(sst_dtime) # Need to make a hard copy
        mask = sst_dtime.mask == False; mask = mask[0,:,:]
        row, col = np.where(mask==True)
        time.data[0, row, col] = time.data[0,row, col] + time_read
    
    # Create output structure
    if ftype == 'l4':
        data = dict(lon=lon,
                    lat=lat,
                    time_read=time_read,
                    time=time,
                    sst=sst,
                    unc=unc,
                    sea_ice_frac=sea_ice_frac,
                    ql=ql,
                    tfill=tfill,
                    sstfill=sstfill,
                    sstao=sstao,
                    sstsf=sstsf)
    elif ftype == 'l3c':
        data = dict(lon=lon,
                    lat=lat,
                    time_read=time_read,
                    time=time,
                    time_bnds=time_bnds,
                    sst=sst,
                    sst_depth=sst_depth,
                    sst_dtime=sst_dtime,
                    sst_depth_dtime=sst_depth_dtime,
                    sses_bias=sses_bias,
                    sses_sd=sses_sd,
                    sst_depth_total_unc=sst_depth_total_unc,
                    l2p_flags=l2p_flags,
                    ql=ql,
                    wind_speed=wind_speed,
                    large_scale_cor_unc=large_scale_cor_unc,
                    synop_cor_unc=synop_cor_unc,
                    uncor_unc=uncor_unc,
                    adj_unc=adj_unc,
                    aerosol_dyn_ind=aerosol_dyn_ind,
                    sens=sens,
                    tfill=tfill,
                    sstfill=sstfill,
                    sstao=sstao,
                    sstsf=sstsf)
    else:
        print 'ftype not recognised or supported'
    
    return data


def read_satellite_xarray(filename, ftype):
    """Read l3c or l4 satellite file
    
    """
    print "Reading %s file: %s" % (ftype, filename)
        
    # Read data - L4 or L3C (note L4 mask and L3C quality level have same array name)
    dataset = xr.open_mfdataset(filename, mask_and_scale=True)
    print 'File has %s variables' % str(len(dataset.keys()))
    
    # Create 2D time field
    if ftype == 'l4':
        time = np.empty((1, 3600, 7200), dtype='datetime64[ns]')
        time = dataset.time.values
    elif ftype == 'l3c':
        mask = np.isnan(dataset.sea_surface_temperature.values) == False
        time = dataset.sst_dtime.values
        time = time + dataset.time.values
    else:
        print 'ftype not recognised or supported'
    
    # Add 2D time field
    dataset["time_2d"] = (['time', 'lat', 'lon'],  time)
    
    return dataset

#ftype    = 'l3c'
#filename = '/gws/nopw/j04/cds_c3s_sst/output/v2.6.0/l3c/AVHRR19_G/2018/03/01/20180301120000-C3S-L3C_GHRSST-SSTskin-AVHRR19_G-ICDR2.0_day-v02.0-fv01.0.nc'
#ftype = 'l4'
#filename = '/gws/nopw/j04/cds_c3s_sst/public/data/ICDR_v2/Analysis/L4/v2.0/2018/01/01/20180101120000-C3S-L4_GHRSST-SSTdepth-OSTIA-GLOB_ICDR2.0-v02.0-fv01.0.nc'
#x = read_satellite_xarray(filename, ftype)

    
