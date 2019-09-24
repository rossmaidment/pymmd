#========================================================================#
# Program to find match-ups between CCI and C3S data and SIRDS reference dataset
#------------------------------------------------------------------------#
# R. Maidment (July 2019)
#========================================================================#

# Import modules
import numpy as np
import itertools
import pandas as pd
from math import radians, cos, sin, asin, sqrt
from datetime import datetime as dt

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    r = 6372.8 # Radius of earth in kilometers. Use 3956 for miles
    return c * r

def get_matchups(idata, sdata, ftype, max_dist, max_time):
    # Create temporal limits including margins
    start_time_sec    = np.min(sdata['time_2d']).values.astype('datetime64[s]')
    end_time_sec      = np.max(sdata['time_2d']).values.astype('datetime64[s]')
    start_time_margin = start_time_sec - np.timedelta64((max_time)*(60*60), 's')
    end_time_margin   = end_time_sec + ((max_time)*(60*60))
    
    # Subset insitu data for this file
    start_time_margin_epoch = np.timedelta64((start_time_margin - np.datetime64(dt(1981,1,1,0,0,0))), 's').astype('int64')
    end_time_margin_epoch   = np.timedelta64((end_time_margin - np.datetime64(dt(1981,1,1,0,0,0))), 's').astype('int64')
    idata_sub = idata[(idata.time_seconds >= start_time_margin_epoch) & (idata.time_seconds <= end_time_margin_epoch)]
    idata_sub = idata_sub.reset_index(drop=True)
    
    # Identify unique in situ platforms for this orbit
    unique_callsign = list(set(idata_sub.callsign))
    unique_callsign.sort()
    print '%s unique callsigns found for this file' % len(unique_callsign)
    
    # Extract each variable (too slow to use xarray dataset directly)
    if ftype == 'l4':
        ""
    else:
        sdata_sea_surface_temperature             = sdata.sea_surface_temperature.values.squeeze()
        sdata_sea_surface_temperature_depth       = sdata.sea_surface_temperature_depth.values.squeeze()
        sdata_sst_dtime                           = sdata.sst_dtime.values.squeeze()
        sdata_sst_depth_dtime                     = sdata.sst_depth_dtime.values.squeeze()
        sdata_sses_bias                           = sdata.sses_bias.values.squeeze()
        sdata_sses_standard_deviation             = sdata.sses_standard_deviation.values.squeeze()
        sdata_sst_depth_total_uncertainty         = sdata.sst_depth_total_uncertainty.values.squeeze()
        sdata_l2p_flags                           = sdata.l2p_flags.values.squeeze()
        sdata_quality_level                       = sdata.quality_level.values.squeeze()
        sdata_wind_speed                          = sdata.wind_speed.values.squeeze()
        sdata_large_scale_correlated_uncertainty  = sdata.large_scale_correlated_uncertainty.values.squeeze()
        sdata_synoptically_correlated_uncertainty = sdata.synoptically_correlated_uncertainty.values.squeeze()
        sdata_uncorrelated_uncertainty            = sdata.uncorrelated_uncertainty.values.squeeze()
        sdata_adjustment_uncertainty              = sdata.adjustment_uncertainty.values.squeeze()
        sdata_aerosol_dynamic_indicator           = sdata.aerosol_dynamic_indicator.values.squeeze()
        sdata_sensitivity                         = sdata.sensitivity.values.squeeze()
        sdata_time_2d                             = sdata.time_2d.values.squeeze()
    
    # Loop over each unqiue callsign (note each callsign will have multiple records)
    matchup = list()
    for idx, item in enumerate(unique_callsign):
        #print idx, item
        idata_sub_cs = idata_sub[idata_sub.callsign == item]
        
        # Calculate in situ location on CCI/C3S grid
        lon = np.array(idata_sub_cs['lon'])
        lat = np.array(idata_sub_cs['lat'])
        slon_i = ((lon + 180) * 20).astype(int) % 7200
        slat_i = ((lat + 90) * 20).astype(int) % 3600
        slon   = sdata.lon.values[slon_i]
        slat   = sdata.lat.values[slat_i]
        itime  = np.array(idata_sub_cs.time)
        rownum = np.array(idata_sub_cs.index)
        
        # Check for valid CCI/C3S data
        if ftype == 'l4':
            idt = np.where(pd.notnull(sdata_time_2d[slat_i,slon_i]) &
                                 (sdata_quality_level[slat_i,slon_i] == 1))[0]
        else:
            idt = np.where(pd.notnull(sdata_time_2d[slat_i,slon_i]) &
                                 (sdata_quality_level[slat_i,slon_i] >= 2))[0]
        
        # Only continue if there are valid matches
        if len(idt) > 0:
            # Reduce to valid pixels
            lon    = lon[idt]
            lat    = lat[idt]
            slon_i = slon_i[idt]
            slat_i = slat_i[idt]
            slon   = slon[idt]
            slat   = slat[idt]
            itime  = itime[idt]
            rownum = rownum[idt]
        
            # Create unique key for each pixel
            num_loc = slon_i*slat_i
        
            # Find duplicate insitu data
            idz     = np.array(np.argsort(num_loc, kind='stable'))
            info    = reduce(lambda lst,items: lst + [(items[0], items[1], sum(map(lambda i: i[1], lst)))],[(key, len(list(it))) for (key, it) in itertools.groupby(num_loc[idz])], [])
            dups    = [x[1] for x in info if x[1] > 1]
            if len(set(num_loc)) == 1:
                first_ind = [-1]
                count     = [0]
            else:
                if len(dups) > 0:
                    first_ind = [x[2] for x in info if x[1] > 1]
                    count     = [x[1] for x in info if x[1] > 1]
                else:
                    first_ind = [-1]
                    count     = [0]
        
            # Remove duplicates
            mask = np.zeros(len(num_loc), dtype='int')
            if first_ind[0] != -1:
                for k, item_count in enumerate(count):
                    index = first_ind[k] + np.arange(0,count[k])
                    index = np.array(idz[index])
                    mask[index] = 1
                    time_diff = sdata_time_2d[slat_i[index[0]],slon_i[index[0]]] - itime[index]
                    tmin = np.min(np.abs(time_diff))
                    loct = np.argmin(np.abs(np.array(time_diff)))
                    mask[index[loct]] = 0
        
            # Keep non-duplicates only
            if any(mask == 0):
                lon    = lon[np.where(mask == 0)]
                lat    = lat[np.where(mask == 0)]
                slon_i = slon_i[np.where(mask == 0)]
                slat_i = slat_i[np.where(mask == 0)]
                slon   = slon[np.where(mask == 0)]
                slat   = slat[np.where(mask == 0)]
                itime  = itime[np.where(mask == 0)]
                rownum = rownum[np.where(mask == 0)]
        
            # Calculate distance
            distance = list()
            for idy in range(0,len(lon)):
                distance.append(haversine(lon[idy],lat[idy],slon[idy],slat[idy]))
            
            distance = np.array(distance)
        
            # Calculate time difference between in situ and satellite
            stime = sdata_time_2d[slat_i,slon_i]
            tdiff = (stime - itime).astype('timedelta64[s]')
        
            # Reject any match-ups outside of predefined limits
            idv = list(np.where((distance < max_dist) & (np.abs(tdiff.astype(int)) <= max_time*3600)))[0]
            #idata_nodups = idata_nodups.loc[(idata_nodups['distance'] < max_dist) & (np.abs(idata_nodups['tdiff_seconds']) <= max_time * 3600)]
        
            # Store any valid match-ups
            if len(idv) > 0:
                idata_valid = idata_sub_cs.loc[rownum[idv]]
                idata_valid['slat']     = slat[idv]
                idata_valid['slon']     = slon[idv]
                idata_valid['slat_i']   = slat_i[idv]
                idata_valid['slon_i']   = slon_i[idv]
                idata_valid['num_loc']  = num_loc[idv]
                idata_valid['distance'] = distance[idv]
                idata_valid['stime']    = stime[idv]
                idata_valid['tdiff']    = tdiff[idv]
            
                if ftype == 'l4':
                    ""
                else:
                    idata_valid['sea_surface_temperature']             = sdata_sea_surface_temperature[idata_valid['slat_i'],idata_valid['slon_i']]
                    idata_valid['sea_surface_temperature_depth']       = sdata_sea_surface_temperature_depth[idata_valid['slat_i'],idata_valid['slon_i']]
                    idata_valid['sst_dtime']                           = sdata_sst_dtime[idata_valid['slat_i'],idata_valid['slon_i']]
                    idata_valid['sst_depth_dtime']                     = sdata_sst_depth_dtime[idata_valid['slat_i'],idata_valid['slon_i']]
                    idata_valid['sses_bias']                           = sdata_sses_bias[idata_valid['slat_i'],idata_valid['slon_i']]
                    idata_valid['sses_standard_deviation']             = sdata_sses_standard_deviation[idata_valid['slat_i'],idata_valid['slon_i']]
                    idata_valid['sst_depth_total_uncertainty']         = sdata_sst_depth_total_uncertainty[idata_valid['slat_i'],idata_valid['slon_i']]
                    idata_valid['l2p_flags']                           = sdata_l2p_flags[idata_valid['slat_i'],idata_valid['slon_i']]
                    idata_valid['quality_level']                       = sdata_quality_level[idata_valid['slat_i'],idata_valid['slon_i']]
                    idata_valid['wind_speed']                          = sdata_wind_speed[idata_valid['slat_i'],idata_valid['slon_i']]
                    idata_valid['large_scale_correlated_uncertainty']  = sdata_large_scale_correlated_uncertainty[idata_valid['slat_i'],idata_valid['slon_i']]
                    idata_valid['synoptically_correlated_uncertainty'] = sdata_synoptically_correlated_uncertainty[idata_valid['slat_i'],idata_valid['slon_i']]
                    idata_valid['uncorrelated_uncertainty']            = sdata_uncorrelated_uncertainty[idata_valid['slat_i'],idata_valid['slon_i']]
                    idata_valid['adjustment_uncertainty']              = sdata_adjustment_uncertainty[idata_valid['slat_i'],idata_valid['slon_i']]
                    idata_valid['aerosol_dynamic_indicator']           = sdata_aerosol_dynamic_indicator[idata_valid['slat_i'],idata_valid['slon_i']]
                    idata_valid['sensitivity']                         = sdata_sensitivity[idata_valid['slat_i'],idata_valid['slon_i']]
                    idata_valid['time_2d']                             = sdata_time_2d[idata_valid['slat_i'],idata_valid['slon_i']]
                
                matchup.append(idata_valid)
    
    # Concatenate dataframe
    matchup = pd.concat(matchup)
    
    return(matchup, len(matchup))


#x = get_matchups(idata, sdata, ftype, max_dist, max_time)



#for idx, item in enumerate(idata_valid):
#    print idata_valid.iloc[:,idx]

#sdata_sub = sdata.isel(lon=idata_valid['slon_i'], lat=idata_valid['slat_i'])

#matchup = list()
#for idx, item in enumerate(unique_callsign[0:1]):
#    print idx, item
#    idata_sub_cs = idata_sub[idata_sub.callsign == item]
    #idata_sub_cs = idata_sub_cs.reset_index(drop=True)
    
    # Calculate in situ location on CCI/C3S grid
#    idata_sub_cs['slon_i'] = ((idata_sub_cs['lon'] + 180) * 20).astype(int) % 7200
#    idata_sub_cs['slat_i'] = ((idata_sub_cs['lat'] + 90) * 20).astype(int) % 3600
#    idata_sub_cs['slon'] = sdata.lon.values[idata_sub_cs['slon_i']]
#    idata_sub_cs['slat'] = sdata.lat.values[idata_sub_cs['slat_i']]
    
    # Check for valid CCI/C3S data
#    if ftype == 'l4':
#        idt = np.where(pd.notnull(sdata.time_2d.values.squeeze()[idata_sub_cs['slat_i'].values,idata_sub_cs['slon_i'].values]) &
#                                 (sdata.quality_level.values.squeeze()[idata_sub_cs['slat_i'].values,idata_sub_cs['slon_i'].values] == 1))[0]
#    else:
#        idt = np.where(pd.notnull(sdata.time_2d.values.squeeze()[idata_sub_cs['slat_i'].values,idata_sub_cs['slon_i'].values]) &
#                                 (sdata.quality_level.values.squeeze()[idata_sub_cs['slat_i'].values,idata_sub_cs['slon_i'].values] >= 2))[0]
    
    # Only continue if there are valid matches
#    if len(idt) > 0:
        # Reduce to valid pixels
#        idata_valid = idata_sub_cs.iloc[idt]
        
        # Create unique key for each pixel
#        idata_valid['num_loc'] = np.array(idata_valid.slon_i * idata_valid.slat_i)
        
        # Find duplicate insitu data
#        idz     = np.array(np.argsort(idata_valid.num_loc, kind='stable'))
#        info    = reduce(lambda lst,items: lst + [(items[0], items[1], sum(map(lambda i: i[1], lst)))],[(key, len(list(it))) for (key, it) in itertools.groupby(idata_valid.num_loc.iloc[idz])], [])
#        dups    = [x[1] for x in info if x[1] > 1]
#        if len(set(idata_valid['num_loc'])) == 1:
#            first_ind = [-1]
#            count     = [0]
#        else:
#            if len(dups) > 0:
#                first_ind = [x[2] for x in info if x[1] > 1]
#                count     = [x[1] for x in info if x[1] > 1]
#            else:
#                first_ind = [-1]
#                count     = [0]
#
        # Remove duplicates
#        mask = np.zeros(len(idata_valid.num_loc), dtype='int')
#        if first_ind[0] != -1:
#            for k, item_count in enumerate(count):
#                index = first_ind[k] + np.arange(0,count[k])
#                index = np.array(idz[index])
#                mask[index] = 1
#                time_diff = sdata.time_2d.values.squeeze()[idata_valid['slat_i'].values[index[0]],idata_valid['slon_i'].values[index[0]]] - idata_valid.time.iloc[index]
#                tmin = np.min(np.abs(time_diff))
#                loct = np.argmin(np.abs(np.array(time_diff)))
#                mask[index[loct]] = 0
        
        # Keep non-duplicates only
#        if any(mask == 0):
#            idata_nodups = idata_valid.iloc[np.where(mask == 0)]
        
        # Calculate distance
#        distances = list()
#        for idy in range(0,len(idata_nodups)):
#            distances.append(haversine(idata_nodups.lon.iloc[idy],idata_nodups.lat.iloc[idy],idata_nodups.slon.iloc[idy],idata_nodups.slat.iloc[idy]))
            
#        idata_nodups['distance'] = distances
        
        # Calculate time difference between in situ and satellite
#        idata_nodups['stime'] = sdata.time_2d.values.squeeze()[idata_nodups['slat_i'].values,idata_nodups['slon_i'].values]
#        idata_nodups['tdiff'] = idata_nodups['stime'] - idata_nodups.time
#        tdiff_hours   = list()
#        tdiff_seconds = list()
#        for time_idx in idata_nodups['tdiff']:
#            tdiff_hours.append(time_idx.total_seconds()/3600)
#            tdiff_seconds.append(time_idx.total_seconds())
        
#        idata_nodups['tdiff_hours']   = tdiff_hours
#        idata_nodups['tdiff_seconds'] = tdiff_seconds
        
        # Reject any match-ups outside of predefined limits
#        idata_nodups = idata_nodups.loc[(idata_nodups['distance'] < max_dist) & (np.abs(idata_nodups['tdiff_seconds']) <= max_time * 3600)]
        
        # Store any valid match-ups
#        if len(idata_nodups) > 0:
            #print '-> %s valid match-up(s) found for callsign %s' % (len(idata_nodups), item)
#            matchup.append(idata_nodups)
        
#    else:
        #print '-> No valid satellite values found for callsign %s' % item
