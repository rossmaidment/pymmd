#========================================================================#
# Program to generate macth-ups between CCI and C3S SST products
# and HadIOD in situ
#------------------------------------------------------------------------#
# R. Maidment
#========================================================================#

# Import modules
import netCDF4
from datetime import datetime as dt
import os
import os.path
import calendar
import fnmatch
import glob
import sys

#from pymmd import read_insitu - not working?
import read_sirds
import read_cci_c3s
import find_matchups

# Command line arguments
date_str = '201803'
sensor   = 'avhrr-n19'
ftype    = 'l3c'
rdac     = 'C3S'
version  = '2.6.0'
cdr      = 'ICDR2.0'
itype    = 'drifter_cmems'
print date_str, sensor, ftype, rdac, version, cdr, itype

print 'Running ESA SST CCI match-up software'
print 'Started at %s' % dt.now().strftime('%Y-%m-%d %H:%M:%S')

# Set time limits for allowable match-ups - in hours. 1 value for L4 for one others
if ftype == 'l4':
    max_time = 12
else:
    max_time = 2

# Create input paths - note CCI and C3S do not use same paths so calculated individually
ystr = date_str[0:4]
mstr = date_str[4:6]
if rdac == 'ESACCI':
    root = '/gws/nopw/j04/esacci_sst/'
elif rdac == 'C3S':
    root = '/gws/nopw/j04/cds_c3s_sst/'
else:
    print 'rdac type not recognised'

# Define paths - move to config file
cci_root    = '/gws/nopw/j04/esacci_sst/'
output_path = '/group_workspaces/jasmin2/nceo_uor/users/rmaidment/validation/output_pymmd/'
temp_path   = '/work/scratch-nompiio/rmaidment/temp_pymmd/'
ecmwf_path  = '/badc/ecmwf-era-interim/data/'

# Create directories
if not os.path.exists(temp_path):
    os.makedirs(temp_path)

if not os.path.exists(output_path):
    os.makedirs(output_path)

# Sensor specific settings
if sensor == 'aatsr-en':
    mmd_str = 'mmd03'
    sen_str = 'AATSR'
    fpath = os.path.join(root, 'output', cdr, 'ATSR', ftype.upper(), 'v'+version, sen_str, ystr, mstr)
elif sensor == 'slstr-s3a':
    mmd_str = 'mmd16'
    sen_str = 'SLSTRA'
    fpath = os.path.join(root, 'output', 'v'+version, ftype.lower(), sen_str, ystr, mstr)
elif sensor == 'ostia-cci':
    mmd_str = 'mmd17'
    sen_str = 'OSTIA-GLOB'
    fpath = os.path.join(root, 'output', cdr, 'Analysis', ftype.upper(), 'v'+version, sen_str, ystr, mstr)
elif sensor == 'avhrr-n19':
    mmd_str = 'mmd04'
    sen_str = 'AVHRR19_G'
    fpath = os.path.join(root, 'output', 'v'+version, ftype.lower(), sen_str, ystr, mstr)
elif sensor == 'avhrr-mta':
    mmd_str = 'mmd04'
    sen_str = 'AVHRRMTA_G'
    fpath = os.path.join(root, 'output', 'v'+version, ftype.lower(), sen_str, ystr, mstr)
else:
    print 'sensor not recongnised'

# Calculate DOY for storing in filename
start_doy = dt.strptime(ystr+mstr+'01', "%Y%m%d").timetuple().tm_yday
end_doy   = dt.strptime(ystr+mstr+"{:02d}".format(calendar.monthrange(int(ystr),int(mstr))[1]), "%Y%m%d").timetuple().tm_yday

# Define output folder - if it doesn't exist then create folder
out_path = os.path.join(output_path, 'validation/mms/mmd/', mmd_str, itype + '-sst_' + sensor)
if not os.path.exists(out_path):
    os.makedirs(out_path)
    print "Creating directory: %s" % out_path

# Define output filename
out_filename = mmd_str + '_sst_' + itype + '-sst_' + sensor + '_' + str(ystr) + '-' + str(start_doy) + '_' + str(ystr) + '-' + str(end_doy) + '_' + ftype + '.nc'
print out_filename

# Input CCI/C3S file types - spatial limit for match-ups defined here
if ftype == 'l2p':
    max_dist = 1.0
    fstr = str(ystr) + str(mstr) + '*-' + rdac + '-' + ftype.upper() + '_GHRSST-SSTskin-' + sen_str + '-' + cdr + '-v02.0-fv01.0.nc'
elif ftype == 'l3c':
    max_dist = 5.0
    fstr = str(ystr) + str(mstr) + '*-' + rdac + '-' + ftype.upper() + '_GHRSST-SSTskin-' + sen_str + '-' + cdr + '_*-v02.0-fv01.0.nc'
elif ftype == 'l4':
    max_dist = 5.0
    fstr = str(ystr) + str(mstr) + '*-' + rdac + '-' + ftype.upper() + '_GHRSST-SSTdepth-' + sen_str + '_' + cdr + '-v02.0-fv01.0.nc'
else:
    print 'ftype not recognised'

# Count files - if none then exit
fmatches = []
for rootdir, dirnames, filenames in os.walk(fpath):
    for filename in fnmatch.filter(filenames, fstr):
        fmatches.append(os.path.join(rootdir, filename))

if len(fmatches) == 0:
    print 'No files found for string %s, now exiting ...' % fstr
    sys.exit()
else:
    print '%s matches found for string %s' % (len(fmatches), fstr)

# Read in situ data - if none then exit
insitu_path = os.path.join(root, 'input/refdata/raw/sirds')
insitu_file = os.path.join(insitu_path, 'SSTCCI2_refdata_' + itype + '_'+ date_str + '.nc')
if os.path.isfile(insitu_file) == False:
    print 'No in situ data for this date'
    sys.exit()

idata = read_sirds.read_insitu(insitu_file)

# Now loop through each satellite file and find matchups
for idx, file in enumerate(fmatches):
    # Read CCI/C3S file
    sdata = read_cci_c3s.read_satellite_xarray(file, ftype)
    
    # Now find match ups
    matches, num_matches = find_matchups.get_matchups(idata, sdata, ftype, max_dist, max_time)
    
    # Only find matches if data exists
    #if sdata.time.values[0] != -1:
        
        # Print the variables found
        #print '-----------------------------'
        #print 'Following variables found:'
        #for name in sdata.items():
        #    print '-> %s ' % name[0]
        #
        #print '-----------------------------'


    
    # Store matches
    











#------------------------------------------------------------------------#
# Import modules
#------------------------------------------------------------------------#
from datetime import datetime as dt
from datetime import timedelta
import os

import time
import subprocess
import uuid

#------------------------------------------------------------------------#
# Define functions
#------------------------------------------------------------------------#
def tidyup_job(command):
    """Removes certain characters from job description.

    Args:
        command (list): Command line arguments

    Returns:
        command (dateTimeObject): Same as input without unwanted characters.

    """
    command=str(command).replace("'","")
    command=str(command).replace(",","")
    command=command[1:-1]
    return command

#------------------------------------------------------------------------#
# Set paths
#------------------------------------------------------------------------#
codedir  = '/group_workspaces/jasmin2/nceo_uor/users/rmaidment/validation/pvir/idl/mmd_gen'
logdir   = '/group_workspaces/jasmin2/nceo_uor/users/rmaidment/validation/log'

#------------------------------------------------------------------------#
# Determine path for log files and create if necessary
#------------------------------------------------------------------------#
today = dt.now().replace(microsecond=0)
logpath_day = os.path.join(logdir, str(today.year), str('%02d' % today.month), str('%02d' % today.day))
if not os.path.exists(logpath_day):
    os.makedirs(logpath_day)

#------------------------------------------------------------------------#
# Specifiy arguments for mmd_gen.sav
#------------------------------------------------------------------------#
#year    = range(2017,2019)
#month   = ['01','02','03','04','05','06','07','08','09','10','11','12']
#sensor  = 'slstr-s3a'
#ftype   = 'l3c'
#rdac    = 'C3S'
#version = '2.6.0'
#cdr     = 'ICDR2.0'
#itype   = 'drifter_cmems'

year    = range(2017,2019)
month   = ['01','02','03','04','05','06','07','08','09','10','11','12']
sensor  = 'avhrr-n19'
ftype   = 'l3c'
rdac    = 'C3S'
version = '2.6.0'
cdr     = 'ICDR2.0'
itype   = 'drifter_cmems'

#------------------------------------------------------------------------#
# Submit job for each year
#------------------------------------------------------------------------#
for yyyy in year:
    for mm in month:
        date_str = str(yyyy)+str(mm)
        print "Generating mmd file for: %s" % date_str
        uniqid = str(uuid.uuid4().fields[-1])[:8]
        uniqid = date_str
        idlcommand = '/apps/exelis/idl85/bin/idl -rt=' + codedir + '/mmd_gen.sav -args', date_str, sensor, ftype, rdac, version, cdr, itype, uniqid
        job = ['bsub',
                '-q', 'short-serial',
                '-W', '4:00',
                '-o', os.path.join(logpath_day,'mmd_gen_'+date_str+'_'+sensor+'_'+ftype+'_'+rdac+'_'+version+'_'+cdr+'_'+itype+'_'+'%J.out'),
                '-e', os.path.join(logpath_day,'mmd_gen_'+date_str+'_'+sensor+'_'+ftype+'_'+rdac+'_'+version+'_'+cdr+'_'+itype+'_'+'%J.err'),
                '-R', 'rusage[mem=60000]',
                '-M', '60000',
                tidyup_job(idlcommand)]
        
        result = subprocess.check_output(tidyup_job(job), shell=True)
        print result
