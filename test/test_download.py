# This file is part of seisdownload, a Python package for seismic event download.
# Copyright (C) 2024 Tianyu Cui (tycuicn@gmail.com)

# import seisdownload module
import seisdownload

### seismic waveform parameter setting ###
client_ini = "IRIS" # search event from this client, but search all station records from all clients
starttime = "2007-09-01 00:00:00"
endtime = "2007-09-01 23:59:59"
channel = ["*Z",]
wave_len = 1800
### event domain setting (circ, rect, catalog domain) ###
evt_domain_type = 'circ'
# 1. [ minlat, maxlat, minlon, maxlon ] for rect domain
# 2. [ center_lat, center_lon, minradius, maxradius ] for circ domain
# 3. arbitrary values for catalog domain
evt_domain_range = [0, 20, 0, 70]   
### station domain setting (circ, rect, globe domain) ###
sta_domain_type = 'circ'
# 1. [ center_lat, center_lon, minradius, maxradius ] for circ domain, 
# if event_circ_center is True, center point:(event_lat, event_lon) for arbitrary center point value
# if event_circ_center is False, center point:(center_lat, center_lon)
# 2. [ minlat, maxlat, minlon, maxlon ] for rect domain
# if rect_distance_limit is True, minradius = sta_min_dis_limit, maxradius = sta_max_dis_limit (default center point:(event_lat, event_lon))
# 3. arbitrary values for globe domain
sta_domain_range = [0, 0, 0, 70] 
### event catalog setting (depth, magnitude) ###
evt_min_dep = 0             # min event depth
evt_max_dep = None          # min event depth
evt_min_mag = 5.0           # min event magnitude
evt_max_mag = None          # max event magnitude
### event catalog file setting  ###
event_catalog = "event_with_PKiKP.csv"
### station catalog setting (array, station name, distance limit) ###
array_name = "*"            # network name
station_name = "*"          # station name
sta_min_dis_limit = 0       # epicentral distance min limit (rect station domain)
sta_max_dis_limit = 180     # epicentral distance max limit (rect station domain)
### other setting ###
event_circ_center = True    # center point is event location or not
rect_distance_limit = False # use distance limit or not (rect station domain)
delete_mseed = True         # delete raw miniseed file or not
remove_response = True      # remove response or not

### download seismic data ###
seisdownload.massdownload_data(starttime, endtime, channel, wave_len, evt_domain_type, evt_domain_range, sta_domain_type, sta_domain_range, 
                    evt_min_dep, evt_max_dep, evt_min_mag, evt_max_mag,event_catalog, client_ini, array_name, station_name, 
                    sta_min_dis_limit, sta_max_dis_limit, event_circ_center, rect_distance_limit, delete_mseed, remove_response)
