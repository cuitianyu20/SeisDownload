# This file is part of seisdownload, a Python package for seismic event download.
# Copyright (C) 2024 Tianyu Cui (tycuicn@gmail.com)


from obspy import UTCDateTime
from obspy.clients.fdsn import Client

'''
Get event catalog from IRIS
rectangular domain: [minlatitude, maxlatitude, minlongitude, maxlongitude] in degree
'''
def event_rect_domain(starttime, endtime, minlat, maxlat, minlon, maxlon, 
                 minmag=0, maxmag=10, min_dep=0, max_dep=None, client_ini="IRIS"):
    # start and end time of the event catalog
    starttime = UTCDateTime(starttime)
    endtime = UTCDateTime(endtime)
    # search for events from client
    client = Client(client_ini)                        # default client: IRIS client
    # allow span more than 180 degree in longitude, rectangular domain
    if maxlon > 180:
        events1 = client.get_events(starttime, endtime, mindepth=min_dep, maxdepth=max_dep, minlatitude=minlat, 
                                   maxlatitude=maxlat, minlongitude=minlon, maxlongitude=180, minmagnitude=minmag, maxmagnitude=maxmag)
        events2 = client.get_events(starttime, endtime, mindepth=min_dep, maxdepth=max_dep, minlatitude=minlat,
                                    maxlatitude=maxlat, minlongitude=-180, maxlongitude=maxlon-360, minmagnitude=minmag, maxmagnitude=maxmag)
        events_cat = events1 + events2
    else:
        events_cat = client.get_events(starttime, endtime, mindepth=min_dep, minlatitude=minlat, maxdepth=max_dep,
                                   maxlatitude=maxlat, minlongitude=minlon, maxlongitude=maxlon, minmagnitude=minmag, maxmagnitude=maxmag)
    return events_cat


'''
Get event catalog from IRIS
circular domain: [latitude, longitude, minradius, maxradius] in degree
'''
def event_circ_domain(starttime, endtime, cen_lat, cen_lon, minradius, maxradius, 
                 minmag=0, maxmag=10, min_dep=0, max_dep=None, client_ini="IRIS"):
    # start and end time of the event catalog
    starttime = UTCDateTime(starttime)
    endtime = UTCDateTime(endtime)
    # search for events from client
    client = Client(client_ini)                        # default client: IRIS client
    # circular domain
    events_cat = client.get_events(starttime, endtime, mindepth=min_dep, maxdepth=max_dep, latitude=cen_lat, longitude=cen_lon,
                                    minradius=minradius, maxradius=maxradius, minmagnitude=minmag, maxmagnitude=maxmag)
    return events_cat
