# This file is part of seisdownload, a Python package for seismic event download.
# Copyright (C) 2024 Tianyu Cui (tycuicn@gmail.com)

import os
import lxml
from obspy.core.event.catalog import read_events
from .mseed2sac import miniseed2sac
from .station_domain import station_dif_domain
from .event_catalog_info import catcsv2xml, save_events_info
from .event_domain import event_rect_domain, event_circ_domain
from obspy.clients.fdsn.mass_downloader import Restrictions, MassDownloader


'''
MassDownload_data: download waveform data !!!
Download Minisedd data and convert to SAC format.

input parameters:
param starttime: start time of waveform data
param endtime: end time of waveform data
param channel: channels (list)
param wave_len: downloaded waveform length in seconds
param evt_domain_type: event domain type (1.circ: circular domain; 
                                          2.rect: rectangular domain; 
                                          3.catalog: read event catalog file (csv or other obspy format))
param evt_domain_range: event domain range (circular domain: [latitude, longitude, minradius, maxradius] in degree;
                                            rectangular domain: [minlatitude, maxlatitude, minlongitude, maxlongitude] in degree)
param sta_domain_type: station domain type (1.circ: circular domain; 
                                            2.rect: rectangular domain; 
                                            3.globe: global domain)
param sta_domain_range: station domain range (circular domain: [latitude, longitude, minradius, maxradius] in degree;
                                                rectangular domain: [minlatitude, maxlatitude, minlongitude, maxlongitude] in degree;
                                                global domain: [])
param evt_min_dep: min event depth in km
param evt_max_dep: max event depth in km
param evt_min_mag: min event magnitude
param evt_max_mag: max event magnitude
param event_catalog: event catalog file name
param client_ini: client name (default: IRIS)
param array_name: array name (default: "*")
param station_name: station name (default: "*")
param sta_min_dis_limit: min distance in degree (default: 0)
param sta_max_dis_limit: max distance in degree (default: 180)
param event_circ_center: if True, use event center as reference point (default: False)
param rect_distance_limit: if True, add distance restriction to the Rectangular domain (default: False)
param delete_mseed: if True, delete corresponding miniseed data if miniseed convert to sac successfully (default: True)
param remove_response: if True, remove instrument response (default: False)

output parameters:
1. waveform_mseed: directory of waveform data (miniseed format)
2. waveform_sac: directory of waveform data (sac format)
3. station_inv: directory of station inventory (stationxml format)
'''
def massdownload_data(starttime, endtime, channel, wave_len, evt_domain_type, evt_domain_range, sta_domain_type, sta_domain_range, 
                      evt_min_dep=0, evt_max_dep=None, evt_min_mag=None, evt_max_mag=None, event_catalog=None, client_ini="IRIS", 
                      array_name="*", station_name="*", sta_min_dis_limit=0, sta_max_dis_limit=180, event_circ_center=True, 
                      rect_distance_limit=False, delete_mseed=True, remove_response=False):
    ######  Module 1: Get event catalog from client  ######
    # convert event domain type to number
    event_domain_number = {'circ': 1, 'rect': 2, 'catalog': 3}
    if evt_domain_type in event_domain_number:
        evt_domain_type = event_domain_number[evt_domain_type]
    else:
        raise SystemExit('Event domain type error!')
    # circular domain: [latitude, longitude, minradius, maxradius] in degree
    if evt_domain_type == 1:
        events_cat = event_circ_domain(starttime, endtime, cen_lat=evt_domain_range[0], cen_lon=evt_domain_range[1], minradius=evt_domain_range[2], maxradius=evt_domain_range[3], 
                                   minmag=evt_min_mag, maxmag=evt_max_mag, min_dep=evt_min_dep, max_dep=evt_max_dep, client_ini=client_ini)
        # save event catalog and plot event map
        save_events_info(events_cat, save_path='.')
    # rectangular domain: [minlatitude, maxlatitude, minlongitude, maxlongitude] in degree
    elif evt_domain_type == 2:
        events_cat = event_rect_domain(starttime, endtime, minlat=evt_domain_range[0], maxlat=evt_domain_range[1], minlon=evt_domain_range[2], maxlon=evt_domain_range[3], 
                                   minmag=evt_min_mag, maxmag=evt_max_mag, min_dep=evt_min_dep, max_dep=evt_max_dep, client_ini=client_ini)
        # save event catalog and plot event map
        save_events_info(events_cat, save_path='.')
    # read event catalog from csv file
    elif evt_domain_type == 3:
        file_format = event_catalog.split('.')[-1]
        if file_format == 'csv':
            events_cat = catcsv2xml(event_catalog, save_path='.')
        else:
            events_cat = read_events(event_catalog)
    else:
        raise SystemExit('Event domain type error!')
    print('Total number of events: %s' % len(events_cat))
    print(events_cat)

    ######  Module 2: Download waveform data by using MassDownloader  ######
    # Define saved data directories
    data_dir = os.getcwd()
    waveform_mseed_dir = os.path.join(data_dir, "waveform_mseed")
    waveform_station_dir = os.path.join(data_dir, "station_inv")
    waveform_sac_dir = os.path.join(data_dir, "waveform_sac")
    # convert station domain type to number
    station_domain_number = {'circ': 1, 'rect': 2, 'globe': 3}
    if sta_domain_type in station_domain_number:
        sta_domain_type = station_domain_number[sta_domain_type]
    else:
        raise SystemExit('Station domain type error!')
    # Download waveform data for each event from all cilents!!!
    for index,event in enumerate(events_cat):
        # Event information.
        event_mag = event.magnitudes[0].mag
        event_time = event.origins[0].time
        event_lat = event.origins[0].latitude
        event_lon = event.origins[0].longitude
        event_dep = event.origins[0].depth
        event_date = '{:04d}-{:02d}-{:02d}-{:02d}-{:02d}-{:02d}'.format(
            event_time.year, event_time.month, event_time.day, event_time.hour,
            event_time.minute, event_time.second)
        # Print the imformation of each event
        print("\n-----------------------------------------")
        print("event,longitude,latitude,magnitude (%s/%s):"%(index+1, len(events_cat)), event_date, event.origins[0].longitude,
              event.origins[0].latitude, event.magnitudes[0].mag)
        # Station data selection for different domain types.
        if sta_domain_type == 1:
            sta_domain = station_dif_domain(sta_domain_type, cen_lat=sta_domain_range[0], cen_lon=sta_domain_range[1], minradius=sta_domain_range[2], maxradius=sta_domain_range[3],
                                        event_lat=event_lat, event_lon=event_lon, event_circ_center=event_circ_center)
        elif sta_domain_type == 2:
            sta_domain = station_dif_domain(sta_domain_type, min_lat=sta_domain_range[0], max_lat=sta_domain_range[1], min_lon=sta_domain_range[2], max_lon=sta_domain_range[3],
                                        event_lat=event_lat, event_lon=event_lon, min_dis_limit=sta_min_dis_limit, max_dis_limit=sta_max_dis_limit, rect_distance_limit=rect_distance_limit)
        elif sta_domain_type == 3:
            sta_domain = station_dif_domain(sta_domain_type)
        # Restrictions for mass_download.
        restrictions = Restrictions(
            # starttime and endtime of waveform data
            starttime=event_time,
            endtime=event_time + wave_len,
            # network and station '*' matches any and can be
            network=array_name,
            station=station_name,
            # If this setting is True, any trace with a gap/overlap will be discarded.
            reject_channels_with_gaps=True,
            # Any trace that is shorter than 95 % of the desired total duration will be discarded.
            minimum_length=0.95,
            sanitize=False,
            minimum_interstation_distance_in_m=0,
            # HH, BH, SH or EH channels.
            channel_priorities=channel,
            # Location codes
            location_priorities=["*"])
        # Define the storage path of waveform data.
        def get_mseed_storage(network, station, location, channel, starttime, endtime):
            # change the format of time
            starttime = starttime.strftime("%Y-%m-%d-%H-%M-%S")
            sac_name = "%s.%s.%s.%s" % (starttime, network, station, channel)
            return os.path.join("%s/%s" % (waveform_mseed_dir, event_date), "%s.mseed" % sac_name)
        try:
            # Download waveform data from all cilents!!!
            print('Downloading waveform data, continue...')
            mdl = MassDownloader(debug=False, configure_logging=False)
            mdl.download(sta_domain, restrictions, mseed_storage=get_mseed_storage, threads_per_client=3,
                         stationxml_storage="%s/{network}.{station}.xml" % (waveform_station_dir))
        except lxml.etree.XMLSyntaxError:
            print('Skipping invalid XML from URL, something have been wrong in one or more stationxml!')
            pass

        ######  Module 3: miniseed2sac: convert miniseed to sac and remove instrument response  ######
        # any mseed file in the event folder
        if not os.path.exists(waveform_mseed_dir) or not os.path.exists(os.path.join(waveform_mseed_dir, event_date)):
            print("\n!!!No miniseed waveform data for event-%s!!!\n" % event_date)
        else:
            mseed_files_exist = any(file.endswith(".mseed") for file in os.listdir(os.path.join(waveform_mseed_dir, event_date)))
            if mseed_files_exist:
                miniseed2sac(waveform_mseed_dir, event_date, waveform_station_dir, waveform_sac_dir, event_lat,
                            event_lon, event_dep, event_mag, delete_mseed=delete_mseed, remove_response=remove_response)
            else:
                print("\n!!!No miniseed waveform data for event-%s!!!\n" % event_date)



