#!/usr/bin/env python
# coding: utf-8

import os
import lxml
import obspy
import shutil
import datetime
import numpy as np
import pandas as pd
from obspy.io import sac
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from obspy.clients.fdsn import Client      
from obspy.geodetics import gps2dist_azimuth
from obspy import read, UTCDateTime, read_inventory
from obspy.clients.fdsn.mass_downloader import (
    CircularDomain,
    RectangularDomain,
    GlobalDomain,
    Restrictions,
    MassDownloader,
)

'''
MassDownload_data: download waveform data from all cilents!!!
Download Minisedd data and convert to SAC format.

Author: Tianyu Cui
E-mail: tycuicn@gmail.com
Date: 2023.09.16
'''
def Massdownload_data(array_name, station_name, domain_type, sta_range, evt_range, ref_lat, ref_lon, evt_mag_range, evt_min_dep, 
                      wave_len, channel, startdate, enddate, min_dis=0, max_dis=180, limit_distance=False, delete_mseed=True):
    # Module 1: Get event catalog from IRIS
    evt_minlat = evt_range[0]
    evt_maxlat = evt_range[1]
    evt_minlon = evt_range[2]
    evt_maxlon = evt_range[3]
    # min and max requested magnitudes
    evt_minmag = evt_mag_range[0]
    evt_maxmag = evt_mag_range[1]
    # start and end time of the event catalog
    starttime = UTCDateTime(startdate)
    endtime = UTCDateTime(enddate)
    # search for events from IRIS
    client = Client("IRIS")                        # IRIS Client
    events = client.get_events(starttime=starttime, endtime=endtime, mindepth=evt_min_dep, minlatitude=evt_minlat,
                               maxlatitude=evt_maxlat, minlongitude=evt_minlon, maxlongitude=evt_maxlon,\
                               minmagnitude=evt_minmag, maxmagnitude=evt_maxmag)
    print("Found %s event(s):" % len(events))
    print(events)
    # store data to dataframe
    info_list = ['Origin Time (UTC)', 'Lat [°]', 'Lon [°]', 'depth [m]',
                 'event_type', 'mag', 'magnitude_type', 'creation_info', 'info']
    df = pd.DataFrame(0, index=np.arange(len(events)), columns=info_list)
    for ii in range(0, len(events)):
        df.loc[ii, ("Origin Time (UTC)")] = events[ii].origins[0].time
        df.loc[ii, ('Lat [°]')] = events[ii].origins[0].latitude
        df.loc[ii, ('Lon [°]')] = events[ii].origins[0].longitude
        df.loc[ii, ('depth [m]')] = events[ii].origins[0].depth
        df.loc[ii, ('event_type')] = events[ii].event_type
        df.loc[ii, ('mag')] = events[ii].magnitudes[0].mag
        df.loc[ii, ('magnitude_type')] = events[ii].magnitudes[0].magnitude_type
        df.loc[ii, ('creation_info')] = str(events[ii].origins[0].creation_info)
        df.loc[ii, ('info')] = events[ii].event_descriptions[0].text
    # save to excel
    df.to_excel("events_info.xlsx", sheet_name="events_info")
    # save fig
    events.plot(projection="global", resolution="h", show=False,
                outfile="events_map.png", method='cartopy')
    
    # Module 2: Download waveform data by using MassDownloader
    # Define saved data directories
    data_dir = os.getcwd()
    waveform_mseed_dir = os.path.join(data_dir, "waveform_mseed")
    waveform_station_dir = os.path.join(data_dir, "station_inv")
    waveform_sac_dir = os.path.join(data_dir, "waveform_sac")
    # Download waveform data for each event from all cilents!!!
    for event in events:
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
        print("event,longitude,latitude,magnitude:", event_date, event.origins[0].longitude, 
              event.origins[0].latitude, event.magnitudes[0].mag)
        # Station data selection for different domain types.
        if domain_type == 0:
            # Circular domain around the epicenter.
            domain = CircularDomain(latitude=ref_lat, longitude=ref_lon,
                                    minradius=sta_range[0], maxradius=sta_range[1])
        elif domain_type == 1:
            # Rectangular domain around the epicenter.
            domain = RectangularDomain(minlatitude=sta_range[0], maxlatitude=sta_range[1],
                                    minlongitude=sta_range[2], maxlongitude=sta_range[3])
            if limit_distance:
                # add distance restriction to the Rectangular domain
                domain_restriction = CircularDomain(latitude=event_lat, longitude=event_lon,
                                                    minradius=min_dis, maxradius=max_dis)
        elif domain_type == 2:
            # Global domain.
            domain = GlobalDomain()
        else:
            raise SystemExit('Domain type error!')
        # Waveform data restrictions.
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
            return os.path.join("%s/%s" % (waveform_mseed_dir,event_date), "%s.mseed" % sac_name)
        try:
            # Download waveform data from all cilents!!!
            print('Downloading waveform data, continue...')
            mdl = MassDownloader(debug=False, configure_logging=False) 
            if limit_distance:
                # Add distance restriction based on the Rectangular domain           
                mdl.download((domain and domain_restriction), restrictions, mseed_storage=get_mseed_storage, threads_per_client=3,
                            stationxml_storage="%s/{network}.{station}.xml" % (waveform_station_dir) )
            else:
                mdl.download(domain, restrictions, mseed_storage=get_mseed_storage, threads_per_client=3,
                        stationxml_storage="%s/{network}.{station}.xml" % (waveform_station_dir) )
        except lxml.etree.XMLSyntaxError:
            print('Skipping invalid XML from URL, something have been wrong in one or more stationxml!')
            pass

        # Module 3: miniseed2sac: convert miniseed to sac and remove instrument response
        # any mseed file in the event folder
        mseed_files_exist = any(file.endswith(".mseed") for file in os.listdir(
            os.path.join(waveform_mseed_dir, event_date)))
        if mseed_files_exist:
            miniseed2sac(waveform_mseed_dir, event_date, waveform_station_dir, waveform_sac_dir, event_lat, 
                         event_lon, event_dep, event_mag, delete_mseed=delete_mseed)
        else:
            print("\n!!!No miniseed waveform data for event-%s!!!\n" % event_date)


'''
Add SAC header values to trace
'''
def mseed_to_sac_header(trace,header_info):
    # set SAC header values
    sacz = sac.SACTrace.from_obspy_trace(trace)
    sacz.stlo = header_info['sta_lon']             # station longitude
    sacz.stla = header_info['sta_lat']             # station latitude
    sacz.stel = header_info['sta_ele']             # station elevation
    sacz.kstnm = header_info['sta']                # station name
    sacz.kcmpnm = header_info['cha']               # channel code
    sacz.knetwk = header_info['net']               # network code
    sacz.khole = header_info['loc']                # location code
    sacz.mag = header_info['eve_mag']              # event magnitude
    sacz.evlo = header_info['eve_lon']             # event longitude
    sacz.evla = header_info['eve_lat']             # event latitude
    sacz.evdp = header_info['eve_dep']/1000        # event depth(km)
    sacz.az = header_info['azi']                   # azimuth
    sacz.baz = header_info['baz']                  # back azimuth
    sacz.dist = header_info['dist']/1000           # distance in kilometers
    sacz.gcarc = header_info['dist']/111190        # distance in degree
    sacz.delta = header_info['delta']              # delta
    sacz.o = 0                                     # set event origin time as reference time
    return sacz


'''
miniseed2sac: convert miniseed to sac and remove instrument response
'''


def miniseed2sac(waveform_mseed, event_date, station_dir, waveform_sac, eve_lat, eve_lon, eve_dep, eve_mag, rotate_sac=False, delete_mseed=True):
    mseed_dir = os.path.join(waveform_mseed, event_date)
    sac_dir = os.path.join(waveform_sac, event_date)
    if not os.path.exists(waveform_sac):
        os.mkdir(waveform_sac)
    if os.path.isdir(sac_dir):
        shutil.rmtree(sac_dir)
    os.mkdir(sac_dir)
    st = read("%s/*.mseed" % mseed_dir)
    for tr in st:
        if np.isnan(np.max(tr.data)) or np.isinf(np.max(tr.data)):
            st.remove(tr)
        net = tr.stats.network
        sta = tr.stats.station
        cha = tr.stats.channel
        loc = tr.stats.location
        station_inv = os.path.join(station_dir, '%s.%s.xml'%(net, sta))
        # get corresponding SAC header values from StationXML
        if not os.path.exists(station_inv):
            current_date = datetime.date.today()
            try:    
                Client('IRIS').get_stations(starttime=UTCDateTime('1990-01-01'), endtime=UTCDateTime(current_date),
                                    network=net, station=sta, channel=cha, location=loc, level='response',
                                    filename=station_inv, format='xml')
            except:
                pass
        if os.path.exists(station_inv):
            remove_instrument = True
            tr_inv = read_inventory(station_inv)
            coordinates = tr_inv.get_coordinates(net + '.' + sta + '.' + loc + '.' + cha)
            sta_lon = coordinates['longitude']
            sta_lat = coordinates['latitude']
            sta_ele = coordinates['elevation']
            # calculate the distance, azimuth and back azimuth
            (dist, azi, baz) = gps2dist_azimuth(eve_lat, eve_lon, sta_lat, sta_lon)
            # SAC header information
            header_info = {'sta_lon': sta_lon, 'sta_lat': sta_lat, 'sta_ele': sta_ele, 'sta': sta, 'cha': cha,
                           'net': net, 'loc': loc, 'eve_mag': eve_mag, 'eve_lon': eve_lon, 'eve_lat': eve_lat, 'eve_dep': eve_dep,
                           'azi': azi, 'baz': baz, 'dist': dist, 'delta': tr.stats.delta}
            # Remove instrument response
            # Notice: instrument response removal by obspy differs with that by SAC software due to water_level !!!
            tr.detrend("demean")
            tr.detrend("linear")
            pre_filt = [0.001, 0.002, 25, 30]
            try:
                # displacement, output unit is meters
                tr.remove_response(inventory=tr_inv, water_level=60, taper=True, 
                                    taper_fraction=0.00001, pre_filt=pre_filt, output="DISP")
                # tr.data = tr.data * 1e9  # convert to nm
            except:
                remove_instrument = False
                print("!!!%s/%s.%s.%s.%s.sac remove response failed!!!" % (sac_dir, event_date, net, sta, cha))
                continue
            # rotate to ZNE, optional
            if rotate_sac:
                tr.rotate(method="->ZNE", inventory=tr_inv)
            sacz = mseed_to_sac_header(tr, header_info)
            if remove_instrument:
                sacz.write("%s/%s.%s.%s.%s.sac" % (sac_dir, event_date, net, sta, cha))
                # Delete miniseed files if miniseed convert to sac successfully
                if delete_mseed and os.path.exists('%s/%s.%s.%s.%s.sac' % (sac_dir, event_date, net, sta, cha)):
                    os.system('rm %s/*%s.%s.%s.mseed' % (mseed_dir, net, sta, cha))
            else:
                unremove_file = os.path.join(sac_dir, 'unremove_sac')
                if not os.path.exists(unremove_file):
                    os.mkdir(unremove_file)
                sacz.write("%s/%s.%s.%s.%s.sac" % (unremove_file, event_date, net, sta, cha))
                # Delete miniseed files if miniseed convert to sac successfully
                if delete_mseed and os.path.exists('%s/%s.%s.%s.%s.sac' % (unremove_file, event_date, net, sta, cha)):
                    os.system('rm %s/*%s.%s.%s.mseed' % (mseed_dir, net, sta, cha))


if __name__ == '__main__':
    '''
    Author: Tianyu Cui
    Date: 2023.09.16

    arrayname: "IU" or "II" or "TA" or "TW" or "IC" or "IU,II,TA,TW,IC" or "*"
    station_name: "ANMO" or "TA01" or "ANMO,TA01" or "*"
    channel: channels (default: ["BHZ", "HHZ", "SHZ", "EHZ"])
    sta_range: 
        domain type:1 (RectangularDomain) sta_range = [sta_lat_min, sta_lat_max, sta_lon_min, sta_lon_max] in degree
                       if limit_distance=True, add distance restriction to the Rectangular domain
                      (RestrictionDomain) [min_dis, max_dis] in degree 
        domain type:2 (CircularDomain) sta_range = [minradius, maxradius] in degree 
                                       mid points: [ref_lat, ref_lon] in degree
        domain type:3 (GlobalDomain) []
    evt_range: [evt_lat_min, evt_lat_max, evt_lon_min, evt_lon_max] in degree
    evt_mag_range: [evt_mag_min, evt_mag_max]
    evt_min_dep: min event depth in km
    wave_len: downloaded waveform length in seconds
    startdate: earthquake catalog start date
    enddate: earthquake catalog end date
    limit_distance: if True, add distance restriction to the Rectangular domain (default: False)
                    min_dis: min distance in degree (default: 0)
                    max_dis: max distance in degree (default: 180)
    delete_mseed: if True, delete corresponding miniseed data if miniseed convert to sac successfully (default: True)
    '''
    Massdownload_data(array_name="*", station_name="*", domain_type=1, sta_range=[0, 60, 40, 180], evt_range=[-10, 60, 40, 180],
                      ref_lat=0, ref_lon=0, evt_mag_range=[5.5, 10], evt_min_dep=50, channel=["BHZ", "HHZ", "SHZ", "EHZ"], wave_len=1800,
                      startdate="2015-01-01 00:00:00", enddate="2015-01-10 21:59:59", max_dis=15, limit_distance=True, delete_mseed=True)

