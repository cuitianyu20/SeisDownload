# This file is part of seisdownload, a Python package for seismic event download.
# Copyright (C) 2024 Tianyu Cui (tycuicn@gmail.com)


import os
import shutil
import datetime
import warnings
import numpy as np
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth
from obspy import read, UTCDateTime, read_inventory
from .add_sac_header import mseed_to_sac_header

'''
miniseed2sac: convert miniseed to sac and remove instrument response
'''
def miniseed2sac(waveform_mseed, event_date, station_dir, waveform_sac, eve_lat, eve_lon, eve_dep, eve_mag, 
                 rotate_sac=False, delete_mseed=True, remove_response=False):
    mseed_dir = os.path.join(waveform_mseed, event_date)
    sac_dir = os.path.join(waveform_sac, event_date)
    if not os.path.exists(waveform_sac):
        os.mkdir(waveform_sac)
    if os.path.isdir(sac_dir):
        shutil.rmtree(sac_dir)
    os.mkdir(sac_dir)
    try:
        st = read("%s/*.mseed" % mseed_dir)
        for tr in st:
            if np.isnan(np.max(tr.data)) or np.isinf(np.max(tr.data)):
                st.remove(tr)
            net = tr.stats.network
            sta = tr.stats.station
            cha = tr.stats.channel
            loc = tr.stats.location
            station_inv = os.path.join(station_dir, '%s.%s.xml' % (net, sta))
            # get corresponding SAC header values from StationXML
            if not os.path.exists(station_inv):
                current_date = datetime.date.today()
                try:
                    print("Downloading station inventory (net:%s, sta:%s)..." % (net,sta))
                    Client('IRIS').get_stations(starttime=UTCDateTime('1990-01-01'), endtime=UTCDateTime(current_date),
                                                network=net, station=sta, channel=cha, location=loc, level='response',
                                                filename=station_inv, format='xml')
                    if os.path.exists(station_inv):
                        print("Download station inventory (net:%s, sta:%s) successfully!" % (net,sta))
                except:
                    pass
            if os.path.exists(station_inv):
                try:
                    tr_inv = read_inventory(station_inv)
                    coordinates = tr_inv.get_coordinates(net + '.' + sta + '.' + loc + '.' + cha)
                    sta_lon = coordinates['longitude']
                    sta_lat = coordinates['latitude']
                    sta_ele = coordinates['elevation']
                    # calculate the distance, azimuth and back azimuth
                    warnings.filterwarnings("ignore", category=DeprecationWarning)
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
                    if remove_response:
                        # displacement, output unit is nm
                        tr.remove_response(inventory=tr_inv, water_level=60, taper=True,
                                        taper_fraction=0.00001, pre_filt=pre_filt, output="DISP")
                        tr.data = tr.data * 1e9  # convert to nm
                    # write SAC header values to trace
                    sacz = mseed_to_sac_header(tr, header_info)
                    # rotate to ZNE, optional
                    if rotate_sac:
                        tr.rotate(method="->ZNE", inventory=tr_inv)
                    # write SAC trace to file
                    sacz.write("%s/%s.%s.%s.%s.sac" % (sac_dir, event_date, net, sta, cha))
                    # Delete miniseed files if miniseed convert to sac successfully
                    if delete_mseed and os.path.exists('%s/%s.%s.%s.%s.sac' % (sac_dir, event_date, net, sta, cha)):
                        os.system('rm %s/*%s.%s.%s.mseed' %(mseed_dir, net, sta, cha))
                except Exception as e:
                    print("!!!%s/%s.%s.%s.%s.sac read inventory failed!!!" % (sac_dir, event_date, net, sta, cha))
                    continue
        try:
            os.rmdir(mseed_dir)
            print(f"The folder:'{mseed_dir}' has been deleted successfully!")
        except OSError as e:
            print("The folder %s : %s" % (mseed_dir, e.strerror))
    except:
        pass
