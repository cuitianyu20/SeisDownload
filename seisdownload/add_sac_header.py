# This file is part of seisdownload, a Python package for seismic event download.
# Copyright (C) 2024 Tianyu Cui (tycuicn@gmail.com)

from obspy.io import sac

'''
Add SAC header values to trace
'''
def mseed_to_sac_header(trace, header_info):
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
    # set event origin time as reference time
    sacz.o = 0
    return sacz