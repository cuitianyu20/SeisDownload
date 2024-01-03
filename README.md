# SeisDownload
![LICENSE](https://img.shields.io/badge/license-MIT-green)
![Author](https://img.shields.io/badge/Author-TianyuCui-blue.svg)


Auto-download mass seismic event data using Obspy.

Functions:
1. You can choose event domain, station domain and perform the epicentral distance restriction.
2. Convert Miniseed format file to SAC format file
3. Automatically remove instrument responce and conventional process (detrend, demean and taper).

***
## Dependencies
#### Tested well:
1. Obspy1.4.0 
2. Numpy 1.25.0,Pandas 1.4.2
***
## Input Parameters
#### Details in codes
```Python
    '''
    Author: Tianyu Cui
    Date: 2023.09.16
    arrayname: "IU" or "II" or "TA" or "TW" or "IC" or "IU,II,TA,TW,IC" or "*"
    station_name: "ANMO" or "TA01" or "ANMO,TA01" or "*"
    channel: channels (default: ["BHZ", "HHZ", "SHZ", "EHZ"])
    sta_range: 
        domain type:0 (CircularDomain) sta_range = [minradius, maxradius] in degree
                        if event_center=True, use event center as reference point
                        mid points: [event_lat, event_lon] in degree
                        if event_center=False, use reference point as reference point 
                        mid points: [ref_lat, ref_lon] in degree
        domain type:1 (RectangularDomain) sta_range = [sta_lat_min, sta_lat_max, sta_lon_min, sta_lon_max] in degree
                       if limit_distance=True, add distance restriction to the Rectangular domain
                      (RestrictionDomain) [min_dis, max_dis] in degree 
        domain type:2 (GlobalDomain) []
    evt_range: [evt_lat_min, evt_lat_max, evt_lon_min, evt_lon_max] in degree
    evt_mag_range: [evt_mag_min, evt_mag_max]
    evt_min_dep: min event depth in km
    ref_lat: reference latitude in degree
    ref_lon: reference longitude in degree
    wave_len: downloaded waveform length in seconds
    startdate: earthquake catalog start date
    enddate: earthquake catalog end date
    min_dis: min distance in degree (default: 0)
    max_dis: max distance in degree (default: 180)
    event_center: if True, use event center as reference point (default: False)
    limit_distance: if True, add distance restriction to the Rectangular domain (default: False)
                    min_dis: min distance in degree (default: 0)
                    max_dis: max distance in degree (default: 180)
    remove_response: if True, remove instrument response (default: False)
    delete_mseed: if True, delete corresponding miniseed data if miniseed convert to sac successfully (default: True)
    '''
```
***
## Demo
##### Successfully Download One Earthquake Data

```
* Network: all.
* Statin: all.
* Channel: 'BHZ' # Z component
* Recorded waveform length: 1800s
* Event magnitude: greater than M5.
* Event depth: greater than 50 km.
* Time limitation: from 2015-01-01T00:00:00 to 2015-01-10T21:59:59.
* Domain: Rectangular domain and add distance restriction (0°-15°)
```
