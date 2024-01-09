# seisdownload
![LICENSE](https://img.shields.io/badge/license-MIT-green)
![Author](https://img.shields.io/badge/Author-TianyuCui-blue.svg)


Auto-download mass seismic event data using Obspy.

<<<<<<< HEAD
Functions:star2::
1. Choose different event domain (circural, rectangular and event catalog file domain), station domain (circural, rectangular and global domain).
=======
Functions:
1. You can choose event domain, station domain and perform the epicentral distance restriction.
>>>>>>> 142f307cb758e362b6da80e05fda991d165d4f7f
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
<<<<<<< HEAD
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
=======
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
>>>>>>> 142f307cb758e362b6da80e05fda991d165d4f7f
```
***
## Example
see test/test_download.py

