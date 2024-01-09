# seisdownload
![LICENSE](https://img.shields.io/badge/license-MIT-green)
![Author](https://img.shields.io/badge/Author-TianyuCui-blue.svg)


Auto-download mass seismic event data using Obspy.
Functions:star2::
1. Choose different event domain (circural, rectangular and event catalog file domain), station domain (circural, rectangular and global domain).
2. Convert Miniseed format file to SAC format file
3. Automatically remove instrument responce and conventional process (detrend, demean and taper).

***
## Dependencies
#### Tested well:
1. Obspy 1.4.0 
2. Numpy 1.25.0, Pandas 1.4.2
***
## 
#### Details in codes
```Python

    MassDownload_data: download waveform data for different domains.
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
```
***
## Example
see [test/test_download.py](test/test_download.py)

