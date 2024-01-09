# This file is part of seisdownload, a Python package for seismic event download.
# Copyright (C) 2024 Tianyu Cui (tycuicn@gmail.com)

import numpy as np
import pandas as pd
from obspy import UTCDateTime
from obspy.core.event import Catalog, Event, Origin, Magnitude


'''
Convert event catalog from csv to xml format and return obspy event catalog
'''
def catcsv2xml(cat_csv, save_path='.'):
    # read csv file, header example: time, lat, lon, depth, mag
    event_data = pd.read_csv(cat_csv, header=0)
    # convert to obspy catalog
    event_cat = Catalog()
    for index, row in event_data.iterrows():
        event = Event()
        event.origins.append(Origin(latitude=row['lat'], longitude=row['lon'], depth=row['depth'], time=UTCDateTime(row['event'])))
        event.magnitudes.append(Magnitude(mag=row['mag']))
        event_cat.append(event)
    # save to xml file
    event_cat.write("%s/events.xml"%save_path, format="QUAKEML")
    return event_cat


'''
Save event catalog to excel and plot event map
'''
def save_events_info(events, save_path='.'):
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
    df.to_excel("%s/events_info.xlsx"%save_path, sheet_name="events_info")
    # save fig
    events.plot(projection="global", resolution="h", show=False,
                outfile="%s/events_map.png"%save_path, method='cartopy')
