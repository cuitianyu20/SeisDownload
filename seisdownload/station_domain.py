# This file is part of seisdownload, a Python package for seismic event download.
# Copyright (C) 2024 Tianyu Cui (tycuicn@gmail.com)

'''
Define station domain for different domain types
domain_type: 0 (CircularDomain) or 1 (RectangularDomain) or 2 (GlobalDomain)
rectangular domain: [minlatitude, maxlatitude, minlongitude, maxlongitude] in degree
rectangular domain with distance restriction: [minlatitude, maxlatitude, minlongitude, maxlongitude, min_dis_limit, max_dis_limit] in degree
circular domain: 
    event center: [event_lat, event_lon, minradius, maxradius] in degree
    reference center: [ref_lat, ref_lon, minradius, maxradius] in degree
'''


from obspy.clients.fdsn.mass_downloader import (
    CircularDomain,
    RectangularDomain,
    GlobalDomain,
)


def station_dif_domain(domain_type, cen_lat=None, cen_lon=None, minradius=None, maxradius=None, event_lat=None, event_lon=None,
                   min_lat=None, max_lat=None, min_lon=None, max_lon=None, min_dis_limit=0, max_dis_limit=180, 
                   event_circ_center=True, rect_distance_limit=False):
    # Circular domain around the epicenter.
    if domain_type == 1 and event_circ_center:
        if minradius is not None and maxradius is not None:
            domain = CircularDomain(latitude=event_lat, longitude=event_lon,
                                    minradius=minradius, maxradius=maxradius)
        else:
            raise SystemExit('Please input minradius and maxradius! (Circular station domain)')
    # Circular domain around a sepcified point (ref_circ_lat, ref_circ_lon).
    elif domain_type == 1 and not event_circ_center:
        value_list = [cen_lat, cen_lon, minradius, maxradius]
        all_not_none = all(value is not None for value in value_list)
        if all_not_none:
            domain = CircularDomain(latitude=cen_lat, longitude=cen_lon,
                                    minradius=minradius, maxradius=maxradius)
        else:
            raise SystemExit('Please input cen_lat, cen_lon, minradius and maxradius! (Circular station domain)')
    # Rectangular domain around the epicenter. 
    elif domain_type == 2 and not rect_distance_limit:
        value_list = [min_lat, max_lat, min_lon, max_lon]
        all_not_none = all(value is not None for value in value_list)
        if all_not_none:
            domain = RectangularDomain(minlatitude=min_lat, maxlatitude=max_lat,
                                        minlongitude=min_lon, maxlongitude=max_lon)
        else:
            raise SystemExit('Please input min_lat, max_lat, min_lon and max_lon! (Rectangular station domain)')
    # Rectangular domain around the epicenter with epicentral distance restriction.
    elif domain_type == 2 and rect_distance_limit:
        value_list = [min_lat, max_lat, min_lon, max_lon, event_lat, event_lon, min_dis_limit, max_dis_limit]
        all_not_none = all(value is not None for value in value_list)
        if all_not_none:
            domain = RectangularDomain(minlatitude=min_lat, maxlatitude=max_lat,
                                        minlongitude=min_lon, maxlongitude=max_lon)
            # add distance restriction to the Rectangular domain
            domain_restriction = CircularDomain(latitude=event_lat, longitude=event_lon,
                                                minradius=min_dis_limit, maxradius=max_dis_limit)
            domain = domain and domain_restriction
        else:
            raise SystemExit('Please input min_lat, max_lat, min_lon, max_lon, event_lat, event_lon, min_dis_limit and max_dis_limit! (Rectangular station domain)')
    # Global domain.
    elif domain_type == 3:
        domain = GlobalDomain()
    else:
        raise SystemExit('Station domain type error!')
    return domain