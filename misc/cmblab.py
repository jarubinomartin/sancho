#!/usr/bin/env python

"""
14/Jun/2025   

@author: jalberto

CMBLAB coordinates and locations for the different CMB experiments

REQUIRES:
astropy

HISTORY:
* 14/Jun/2025 - original version. 

"""



import numpy as np
import astropy as ap
import astropy.coordinates
import astropy.units as u


# Locations of CMB telescopes at Teide Observatory
def cmblab_locations():
    return {
        "QT1": {
            "location": ap.coordinates.EarthLocation(lat=28.300340*u.deg, lon=-16.5101*u.deg, height=2395*u.m),
            "description": "First QUIJOTE telescope"
        },
        "QT2": {
            "location": ap.coordinates.EarthLocation(lat=28.300275*u.deg, lon=-16.5101*u.deg, height=2395*u.m),
            "description": "Second QUIJOTE telescope"
        },
        "Tenerife": {
            "location": ap.coordinates.EarthLocation(lat=28.2925*u.deg, lon=-16.5*u.deg, height=2400*u.m),
            "description": "Tenerife Experiment"
        },
        "TMS": {
            "location": ap.coordinates.EarthLocation(lat=28.300571*u.deg,lon=-16.510272*u.deg,height=2395*u.m),
            "description": "Tenerife Microwave Spectrometer"
        },
        "VSA": {
            "location": ap.coordinates.EarthLocation(lat=28.300643*u.deg,lon=-16.510279*u.deg,height=2395*u.m),
            "description": "Very Small Array"
        }
    }

