# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 16:21:35 2024

@author: z5379427
"""

import numpy as np


def ffdi(df,t,rh,u):    
    ffdi = 2*np.exp(-0.45+0.987*np.log(df)-0.0345*rh+0.0338*t+0.0234*u)
    return ffdi

def u_ffdi(ffdi,df,t,rh):    
    u_ffdi = (np.log(ffdi/2) + 0.45 - 0.987*np.log(df) + 0.0345*rh - 0.0338*t)/0.0234    
    return u_ffdi

def u_log_abl(utau,k,z,z0):    
    u_log_abl = utau/k * np.log(z/z0)    
    return u_log_abl

def utau_log_abl(u,k,z,z0):    
    utau_log_abl = (k*u) / np.log(z/z0)    
    return utau_log_abl




kappa = 0.41
z0 = 0.003
can_h = 17

d_factor = 10;
temp = 35;
rel_hum = 10;
FFDI = 80;

u_10_ms = u_ffdi(FFDI, d_factor, temp, rel_hum)*1000/3600 # defined as km/h at 10 m
u_tau = utau_log_abl(u_10_ms, kappa, 10, z0) # based on z = 10 m
u_4h = u_log_abl(u_tau, kappa, 4*can_h, z0)
