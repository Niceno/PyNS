# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 14:12:35 2017

@author: kerstin.cramer
"""

# Standard Python modules
from pyns.standard import *

import numpy as np

#--------------------------------------------------------------------------
def p_v_sat(t):
#--------------------------------------------------------------------------
  
  t = t + 273.15 # convert to K
  
  if (isinstance(t,int) or isinstance(t,float)):
    if   (t >334 and t<=363.15):
      p_v= 1E+5 * np.power(10,5.0768  - 1659.793/(t-45.854))
    elif (t >304 and t<=334):
      p_v= 1E+5 * np.power(10,5.20389 - 1733.926/(t-39.485))
    elif (t >273 and t<=304):
      p_v= 1E+5 * np.power(10,5.40221 - 1838.675/(t-31.737))
    else:
      p_v= float('nan')
 
  else:
    dim=t.ndim
    p_v=zeros(t.shape)
  
    if dim==1:
      for ii in range(0,len(t)):
        if   (t[ii] >334 and t[ii]<=363.15):
          p_v[ii]= 1E+5 * np.power(10,5.0768  - 1659.793/(t[ii]-45.854))
        elif (t[ii] >304 and t[ii]<=334):
          p_v[ii]= 1E+5 * np.power(10,5.20389 - 1733.926/(t[ii]-39.485))
        elif (t[ii] >273 and t[ii]<=304):
          p_v[ii]= 1E+5 * np.power(10,5.40221 - 1838.675/(t[ii]-31.737))
          
    elif dim==2:
      for ii in range(0,t.shape[0]):
        for jj in range(0,t.shape[1]):
          if   (t[ii,jj] >334 and t[ii,jj]<=363.15):
            p_v[ii,jj]= 1E+5 * np.power(10,5.0768  - 1659.793/(t[ii,jj]-45.854))
          elif (t[ii,jj] >304 and t[ii,jj]<=334):
            p_v[ii,jj]= 1E+5 * np.power(10,5.20389 - 1733.926/(t[ii,jj]-39.485))
          elif (t[ii,jj] >273 and t[ii,jj]<=304):
            p_v[ii,jj]= 1E+5 * np.power(10,5.40221 - 1838.675/(t[ii,jj]-31.737))
            
    elif dim==3:
      for ii in range(0,t.shape[0]):
        for jj in range(0,t.shape[1]):
          for kk in range(0,t.shape[2]):
            if   (t[ii,jj,kk] >334 and t[ii,jj,kk]<=363.15):
              p_v[ii,jj,kk]= 1E+5 * np.power(10,5.0768  - 1659.793/(t[ii,jj,kk]-45.854))
            elif (t[ii,jj,kk] >304 and t[ii,jj,kk]<=334):
              p_v[ii,jj,kk]= 1E+5 * np.power(10,5.20389 - 1733.926/(t[ii,jj,kk]-39.485))
            elif (t[ii,jj,kk] >273 and t[ii,jj,kk]<=304):
              p_v[ii,jj,kk]= 1E+5 * np.power(10,5.40221 - 1838.675/(t[ii,jj,kk]-31.737))
  
    p_v[p_v==0]= float('nan')
  
  return p_v # end of function
  
#--------------------------------------------------------------------------
def p_v_sat_salt(t,a,M_salt):
#--------------------------------------------------------------------------
  
  p_v=p_v_sat(t)
  
  p_v = p_v *np.exp(-1.17444*np.power(a/M_salt,0.5) / (1+np.power(a/M_salt,0.5)))
  
  return p_v # end of function
 
#--------------------------------------------------------------------------
def t_sat(p_v):
#--------------------------------------------------------------------------
  
  if (isinstance(p_v,int) or isinstance(p_v,float)):
    if (p_v >2.0725e+04 and p_v <=7.010428e+04):
      t= 45.854 + 1659.793 /(5.0768  - np.log10(p_v * 1E-5))
    elif (p_v >4.4556e+03 and p_v<=2.0725e+04):
      t= 39.485 + 1733.926 /(5.20389 - np.log10(p_v * 1E-5))
    elif (p_v >604.1854 and p_v<=4.4556e+03):
      t = 31.737 + 1838.675 /(5.40221 - np.log10(p_v * 1E-5))
    else:
      t= float('nan')
 
  else:
    dim=p_v.ndim
    t=zeros(p_v.shape)
  
    if dim==1:
      for ii in range(0,len(p_v)):
        if   (p_v[ii] >2.0725e+04 and p_v[ii]<=7.010428e+04):
          t[ii]= 45.854 + 1659.793 /(5.0768  - np.log10(p_v[ii] * 1E-5))		
        elif (p_v[ii] >4.4556e+03 and p_v[ii]<=2.0725e+04):
          t[ii]= 39.485 + 1733.926 /(5.20389 - np.log10(p_v[ii] * 1E-5))	
        elif (p_v[ii] >604.1854 and p_v[ii]<=4.4556e+03):
          t[ii]= 31.737 + 1838.675 /(5.40221 - np.log10(p_v[ii] * 1E-5))	
          
    elif dim==2:
      for ii in range(0,t.shape[0]):
        for jj in range(0,t.shape[1]):
          if   (p_v[ii,jj] >2.0725e+04 and p_v[ii,jj]<=7.010428e+04):
            t[ii,jj]= 45.854 + 1659.793 /(5.0768  - np.log10(p_v[ii,jj] * 1E-5))
          elif (p_v[ii,jj] >4.4556e+03 and p_v[ii,jj]<=2.0725e+04):
            t[ii,jj]= 39.485 + 1733.926 /(5.20389 - np.log10(p_v[ii,jj] * 1E-5))
          elif (p_v[ii,jj] >604.1854 and p_v[ii,jj]<=4.4556e+03):
            t[ii,jj]= 31.737 + 1838.675 /(5.40221 - np.log10(p_v[ii,jj] * 1E-5))
            
    elif dim==3:
      for ii in range(0,p_v.shape[0]):
        for jj in range(0,p_v.shape[1]):
          for kk in range(0,p_v.shape[2]):
            if   (p_v[ii,jj,kk] >2.0725e+04 and p_v[ii,jj,kk]<=7.010428e+04):
              t[ii,jj,kk]= 45.854 + 1659.793 /(5.0768  - np.log10(p_v[ii,jj,kk] * 1E-5))
            elif (p_v[ii,jj,kk] >4.4556e+03 and p_v[ii,jj,kk]<=2.0725e+04):
              t[ii,jj,kk]= 39.485 + 1733.926 /(5.20389 - np.log10(p_v[ii,jj,kk] * 1E-5))
            elif (p_v[ii,jj,kk] >604.1854 and p_v[ii,jj,kk]<=4.4556e+03):
              t[ii,jj,kk]= 31.737 + 1838.675 /(5.40221 - np.log10(p_v[ii,jj,kk] * 1E-5))
  
    t[t==0]= float('nan')
  
  t = t - 273.15    # convert to °C
    
  return t # end of function