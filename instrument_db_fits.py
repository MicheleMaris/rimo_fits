import sys

__DESCRIPTION__="""
The instrument_db_fits library is used to upload LFI and HFI instrument_db.
"""

__VERSION__="""Version 1 - 2014 July 31 - By M. Maris, INAF / Trieste Astronomical Observatory and Planck/LFI DPC"""

from rimo_fits import *
from csv_table import csv_table

class instrument_db_fits(csv_table) :
   def __init__(self,fitsfilename,ihdu=1,autoComplete=True) :
      """
      ihdu=index of hdu to be readed (default 1)
      herits from csv_table
      autoComplete=True (default) introduces extra columns"""
      import pyfits
      # creates an empty csv_table
      csv_table.__init__(self,None)
      self.__info__['primary_header']=None
      self.__info__['header']=None
      if fitsfilename == None : return
      try :
         o=pyfits.open(fitsfilename)
      except :
         return
      self.__fitsfilename__='hdu='+str(ihdu)+'@file='+fitsfilename
      self.__info__['primary_header']=o[0].header
      self.__info__['header']=o[ihdu].header
      self.get_from_hdu(o[ihdu])
      o.close()
      if autoComplete : self.__complete()
   #def get_from_hdu(self,hdu) :
      #"gets from hdu the instrument_db"
      #for it in range(1,hdu.header['TFIELDS']+1) :
         #n=hdu.header['TTYPE%d'%it]
         #u=hdu.header['TUNIT%d'%it]
         #f=hdu.header['TFORM%d'%it]
         #self.newcolumn(n.capitalize(),hdu.data.field(n),unit=u,shape=f)
   def __complete(self) :
      "method to allow introduction of autocompletion columns"
      import numpy as np
      inst = []
      fh = []
      arm = []
      for k in self.Radiometer :
         inst.append(k[0:3].upper())
         if k[0:3].upper() == 'LFI' :
            fh.append(int(k[3:5]))
            arm.append(k[5:])
         else :
            fh.append(-1)
            arm.append('U')
      fh=np.array(fh)
      ch=70*(18<=fh)*(fh<= 23)
      ch+=44*(24<=fh)*(fh<= 26)
      ch+=30*(27<=fh)*(fh<= 28)
      self.newcolumn('instrument',np.array(inst),unit='',shape='1d',description='instrument')
      self.newcolumn('fh',fh,unit='',shape='1d',description='feed-horn')
      self.newcolumn('arm',np.array(arm),unit='',shape='1d',description='arm')
      self.newcolumn('ch',np.array(ch),unit='',shape='1d',description='frequency channel')
   def header(self,primary=False) :
      "returns the fits header"
      if primary :
         return self.__info__['primary_header']
      return self.__info__['header']
