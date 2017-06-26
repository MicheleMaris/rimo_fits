import sys

__DESCRIPTION__="""
The rimofits library is used to upload LFI and HFI RIMO.
From version 1.3 includes Fourier transform of bandpasses and fourier filtering for high frequencies

Added possibility to have a table of statistics for LFI_RIMO
"""

__VERSION__="""Version 1.8 - 29 May 2012 - 2012 Jun 25 - 2012 Nov 21 - 2013 Jan 24 - 2013 Jan 25 - 2013 Apr 18 - 2014 Mar 20 - 2014 May 12 - 2014 Oct 20 - 2017 Jun 22 - By M. Maris, INAF / Trieste Astronomical Observatory and Planck/LFI DPC
"""

#class SpectraEnergyDistributionCollection :
   #class sed :
      #def __init__(self,name,func,fmin,fmax,underflow,overflow) :
         #self.name=name
         #self.fmin=fmin
         #self.fmax=fmax
         #self.underflow=underflow
         #self.overflow=overflow
         #self.func = func
      #def __calc__(self,Freq) :
         #return self.func(Freq)
   #def __init__(self) : pass
   #def etaDeltaT(self,FreqGHz) :
      #from blackbody_array import CMB
      #out=CMB.etaDeltaT(self,FreqGHz)
      #outN=((outN[1:]+outN[:-1])*(FreqGHz[1:]-FreqGHz[:-1])*0.5).sum()
      #return out/outN
   #def 

def Detector2Arm(arm) :
   "Convert detector code into all the possible codes for detector or arm"
   _arm=str(arm).strip().lower()
   if _arm=='s' or _arm=='x' or _arm=='1' : return ['S','X',1]
   if _arm=='m' or _arm=='y' or _arm=='0' : return ['M','Y',0]
   return None

class BandPassStatistics :
   def __init__(self,BP,FcentGHz=None) :
      """Computes a number of usefull band pass statistics"""
      from blackbody import BlackBody,CMB
      self.Fcent = BP.Calc_Fcent() if FcentGHz == None else FcentGHz
      self.MAX = BP.Trans.max()
      self.Bcmb=BlackBody.bbn_cgs(self.Fcent*1e9,CMB.Tcmb)/1e-23/1e6
      self.Brj1=BlackBody.bbn_rj_cgs(self.Fcent*1e9,1)/1e-23/1e6
      self.dBdT_cmb=BlackBody.bbn_diff(self.Fcent,CMB.Tcmb)['bbn_diff']
   def __getitem__(self,k) :
      return self.__dict__[k]
   def Tb(self,Snu,MJySr=False) :
      "returns brightness temperature"
      from blackbody import BlackBody
      return BlackBody.Tb(self.Fcent*1e9,Snu,MJySr=MJySr)
   def keys(self) :
      return self.__dict__.keys()

class BandPass :
   """class to handle a BandPass"""
   def __init__(self) :
      self.Wavenumber=None
      self.Trans=None
      self.Error=None
      self.WeightedBy=None
   def __len__(self) :
      try :
         return len(self.Wavenumber)
      except :
         return 0
   def Statistics(self) :
      return BandPassStatistics(self)
   def Calc_Fcent(self) :
      """computes the central frequency for the given bandpass"""
      x=self.Trans*self.FreqGHz
      y=self.Trans
      r=0.5*(x[1:]+x[:-1]).sum()
      r1=0.5*(y[1:]+y[:-1]).sum()
      return r/r1
   def keys(self,onlyArray=False) :
      if onlyArray :
         return ['Wavenumber','Trans','Error','Flag','FreqGHz']
      return self.__dict__.keys()
   def __getitem__(self,k) :
      return self.__dict__[k]
   def __setitem__(self,k,that) :
      self.__dict__[k]=that
   def __str__(self) :
      line=[]
      for k in ['EXTNAME','Chain','step'] :
         try :
            line.append(k+' = '+str(self[k]))
         except : 
            line.append(k+' = unknown')
            pass
      line.append('ELEMENTS = '+str(len(self)))
      return "\n".join(line)
   def weightedBandPass(self,name,SED) :
      """returns a bandpass weigthed by a specific SED
         SED is a function returning a weight at a given FreqGHz
      """
      out=self.copy()
      out.Trans*=SED(out.FreqGHz)
      out['WeightedBy']=name
      norm=((out.Trans[1:]+out.Trans[:-1])*(out.FreqGHz[1:]-out.FreqGHz[:-1])/2.).sum()
      out.Trans*=1./norm
      return out
   def copy(self) :
      import copy
      return copy.deepcopy(self)
   def save(self,pickle_file) :
      "save as pickle file"
      import pickle
      if type(pickle_file)==type('') :
         self.filename=pickle_file
         try :
            pickle.dump(self.__dict__,open(pickle_file,'w'))
         except :
            return False
      else :
         try :
            pickle.dump(self.__dict__,pickle_file)
         except :
            return False
      return True
   def load(self,pickle_file) :
      "load from pickle file"
      import pickle
      if type(pickle_file)==type('') :
         self.filename=pickle_file
         try :
            self.__dict__=pickle.load(open(pickle_file,'r'))
         except :
            return False
      else :
         try :
            self.__dict__=pickle.load(pickle_file)
         except :
            return False
      return True
   def fft(self,fftFilt,forceToZero=True,renorm=True,asBP=False) :
      """returns the fft of the bandpass, fftFilt=None, no filtering, fftFilt<0 highpass filter, fftFilt>0, lowpass filter
       negative samples are forced to zero
      """
      import numpy.fft as FFT
      import numpy as np
      step=(self['FreqGHz'][1:]-self['FreqGHz'][0:-1]).mean()
      fbaseline=FFT.fftfreq(len(self['FreqGHz']))/step
      mT=self['Trans'].mean()
      ft=FFT.fft(self['Trans']-mT)
      if fftFilt == None :
         return {'ft-zero':mT,'ft-step':step,'ft-Base':fbaseline/step,'ft-Trans':ft,'ft-Trans-Filtered':None,'Trans-Filtered':None,'ft-Filter-Threshold':None}
      else :
         if fftFilt >= 0 :
            idx=np.where(abs(fbaseline)>abs(fftFilt))[0]
         else :
            idx=np.where(abs(fbaseline)<abs(fftFilt))[0]
         print "Applied FFT Filter ",fftFilt,' over ',len(idx),' elements '
         ff=ft*1
         if len(idx) > 0 : ff[idx]=0.
         TFilter=np.real(FFT.ifft(ff))+mT
         if forceToZero :
            idx=np.where(TFilter < 0)[0]
            if len(idx)>0 : 
               print "Forced to 0 negative values in ",len(idx)," samples "
               TFilter[idx]=0
         if renorm :
            print "Renormalization constant ",TFilter.sum()*step
            TFilter*=1/(TFilter.sum()*step)
         if asBP :
            out=self.copy()
            out.Trans=TFilter
            out.ft_Filter_Threshold=fftFilt
            return out
         else :
            return {'ft-zero':mT,'ft-step':step,'ft-Base':fbaseline/step,'ft-Trans':ft,'ft-Trans-Filtered':ff,'Trans-Filtered':TFilter,'ft-Filter-Threshold':fftFilt}
            
class _BandPass_Fits_Base :
   """BaseClass to handle a fits file"""
   def __init__(self,FitsName,__path__='.') :
      import numpy as np
      import pyfits
      self.DefaultPath="/".join(__path__)
      self.FitsName=FitsName+''
      self.FGHz_2_WN_cm=1./29.9792458
      try :
         self.f=pyfits.open(FitsName)
      except :
         self.f = None
         return
      self.fmap = {}
      count=0
      for k in self.f : 
         try :
            n=k.header['extname']
         except :
            n='PRIMARY'
         self.fmap[n] = count
         count+=1
      self.FMapParameters = self.Get_Frequency_Map_Parameters()
      self.FChanParameters = self.Get_Channel_Parameters()
   def __del__(self) :
      self.close()
   def find(self,Detector) :
      pass
   ##def __len__(self) :
      ##return len(self.
   #def tocsv(self,filename) :
      #try :
         #f=open(filename,'w')
      #except :
         #print "Can not open %s"%filename
         #return
      #f.write('FreqGHz,Trans\n')
      #for l in range(
      
      #f.close()
   def Get_Frequency_Map_Parameters(self,HDUname='FREQUENCY_MAP_PARAMETERS',HDUnumber=2) :
      HDUid=HDUname if HDUnumber == None else HDUnumber 
      try :
         L={}
         L['Header'] = str(self.f[HDUid].header).split('\n')
         for it in range(1,self.f[HDUid].header['TFIELDS']+1) :
            n=self.f[HDUid].header['TTYPE%d'%it]
            L[n.capitalize()]=self.f[HDUid].data.field(n)
         return L
      except :
         return None
   def find_on_frequency_map(self,frequency) :
      from numpy import array,where
      if self.FMapParameters == None :
         return array([])
      return where(self.FMapParameters['Frequency'] == ('%03d'%int(frequency)))[0]
   def Get_Channel_Parameters(self,HDUname='CHANNEL_PARAMETERS',HDUnumber=1) :
      HDUid=HDUname if HDUnumber == None else HDUnumber 
      try :
         L={}
         L['Header'] = str(self.f[HDUid].header).split('\n')
         for it in range(1,self.f[HDUid].header['TFIELDS']+1) :
            n=self.f[HDUid].header['TTYPE%d'%it]
            L[n.capitalize()]=self.f[HDUid].data.field(n)
         return L
      except :
         return None
   def find_on_channel_parameters(self,detector) :
      from numpy import array,where
      if self.FChanParameters == None :
         return array([])
      return where(self.FChanParameters['Detector'] == detector)[0]
   def Fields(self,ihdu) :
      """Return a dictionary which allows to identify the columns"""
      import numpy
      try :
         nf=self.f[ihdu].header['TFIELDS']
      except :
         return None
      if nf <= 0 :
         return None
      a={}
      for ifield in range(1,nf+1) :
         _N=self.f[ihdu].header['TTYPE%d'%ifield]
         _U=self.f[ihdu].header['TUNIT%d'%ifield]
         _F=self.f[ihdu].header['TFORM%d'%ifield]
         a[_N]={'column':ifield,'unit':_U,'form':_F.strip()}
      return a
   def close(self) :
      try :
         self.f.close()
      except :
         pass
   def __call__(self,Detector,skip_NullWavenumbers=True,skip_NullTransmission=False,skip_BadFlag=True,fft=False,fftFilt=None,useFmap=True,detectorIndex=None) :
      """The "wavenumber" column is converted to a Frequency"""
      import numpy as np
      if detectorIndex == None :
         if useFmap : 
            n=self.fmap[Detector]
         else :
            n=self.find(Detector)
         if Detector[0]=='F' :
            idx = self.find_on_frequency_map(Detector[1:])
            try :
               RefFreqGHz =self.FMapParameters['Centralfreq'][idx][0]
            except :
               "this is just to fix a bug in the HFI RIMO"
               RefFreqGHz =float(self.FMapParameters['Frequency'][idx][0])
         else :
            RefFreqGHz=None
      else :
         RefFreqGHz=None
         n=int(detectorIndex)
      header = str(self.f[n].header)
      extname=self.f[n].header['EXTNAME']
      x = self.f[n].data.field('WAVENUMBER')
      y = self.f[n].data.field('TRANSMISSION')
      try :
         err = self.f[n].data.field('UNCERTAINITY')
      except :
         print "column UNCERTAINTY not found, try with UNCERTAINTY"
         try :
            err = self.f[n].data.field('UNCERTAINITY')
         except :
            print "column UNCERTAINITY not found, error filled with NAN"
            err = x+np.nan
      try :
         flag = self.f[n].data.field('FLAG').strip()=='T'
      except :
         flag = self.f[n].data.field('FLAG')
      ok=np.ones(len(x))
      if skip_NullWavenumbers : ok*=x>0
      if skip_NullTransmission : ok*=y>0
      #if skip_BadFlag : ok*=flag==False
      idx = np.where(ok)[0]
      fld=self.Fields(n)
      history=''
      if fld['WAVENUMBER']['unit'].strip().lower()=='ghz' :
         w=x*self.FGHz_2_WN_cm
         f=x
         history='WAVENUMBER converted from GHz to 1/cm'
      elif fld['WAVENUMBER']['unit'].strip().lower()=='1/cm' or fld['WAVENUMBER']['unit'].strip().lower()=='cm^-1':
         f=x/self.FGHz_2_WN_cm
         w=x
         history='WAVENUMBER converted from 1/cm to GHz'
      else :
         f=x
         w=x
         history='WAVENUMBER taken as it is'
         print "Warning: it is not possible to interprete TUNIT for wavenumber"
         print fld['WAVENUMBER']['unit'].strip().lower()
      cml=y[idx].cumsum()
      out = BandPass()
      out['EXTNAME']=extname
      out['header']=header
      out['FreqGHz']=f[idx]
      out['Wavenumber']=w[idx]
      out['Trans']=y[idx]
      out['Error']=err[idx]
      out['Flag']=flag[idx]
      out['Chain']=Detector if Detector!= None else extname
      out['history']=history
      out['RefFreqGHz']=RefFreqGHz 
      out['Cumulant']=cml/cml.max()
      try :
         out['step']=(out.FreqGHz.max()-out.FreqGHz.min())/float(len(out)-1)
      except :
         out['step']=0
      #out={'EXTNAME':extname,'header':header,'FreqGHz':f[idx],'Wavenumber':w[idx],'Trans':y[idx],'Error':err[idx],'Flag':flag[idx],'Chain':Detector,'history':history,'RefFreqGHz':RefFreqGHz,'Cumulant':cml/cml.max()}
      if fft==True or fftFilt!=None:
         a=out.fft(fftFilt)
         for k in a.keys() : out[k]=a[k]
      return out
   #
   def StatsTable(self,fhlist=[18,19,20,21,22,23,24,25,26,27,28],armlist=['S','M']) :
      import numpy as np
      from SmartTable import base_table
      #
      FMapParameters=self.FMapParameters.copy()
      #
      FChanParameters=self.FChanParameters.copy()
      dtc=[]
      for k in FChanParameters['Detector'] : dtc.append(k[0].strip())
      FChanParameters['Detector']=np.array(dtc)
      #
      out=base_table()
      for k in "fh,arm,alt_arm,polarization,radiometer,Fcent,MAX,Bcmb,Brj1,dBdT_cmb".split(',') :
         out.newcolumn(k,[])
      #
      for k in "rimo_detector,rimo_fwhm_deg,rimo_fsamp".split(',') :
         out.newcolumn(k,[])
      #
      for k in "map_fwhm_deg,map_noise,map_central_freq".split(',') :
         out.newcolumn(k,[])
      #
      for fh in fhlist :
         for arm in armlist :
            alt_arm = Detector2Arm(arm)
            BP=self(str(fh)+arm)
            ST=BP.Statistics()
            out['fh'].append(fh)
            out['arm'].append(alt_arm[0])
            out['alt_arm'].append(str(alt_arm[0])+str(alt_arm[2])+str(alt_arm[1]))
            out['polarization'].append(alt_arm[1])
            out['radiometer'].append(alt_arm[2])
            out['Fcent'].append(ST['Fcent'])
            out['MAX'].append(ST['MAX'])
            out['Bcmb'].append(ST['Bcmb'])
            out['Brj1'].append(ST['Brj1'])
            out['dBdT_cmb'].append(ST['dBdT_cmb'])
            #
            kk='LFI'+str(fh)+arm
            i=np.where(FChanParameters['Detector']==kk)[0][0]
            out['rimo_detector'].append(FChanParameters['Detector'][i][0])
            out['rimo_fwhm_deg'].append(FChanParameters['Fwhm'][i]/60.)
            out['rimo_fsamp'].append(FChanParameters['F_samp'][i])
            #
            if fh < 24 :
               out['map_fwhm_deg'].append(FMapParameters['Fwhm'][2]/60.)
               out['map_noise'].append(FMapParameters['Noise'][2])
               out['map_central_freq'].append(FMapParameters['Centralfreq'][2])
            elif fh < 27 :
               out['map_fwhm_deg'].append(FMapParameters['Fwhm'][1]/60.)
               out['map_noise'].append(FMapParameters['Noise'][1])
               out['map_central_freq'].append(FMapParameters['Centralfreq'][1])
            else :
               out['map_fwhm_deg'].append(FMapParameters['Fwhm'][0]/60.)
               out['map_noise'].append(FMapParameters['Noise'][0])
               out['map_central_freq'].append(FMapParameters['Centralfreq'][0])
      #
      for k in out.keys() : out[k]=np.array(out[k])
      return out
      
class LFI_BandPass_Fits(_BandPass_Fits_Base) :
   """Class to handle an LFI bandpass
   
      Usage:
         LFI_BP=LFI_BandPass_Fits()
      uses the default RIMO included in the library distribution.
      
      The current Default RIMO can be discovered with
         print LFI_BandPass_Fits('').defaultRIMO()
      
      Alternativelly
      
         LFI_BP=LFI_BandPass_Fits('LFI_RIMO_PINCO.fits')
         
      uses the LFI RIMO contained into the file 'LFI_RIMO_PINCO.fits'.
   """
   def __init__(self,*arg) :
      if len(arg) == 0 :
         _BandPass_Fits_Base.__init__(self,self.defaultRIMO())
      else :
         _BandPass_Fits_Base.__init__(self,arg[0])
   def defaultRIMO(self) :
      """Returns the name and path of the default RIMO"""
      return _BandPass_Fits_Base('').DefaultPath+'/'+'LFI_RIMO_18092012_DX9d.fits'
   def find(self,Detector) :
      if Detector[0]=='F' :
         a='BANDPASS_'+Detector
         try :
            return self.fmap[a]
         except :
            return None
      a='BANDPASS_'
      try :
         n=int(Detector[0:2])
      except :
         return None
      if 27 <= n and n <= 28 :
         a+='030'
      elif 24 <= n and n <= 26 :
         a+='044'
      elif 18 <= n and n <= 23 :
         a+='070'
      else :
         return None
      a+='-'+Detector
      try :
         return self.fmap[a]
      except :
         return None
   def detector2horn(self,Detector) :
      n=int(Detector[0:2])
      if 27 <= n and n <= 28 :
         return 30
      elif 24 <= n and n <= 26 :
         return 44
      elif 18 <= n and n <= 23 :
         return 70
      else :
         return None
   def __call__(self,Detector,skip_NullWavenumbers=True,skip_NullTransmission=False,skip_BadFlag=True,fft=False,fftFilt=None,detectorIndex=None,justBP=False) :
      idx=self.find(Detector)
      if idx == None : return
      BP = _BandPass_Fits_Base.__call__(self,None,skip_NullWavenumbers=skip_NullWavenumbers,skip_NullTransmission=skip_NullTransmission,skip_BadFlag=skip_BadFlag,fft=fft,fftFilt=fftFilt,detectorIndex=idx)
      if justBP : return BP
      n=int(Detector[0:2])
      horn = self.detector2horn(Detector)
      idx = self.find_on_frequency_map('0'+str(horn))
      BP.RefFreqGHz =self.FMapParameters['Centralfreq'][idx][0]
      return BP

class HFI_BandPass_Fits(_BandPass_Fits_Base) :
   """Class to handle an HFI bandpass
   
      Usage:
         HFI_BP=HFI_BandPass_Fits()
      uses the default RIMO included in the library distribution.
      
      The current Default RIMO can be discovered with
         print HFI_BandPass_Fits('').defaultRIMO()
      
      Alternativelly
      
         HFI_BP=HFI_BandPass_Fits('HFI_RIMO_PALLO.fits')
         
      uses the HFI RIMO contained into the file 'HFI_RIMO_PALLO.fits'.
   """
   def __init__(self,*arg) :
      if len(arg) == 0 :
         _BandPass_Fits_Base.__init__(self,self.defaultRIMO())
      else :
         _BandPass_Fits_Base.__init__(self,arg[0])
   def defaultRIMO(self) :
      """Returns the name and path of the default RIMO"""
      return _BandPass_Fits_Base('').DefaultPath+'/'+'HFI-RIMO-20111017_2_71.fits'
   def find(self,Detector) :
      a='BANDPASS_'+Detector
      try :
         return self.fmap[a]
      except :
         return None
         

if __name__=='__main__' :
   import sys
   if len(sys.argv) < 2 :
      print """
Tabulates a list of statistics for the RIMO

Usage :
   > python rimo_fits.py <rimo_fits_file.fits>

Units are K, GHz, MJy/Sr
"""
      sys.exit(0)
   print "Reading ",sys.argv[1]
   RIMO=rimo_fits.LFI_BandPass_Fits(sys.argv[1])
   if RIMO.f == None :
      print "Error, file not found"
      sys.exit(0)

   print "fh,arm,polarization,radiometer,Fcent,MAX,Bcmb,Brj1,dBdT_cmb"
   for fh in range(18,29) :
      for arm in ['S','M'] :
         alt_arm = rimo_fits.Detector2Arm(arm)
         BP=RIMO(str(fh)+arm)
         ST=BP.Statistics()
         line=str(fh)
         line+=','+str(alt_arm[0])
         line+=','+str(alt_arm[1])
         line+=','+str(alt_arm[2])
         line+=',%16.12e'%ST['Fcent']
         line+=',%16.12e'%ST['MAX']
         line+=',%16.12e'%ST['Bcmb']
         line+=',%16.12e'%ST['Brj1']
         line+=',%16.12e'%ST['dBdT_cmb']
         print line
