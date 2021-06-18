# -*- coding: utf-8 -*-
"""
"""

#__________________________________________________________________________________
#
# data paths and probe settings 


#meanXYFlag = 1 # 1= mean xy, 0= max
#MaxAspectRatio = 10 # Airspeed diameter / array diameter
#Streak = 5 # Exclude particles 5 long and 1 wide
#BiggestParticle  # 0 = BBox, 1 = largest particle
#ThresholdSize = Size to switch between stereo and standard psd
#ThresholdDeltaDiameterY = um allowed difference in y diameter for stereo, -1 no threshold


import numpy as np


def GetFlightInfo2DS():

    Info2DS = {}

    # Probe parameters
    Info2DS['Lambda_um'] = 0.785	# in um 
    Info2DS['PixelSize']= 10
    Info2DS['ArrayElements'] = 128
    Info2DS['c']=8
    
    # Acceptance parameters
    #Info2DS['IATmin'] =1E-5
    Info2DS['Streak'] = 5
    Info2DS['MaxAspectRatio'] =10 
    
    #Flight and data directory    
    Info2DS['C078','Path2DS']= 'D:/PICASSO/rawdata/C078/2DS/Output/'
    Info2DS['C078','Path2DSsave']= 'D:/PICASSO/rawdata/C078/2DS/Output/Diagnostics/'   
    Info2DS['C078', 'IAT_threshold'] =1E-5
    Info2DS['C078','TAS']=100 #m/s
    Info2DS['C078','ArmSep'] = 63

    Info2DS['B895','Path2DS']= 'C:/Users/mbexjso2/OneDrive - The University of Manchester/Documents/Work/B895/2DS/'
    Info2DS['B895','Path2DSsave']= r"C:/Users/mbexjso2/OneDrive - The University of Manchester/Documents/Work/B895/2DS/Colocation/"
    Info2DS['B895','FlightDate'] = np.datetime64('2015-03-13 00:00:00')
    Info2DS['B895', 'IAT_threshold'] =1E-5
    Info2DS['B895', 'ColocationThreshold'] = 5E-7
    Info2DS['B895','TAS']=100 #m/s
    Info2DS['B895','ThresholdDeltaDiameterY']= -1 # um allowed difference in y diameter for stereo, -1 not threshold
    Info2DS['B895','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['B895','FlightNumber'] = 'b895'
    Info2DS['B895','MeanXYFlag'] = 0 # 1= mean xy, 0= max
    Info2DS['B895','ArmSep'] = 63
    
    Info2DS['C031','Path2DS']= 'C:/Users/mbexjso2/OneDrive - The University of Manchester/Documents/Work/Clarify/C031/'
    Info2DS['C031','Path2DSsave']= 'C:/Users/mbexjso2/OneDrive - The University of Manchester/Documents/Work/Clarify/C031/Colocation/'
    Info2DS['C031','FlightDate'] = np.datetime64('2015-08-17 00:00:00')
    Info2DS['C031', 'IAT_threshold'] =1E-6
    Info2DS['C031', 'ColocationThreshold'] = 5E-7
    Info2DS['C031','TAS']=100 #m/s
    Info2DS['C031','ThresholdDeltaDiameterY']= -1 # um allowed difference in y diameter for stereo, -1 not threshold
    Info2DS['C031','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['C031','FlightNumber'] = 'c031'
    Info2DS['C031','ArmSep'] = 63
    
    Info2DS['B895mphys','Path2DS']= 'C:/Users/mbexjso2/OneDrive - The University of Manchester/Documents/Work/Mphys/Cirrcrex/B895/'
    Info2DS['B895mphys','Path2DSsave']= 'C:/Users/mbexjso2/OneDrive - The University of Manchester/Documents/Work/Mphys/Cirrcrex/B895/Colocation/'
    Info2DS['B895mphys','FlightDate'] = np.datetime64('2015-03-13 00:00:00')
    Info2DS['B895mphys', 'IAT_threshold'] =1E-5
    Info2DS['B895mphys', 'ColocationThreshold'] = 5E-7
    Info2DS['B895mphys','TAS']=100 #m/s
    Info2DS['B895mphys','ArmSep'] = 63
    
    Info2DS['B895_dataPC','Path2DS']= 'D:/CIRCCREX/B895/'
    Info2DS['B895_dataPC','Path2DSsave']= 'D:/CIRCCREX/B895/Colocation/'
    Info2DS['B895_dataPC','FlightDate'] = np.datetime64('2015-03-13 00:00:00')
    Info2DS['B895_dataPC', 'IAT_threshold'] =1E-5
    Info2DS['B895_dataPC', 'ColocationThreshold'] = 5E-7
    Info2DS['B895_dataPC','TAS']=100 #m/s
    Info2DS['B895_dataPC','ThresholdDeltaDiameterY']= -1 # um allowed difference in y diameter for stereo, -1 not threshold
    Info2DS['B895_dataPC','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['B895_dataPC','FlightNumber'] = 'b895'
    Info2DS['B895_dataPC','ArmSep'] = 63
    
    Info2DS['B894_dataPC','Path2DS']= 'D:/CIRCCREX/B894/2ds/OasisOutput/'
    Info2DS['B894_dataPC','Path2DSsave']= 'D:/CIRCCREX/B894/2ds/OasisOutput/Colocation/'
    Info2DS['B894_dataPC','FlightDate'] = np.datetime64('2015-03-11 00:00:00')
    Info2DS['B894_dataPC', 'IAT_threshold'] =1E-5
    Info2DS['B894_dataPC', 'ColocationThreshold'] = 5E-7
    Info2DS['B894_dataPC','TAS']=100 #m/s
    Info2DS['B894_dataPC','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 not threshold
    Info2DS['B894_dataPC','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['B894_dataPC','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['B894_dataPC','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['B894_dataPC','FlightNumber'] = 'b894'
    Info2DS['B894_dataPC','ArmSep'] = 63
    
    
    Info2DS['C174_dataPC','Path2DS']= 'D:/PICASSO/rawdata/C174/2ds/Output/'
    Info2DS['C174_dataPC','Path2DSsave']= 'D:/PICASSO/rawdata/C174/2ds/Output/Colocation/'
    Info2DS['C174_dataPC','FlightDate'] = np.datetime64('2019-05-23 00:00:00')
    Info2DS['C174_dataPC', 'IAT_threshold'] =1E-5
    Info2DS['C174_dataPC', 'ColocationThreshold'] = 5E-6
    Info2DS['C174_dataPC','TAS']=100 #m/s
    Info2DS['C174_dataPC','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['C174_dataPC','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['C174_dataPC','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['C174_dataPC','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['C174_dataPC','FlightNumber'] = 'c174'
    Info2DS['C174_dataPC','ArmSep'] = 63
    Info2DS['C174_dataPC','CorePath'] ='D:/PICASSO/rawdata/FAAM_Data/AllCoreFaam/'
    Info2DS['C174_dataPC','CoreFileName'] ='core_faam_20190523_v004_r0_c174_1hz.nc'

    Info2DS['C174_3V','Path2DS']= 'D:/PICASSO/rawdata/C174/3vcpi/2DS/Output/'
    Info2DS['C174_3V','Path2DSsave']= 'D:/PICASSO/rawdata/C174/3vcpi/2DS/Output/Colocation/'
    Info2DS['C174_3V','FlightDate'] = np.datetime64('2019-05-23 00:00:00')
    Info2DS['C174_3V', 'IAT_threshold'] =1E-5
    Info2DS['C174_3V', 'ColocationThreshold'] = 1E-6
    Info2DS['C174_3V','TAS']=100 #m/s
    Info2DS['C174_3V','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['C174_3V','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['C174_3V','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['C174_3V','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['C174_3V','ArmSep'] = 50
    Info2DS['C174_3V','FlightNumber'] = 'c174'
    Info2DS['C174_3V','CorePath'] ='D:/PICASSO/rawdata/FAAM_Data/AllCoreFaam/'
    Info2DS['C174_3V','CoreFileName'] ='core_faam_20190523_v004_r0_c174_1hz.nc'

    Info2DS['C172_dataPC','Path2DS']= 'D:/PICASSO/rawdata/C172/2ds/Output/'
    Info2DS['C172_dataPC','Path2DSsave']= 'D:/PICASSO/rawdata/C172/2ds/Output/Colocation/'
    Info2DS['C172_dataPC','FlightDate'] = np.datetime64('2019-05-07 00:00:00')
    Info2DS['C172_dataPC', 'IAT_threshold'] =1E-5
    Info2DS['C172_dataPC', 'ColocationThreshold'] = 5E-6
    Info2DS['C172_dataPC','TAS']=100 #m/s
    Info2DS['C172_dataPC','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['C172_dataPC','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['C172_dataPC','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['C172_dataPC','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['C172_dataPC','FlightNumber'] = 'c172'
    Info2DS['C172_dataPC','ArmSep'] = 63
    Info2DS['C172_dataPC','CorePath'] ='D:/PICASSO/rawdata/FAAM_Data/AllCoreFaam/'
    Info2DS['C172_dataPC','CoreFileName'] ='core_faam_20190507_v004_r0_c172_1hz.nc'

    Info2DS['C171_dataPC','Path2DS']= 'D:/PICASSO/rawdata/C171/2ds/Output/'
    Info2DS['C171_dataPC','Path2DSsave']= 'D:/PICASSO/rawdata/C171/2ds/Output/Colocation/'
    Info2DS['C171_dataPC','FlightDate'] = np.datetime64('2019-05-01 00:00:00')
    Info2DS['C171_dataPC', 'IAT_threshold'] =1E-5
    Info2DS['C171_dataPC', 'ColocationThreshold'] = 5E-6
    Info2DS['C171_dataPC','TAS']=100 #m/s
    Info2DS['C171_dataPC','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['C171_dataPC','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['C171_dataPC','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['C171_dataPC','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['C171_dataPC','FlightNumber'] = 'c171'
    Info2DS['C171_dataPC','ArmSep'] = 63
    Info2DS['C171_dataPC','CorePath'] ='D:/PICASSO/rawdata/FAAM_Data/AllCoreFaam/'
    Info2DS['C171_dataPC','CoreFileName'] ='core_faam_20190501_v004_r0_c171_1hz.nc'
    
    Info2DS['C171_3V','Path2DS']= 'D:/PICASSO/rawdata/C171/3V-CPI/Output/'
    Info2DS['C171_3V','Path2DSsave']= 'D:/PICASSO/rawdata/C171/3V-CPI/Output/Colocation/'
    Info2DS['C171_3V','FlightDate'] = np.datetime64('2019-05-01 00:00:00')
    Info2DS['C171_3V', 'IAT_threshold'] =1E-5
    Info2DS['C171_3V', 'ColocationThreshold'] = 1E-6
    Info2DS['C171_3V','TAS']=100 #m/s
    Info2DS['C171_3V','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 not threshold
    Info2DS['C171_3V','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['C171_3V','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['C171_3V','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['C171_3V','ArmSep'] = 50
    Info2DS['C171_3V','FlightNumber'] = 'c171'
    Info2DS['C171_3V','CorePath'] ='D:/PICASSO/rawdata/FAAM_Data/AllCoreFaam/'
    Info2DS['C171_3V','CoreFileName'] ='core_faam_20190501_v004_r0_c171_1hz.nc'

    Info2DS['C170_dataPC','Path2DS']= 'D:/PICASSO/rawdata/C170/2ds/Output/'
    Info2DS['C170_dataPC','Path2DSsave']= 'D:/PICASSO/rawdata/C170/2ds/Output/Colocation/'
    Info2DS['C170_dataPC','FlightDate'] = np.datetime64('2019-04-26 00:00:00')
    Info2DS['C170_dataPC', 'IAT_threshold'] =1E-5
    Info2DS['C170_dataPC', 'ColocationThreshold'] = 5E-6
    Info2DS['C170_dataPC','TAS']=100 #m/s
    Info2DS['C170_dataPC','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['C170_dataPC','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['C170_dataPC','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['C170_dataPC','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['C170_dataPC','FlightNumber'] = 'c170'
    Info2DS['C170_dataPC','ArmSep'] = 63
    Info2DS['C170_dataPC','CorePath'] ='D:/PICASSO/rawdata/FAAM_Data/AllCoreFaam/'
    Info2DS['C170_dataPC','CoreFileName'] ='core_faam_20190426_v004_r0_c170_1hz.nc'


    Info2DS['C169_dataPC','Path2DS']= 'D:/PICASSO/rawdata/C169/2ds/Output/'
    Info2DS['C169_dataPC','Path2DSsave']= 'D:/PICASSO/rawdata/C169/2ds/Output/Colocation/'
    Info2DS['C169_dataPC','FlightDate'] = np.datetime64('2019-04-16 00:00:00')
    Info2DS['C169_dataPC', 'IAT_threshold'] =1E-5
    Info2DS['C169_dataPC', 'ColocationThreshold'] = 5E-6
    Info2DS['C169_dataPC','TAS']=100 #m/s
    Info2DS['C169_dataPC','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['C169_dataPC','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['C169_dataPC','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['C169_dataPC','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['C169_dataPC','FlightNumber'] = 'c169'
    Info2DS['C169_dataPC','ArmSep'] = 63
    Info2DS['C169_dataPC','CorePath'] ='D:/PICASSO/rawdata/FAAM_Data/AllCoreFaam/'
    Info2DS['C169_dataPC','CoreFileName'] ='core_faam_20190416_v004_r0_c169_1hz.nc'

    Info2DS['C098_dataPC','Path2DS']= 'D:/PICASSO/rawdata/C098/2ds/Output/'
    Info2DS['C098_dataPC','Path2DSsave']= 'D:/PICASSO/rawdata/C098/2ds/Output/Colocation/'
    Info2DS['C098_dataPC','FlightDate'] = np.datetime64('2018-04-24 00:00:00')
    Info2DS['C098_dataPC', 'ColocationThreshold'] = 1E-6
    Info2DS['C098_dataPC', 'IAT_threshold'] =1E-5
    Info2DS['C098_dataPC','TAS']=100 #m/s
    Info2DS['C098_dataPC','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['C098_dataPC','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['C098_dataPC','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['C098_dataPC','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['C098_dataPC','FlightNumber'] = 'c098'
    Info2DS['C098_dataPC','ArmSep'] = 63
    Info2DS['C098_dataPC','CorePath'] ='D:/PICASSO/rawdata/FAAM_Data/AllCoreFaam/'
    Info2DS['C098_dataPC','CoreFileName'] ='core_faam_20180424_v004_r0_c098_1hz.nc'

    Info2DS['C097_dataPC','Path2DS']= 'D:/PICASSO/rawdata/C097/2ds/Output/'
    Info2DS['C097_dataPC','Path2DSsave']= 'D:/PICASSO/rawdata/C097/2ds/Output/Colocation/'
    Info2DS['C097_dataPC','FlightDate'] = np.datetime64('2018-04-23 00:00:00')
    Info2DS['C097_dataPC', 'ColocationThreshold'] = 1E-6
    Info2DS['C097_dataPC', 'IAT_threshold'] =1E-5
    Info2DS['C097_dataPC','TAS']=100 #m/s
    Info2DS['C097_dataPC','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['C097_dataPC','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['C097_dataPC','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['C097_dataPC','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['C097_dataPC','FlightNumber'] = 'c097'
    Info2DS['C097_dataPC','ArmSep'] = 63
    Info2DS['C097_dataPC','CorePath'] ='D:/PICASSO/rawdata/FAAM_Data/AllCoreFaam/'
    Info2DS['C097_dataPC','CoreFileName'] ='core_faam_20180423_v004_r0_c097_1hz.nc'


    Info2DS['TwinOtter354','Path2DS']= 'D:/EUREC4A/TwinOtter/354/2DS/Output/'
    Info2DS['TwinOtter354','Path2DSsave']= 'D:/EUREC4A/TwinOtter/354/2DS/Output/Colocation/'
    Info2DS['TwinOtter354','FlightDate'] = np.datetime64('2020-02-15 00:00:00')
    Info2DS['TwinOtter354', 'IAT_threshold'] =1E-5
    Info2DS['TwinOtter354', 'ColocationThreshold'] = 5E-6
    Info2DS['TwinOtter354','TAS']=100 #m/s

    Info2DS['TwinOtter353','Path2DS']= 'D:/EUREC4A/TwinOtter/353/2DS/Output/'
    Info2DS['TwinOtter353','Path2DSsave']= 'D:/EUREC4A/TwinOtter/353/2DS/Output/Colocation/'
    Info2DS['TwinOtter353','FlightDate'] = np.datetime64('2020-02-15 00:00:00')
    Info2DS['TwinOtter353', 'IAT_threshold'] =1E-5
    Info2DS['TwinOtter353', 'ColocationThreshold'] = 5E-6
    Info2DS['TwinOtter353','TAS']=100 #m/s

    Info2DS['B926','Path2DS']= 'D:/ICE_D/B926/2DS/Output/'
    Info2DS['B926','Path2DSsave']= 'D:/ICE_D/B926/2DS/Output/Colocation/'
    Info2DS['B926','FlightDate'] = np.datetime64('2015-08-14 00:00:00')
    Info2DS['B926', 'IAT_threshold'] =1E-5
    Info2DS['B926', 'ColocationThreshold'] = 5E-7
    Info2DS['B926','TAS']=100 #m/s

    Info2DS['C050','Path2DS']= 'D:/CLARIFY/C050/2DS/Output/'
    Info2DS['C050','Path2DSsave']= 'D:/CLARIFY/C050/2DS/Output/Colocation/'
    Info2DS['C050','FlightDate'] = np.datetime64('2017-09-04 00:00:00')
    Info2DS['C050', 'IAT_threshold'] =1E-5
    Info2DS['C050', 'ColocationThreshold'] = 5E-7
    Info2DS['C050','TAS']=100 #m/s
    Info2DS['C050','ThresholdDeltaDiameterY']= -1 # um allowed difference in y diameter for stereo, -1 not threshold
    Info2DS['C050','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['C050','FlightNumber'] = 'c050'

    Info2DS['CLACE','Path2DS']= 'D:/CLACE2017/3V2DS/Output/'
    Info2DS['CLACE','Path2DSsave']= 'D:/CLACE2017/3V2DS/Output/Colocation/'
    Info2DS['CLACE','FlightDate'] = np.datetime64('2019-02-12 00:00:00')
    Info2DS['CLACE', 'IAT_threshold'] =3E-5
    Info2DS['CLACE', 'ColocationThreshold'] = 8E-6
    Info2DS['CLACE','TAS']=20 #m/s
    Info2DS['CLACE','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['CLACE','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['CLACE','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['CLACE','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['CLACE','ArmSep'] = 63
    
    Info2DS['C164','Path2DS']= 'D:/PIKNMIX_2019/C164/2ds/Output/'
    Info2DS['C164','Path2DSsave']= 'D:/PIKNMIX_2019/C164/2ds/Output/Colocation/'
    Info2DS['C164','FlightDate'] = np.datetime64('2019-03-25 00:00:00')
    Info2DS['C164', 'ColocationThreshold'] = 1E-6
    Info2DS['C164', 'IAT_threshold'] =1E-5
    Info2DS['C164','TAS']=100 #m/s
    Info2DS['C164','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['C164','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['C164','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['C164','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['C164','FlightNumber'] = 'C164'
    Info2DS['C164','ArmSep'] = 63
    Info2DS['C164','CorePath'] ='D:/PIKNMIX_2019/FAAM_data/c164-mar-25/core_processed/'
    Info2DS['C164','CoreFileName'] ='core_faam_20190325_v004_r0_c164_1hz.nc'
    
    Info2DS['C168','Path2DS']= 'D:/PIKNMIX_2019/C168/2DS/Output/'
    Info2DS['C168','Path2DSsave']= 'D:/PIKNMIX_2019/C168/2DS/Output/Colocation/'
    Info2DS['C168','FlightDate'] = np.datetime64('2019-03-28 00:00:00')
    Info2DS['C168', 'ColocationThreshold'] = 1E-6
    Info2DS['C168', 'IAT_threshold'] =1E-6
    Info2DS['C168','TAS']=100 #m/s
    Info2DS['C168','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['C168','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['C168','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['C168','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['C168','FlightNumber'] = 'C168'
    Info2DS['C168','ArmSep'] = 63
    Info2DS['C168','CorePath'] ='D:/PIKNMIX_2019/FAAM_data/c168-mar-28/core_processed/'
    Info2DS['C168','CoreFileName'] ='core_faam_20190328_v004_r0_c168_1hz.nc'
    
    Info2DS['MAC_218','Path2DS']= 'D:/MAC/218/2DS/Output/'
    Info2DS['MAC_218','Path2DSsave']= 'D:/MAC/218/2DS/Output/Colocation/'
    Info2DS['MAC_218','FlightDate'] = np.datetime64('2015-11-27 00:00:00')
    Info2DS['MAC_218', 'ColocationThreshold'] = 1E-6
    Info2DS['MAC_218', 'IAT_threshold'] =1E-6
    Info2DS['MAC_218','TAS']=60 #m/s
    Info2DS['MAC_218','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['MAC_218','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['MAC_218','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['MAC_218','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['MAC_218','FlightNumber'] = '218'
    Info2DS['MAC_218','ArmSep'] = 63
    #Info2DS['MAC_218','CorePath'] ='D:/PICASSO/rawdata/FAAM_Data/AllCoreFaam/'
    #Info2DS['MAC_218','CoreFileName'] ='core_faam_20180423_v004_r0_c097_1hz.nc'
    
    Info2DS['MAC_219','Path2DS']= 'D:/MAC/219/2d-s/OasisOutput/'
    Info2DS['MAC_219','Path2DSsave']= 'D:/MAC/219/2d-s/OasisOutput/Colocation/'
    Info2DS['MAC_219','FlightDate'] = np.datetime64('2015-11-27 00:00:00')
    Info2DS['MAC_219', 'ColocationThreshold'] = 1E-6
    Info2DS['MAC_219', 'IAT_threshold'] =1E-6
    Info2DS['MAC_219','TAS']=60 #m/s
    Info2DS['MAC_219','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['MAC_219','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['MAC_219','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['MAC_219','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['MAC_219','FlightNumber'] = '219'
    Info2DS['MAC_219','ArmSep'] = 63

    Info2DS['MAC_220','Path2DS']= 'D:/MAC/220/2d-s/OasisOutput/'
    Info2DS['MAC_220','Path2DSsave']= 'D:/MAC/220/2d-s/OasisOutput/Colocation/'
    Info2DS['MAC_220','FlightDate'] = np.datetime64('2015-11-28 00:00:00')
    Info2DS['MAC_220', 'ColocationThreshold'] = 1E-6
    Info2DS['MAC_220', 'IAT_threshold'] =1E-6
    Info2DS['MAC_220','TAS']=60 #m/s
    Info2DS['MAC_220','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['MAC_220','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['MAC_220','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['MAC_220','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['MAC_220','FlightNumber'] = '220'
    Info2DS['MAC_220','ArmSep'] = 63
    
    Info2DS['MAC_221','Path2DS']= 'D:/MAC/221/2d-s/OasisOutput/'
    Info2DS['MAC_221','Path2DSsave']= 'D:/MAC/221/2d-s/OasisOutput/Colocation/'
    Info2DS['MAC_221','FlightDate'] = np.datetime64('2015-11-29 00:00:00')
    Info2DS['MAC_221', 'ColocationThreshold'] = 1E-6
    Info2DS['MAC_221', 'IAT_threshold'] =1E-6
    Info2DS['MAC_221','TAS']=60 #m/s
    Info2DS['MAC_221','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['MAC_221','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['MAC_221','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['MAC_221','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['MAC_221','FlightNumber'] = '221'
    Info2DS['MAC_221','ArmSep'] = 63
    
    Info2DS['MAC_222','Path2DS']= 'D:/MAC/222/2d-s/OasisOutput/'
    Info2DS['MAC_222','Path2DSsave']= 'D:/MAC/222/2d-s/OasisOutput/Colocation/'
    Info2DS['MAC_222','FlightDate'] = np.datetime64('2015-11-30 00:00:00')
    Info2DS['MAC_222', 'ColocationThreshold'] = 1E-6
    Info2DS['MAC_222', 'IAT_threshold'] =1E-6
    Info2DS['MAC_222','TAS']=60 #m/s
    Info2DS['MAC_222','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['MAC_222','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['MAC_222','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['MAC_222','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['MAC_222','FlightNumber'] = '222'
    Info2DS['MAC_222','ArmSep'] = 63
    
    Info2DS['MAC_223','Path2DS']= 'D:/MAC/223/2d-s/OasisOut/'
    Info2DS['MAC_223','Path2DSsave']= 'D:/MAC/223/2d-s/OasisOut/Colocation/'
    Info2DS['MAC_223','FlightDate'] = np.datetime64('2015-12-03 00:00:00')
    Info2DS['MAC_223', 'ColocationThreshold'] = 1E-6
    Info2DS['MAC_223', 'IAT_threshold'] =1E-6
    Info2DS['MAC_223','TAS']=60 #m/s
    Info2DS['MAC_223','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['MAC_223','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['MAC_223','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['MAC_223','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['MAC_223','FlightNumber'] = '223'
    Info2DS['MAC_223','ArmSep'] = 63
    
    Info2DS['MAC_224','Path2DS']= 'D:/MAC/224/2d-s/OasisOut/'
    Info2DS['MAC_224','Path2DSsave']= 'D:/MAC/224/2d-s/OasisOut/Colocation/'
    Info2DS['MAC_224','FlightDate'] = np.datetime64('2015-12-06 00:00:00')
    Info2DS['MAC_224', 'ColocationThreshold'] = 1E-6
    Info2DS['MAC_224', 'IAT_threshold'] =1E-6
    Info2DS['MAC_224','TAS']=60 #m/s
    Info2DS['MAC_224','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['MAC_224','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['MAC_224','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['MAC_224','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['MAC_224','FlightNumber'] = '224'
    Info2DS['MAC_224','ArmSep'] = 63
    
    Info2DS['MAC_225','Path2DS']= 'D:/MAC/225/2d-s/OasisOut/'
    Info2DS['MAC_225','Path2DSsave']= 'D:/MAC/225/2d-s/OasisOut/Colocation/'
    Info2DS['MAC_225','FlightDate'] = np.datetime64('2015-12-07 00:00:00')
    Info2DS['MAC_225', 'ColocationThreshold'] = 3E-6
    Info2DS['MAC_225', 'IAT_threshold'] =2E-6
    Info2DS['MAC_225','TAS']=60 #m/s
    Info2DS['MAC_225','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['MAC_225','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['MAC_225','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['MAC_225','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['MAC_225','FlightNumber'] = '225'
    Info2DS['MAC_225','ArmSep'] = 63
    
    Info2DS['MAC_226','Path2DS']= 'D:/MAC/226/2d-s/OasisOutput/'
    Info2DS['MAC_226','Path2DSsave']= 'D:/MAC/226/2d-s/OasisOutput/Colocation/'
    Info2DS['MAC_226','FlightDate'] = np.datetime64('2015-12-07 00:00:00')
    Info2DS['MAC_226', 'ColocationThreshold'] = 2E-6
    Info2DS['MAC_226', 'IAT_threshold'] =1E-6
    Info2DS['MAC_226','TAS']=60 #m/s
    Info2DS['MAC_226','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['MAC_226','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['MAC_226','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['MAC_226','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['MAC_226','FlightNumber'] = '226'
    Info2DS['MAC_226','ArmSep'] = 63
    
    Info2DS['MAC_227','Path2DS']= 'D:/MAC/227/2d-s/OasisOutput/'
    Info2DS['MAC_227','Path2DSsave']= 'D:/MAC/227/2d-s/OasisOutput/Colocation/'
    Info2DS['MAC_227','FlightDate'] = np.datetime64('2015-12-08 00:00:00')
    Info2DS['MAC_227', 'ColocationThreshold'] = 2E-6
    Info2DS['MAC_227', 'IAT_threshold'] =1E-6
    Info2DS['MAC_227','TAS']=60 #m/s
    Info2DS['MAC_227','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['MAC_227','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['MAC_227','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['MAC_227','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['MAC_227','FlightNumber'] = '227'
    Info2DS['MAC_227','ArmSep'] = 63
    
    Info2DS['MAC_228','Path2DS']= 'D:/MAC/228/2d-s/OasisOutput/'
    Info2DS['MAC_228','Path2DSsave']= 'D:/MAC/228/2d-s/OasisOutput/Colocation/'
    Info2DS['MAC_228','FlightDate'] = np.datetime64('2015-12-09 00:00:00')
    Info2DS['MAC_228', 'ColocationThreshold'] = 2E-6
    Info2DS['MAC_228', 'IAT_threshold'] =1E-6
    Info2DS['MAC_228','TAS']=60 #m/s
    Info2DS['MAC_228','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['MAC_228','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['MAC_228','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['MAC_228','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['MAC_228','FlightNumber'] = '228'
    Info2DS['MAC_228','ArmSep'] = 63
    
    Info2DS['MAC_230','Path2DS']= 'D:/MAC/230/2d-s/OasisOutput/'
    Info2DS['MAC_230','Path2DSsave']= 'D:/MAC/230/2d-s/OasisOutput/Colocation/'
    Info2DS['MAC_230','FlightDate'] = np.datetime64('2015-12-10 00:00:00')
    Info2DS['MAC_230', 'ColocationThreshold'] = 2E-6
    Info2DS['MAC_230', 'IAT_threshold'] =1E-6
    Info2DS['MAC_230','TAS']=60 #m/s
    Info2DS['MAC_230','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['MAC_230','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['MAC_230','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['MAC_230','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['MAC_230','FlightNumber'] = '230'
    Info2DS['MAC_230','ArmSep'] = 63
    
    Info2DS['MAC_231','Path2DS']= 'D:/MAC/231/2d-s/OasisOut/'
    Info2DS['MAC_231','Path2DSsave']= 'D:/MAC/231/2d-s/OasisOut/Colocation/'
    Info2DS['MAC_231','FlightDate'] = np.datetime64('2015-12-11 00:00:00')
    Info2DS['MAC_231', 'ColocationThreshold'] = 2E-6
    Info2DS['MAC_231', 'IAT_threshold'] =1E-6
    Info2DS['MAC_231','TAS']=60 #m/s
    Info2DS['MAC_231','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['MAC_231','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['MAC_231','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['MAC_231','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['MAC_231','FlightNumber'] = '231'
    Info2DS['MAC_231','ArmSep'] = 63
    
    Info2DS['MAC_232','Path2DS']= 'D:/MAC/232/2d-s/OasisOutput/'
    Info2DS['MAC_232','Path2DSsave']= 'D:/MAC/232/2d-s/OasisOutput/Colocation/'
    Info2DS['MAC_232','FlightDate'] = np.datetime64('2015-12-11 00:00:00')
    Info2DS['MAC_232', 'ColocationThreshold'] = 2E-6
    Info2DS['MAC_232', 'IAT_threshold'] =1E-6
    Info2DS['MAC_232','TAS']=60 #m/s
    Info2DS['MAC_232','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['MAC_232','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['MAC_232','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['MAC_232','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['MAC_232','FlightNumber'] = '232'
    Info2DS['MAC_232','ArmSep'] = 63
    
    Info2DS['MAC_233','Path2DS']= 'D:/MAC/233/2d-s/OasisOutput/'
    Info2DS['MAC_233','Path2DSsave']= 'D:/MAC/233/2d-s/OasisOutput/Colocation/'
    Info2DS['MAC_233','FlightDate'] = np.datetime64('2015-12-11 00:00:00')
    Info2DS['MAC_233', 'ColocationThreshold'] = 2E-6
    Info2DS['MAC_233', 'IAT_threshold'] =1E-6
    Info2DS['MAC_233','TAS']=60 #m/s
    Info2DS['MAC_233','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['MAC_233','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['MAC_233','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['MAC_233','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['MAC_233','FlightNumber'] = '233'
    Info2DS['MAC_233','ArmSep'] = 63
    
    Info2DS['MAC_234','Path2DS']= 'D:/MAC/234/2d-s/OasisOut/'
    Info2DS['MAC_234','Path2DSsave']= 'D:/MAC/234/2d-s/OasisOut/Colocation/'
    Info2DS['MAC_234','FlightDate'] = np.datetime64('2015-12-13 00:00:00')
    Info2DS['MAC_234', 'ColocationThreshold'] = 2E-6
    Info2DS['MAC_234', 'IAT_threshold'] =1E-6
    Info2DS['MAC_234','TAS']=60 #m/s
    Info2DS['MAC_234','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['MAC_234','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['MAC_234','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['MAC_234','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['MAC_234','FlightNumber'] = '234'
    Info2DS['MAC_234','ArmSep'] = 63
    
    Info2DS['MAC_235','Path2DS']= 'D:/MAC/235/2d-s/OasisOut/'
    Info2DS['MAC_235','Path2DSsave']= 'D:/MAC/235/2d-s/OasisOut/Colocation/'
    Info2DS['MAC_235','FlightDate'] = np.datetime64('2015-12-14 00:00:00')
    Info2DS['MAC_235', 'ColocationThreshold'] = 2E-6
    Info2DS['MAC_235', 'IAT_threshold'] =1E-6
    Info2DS['MAC_235','TAS']=60 #m/s
    Info2DS['MAC_235','ThresholdDeltaDiameterY']= 40 # um allowed difference in y diameter for stereo, -1 no threshold
    Info2DS['MAC_235','ThresholdSize'] = 300 # Size to switch between stereo and standard psd
    Info2DS['MAC_235','MeanXYFlag'] = 1 # 1= mean xy, 0= max
    Info2DS['MAC_235','BiggestParticle'] = 1 # #BiggestParticle  # 0 = BBox, 1 = largest particle
    Info2DS['MAC_235','FlightNumber'] = '235'
    Info2DS['MAC_235','ArmSep'] = 63
    
    
    return Info2DS



#__________________________________________________________________________________
#