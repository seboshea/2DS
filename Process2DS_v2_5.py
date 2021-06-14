# -*- coding: utf-8 -*-

#Procedures to identify stereo 2DS particles and calculate PSDs

#Steps:
#1) Rawfile processing using OASIS and Image Output as .h5 file. 
#2) Add data paths and probes settings to GetFlightInfo2DS().
#3) Info2DS = GetFlightInfo2DS()
#4) BatchFindStereo(Info2DS,FlightNumberStr), find stereo particles and create PSDs using stereo and traditional.
#5) HybridStereoProcessing(Info2DS,FlightNumberStr), combine traditional and stereo psds for all files in folder.

#v1.0 - 25/08/20
# original. Duplicates the procedures from 'Colocation2DSV2.py'

#v2.0 - 10/09/20
# restructured and simplified FindParticlesOnBothChannels
# save stereo particle diagnostic graphs

#v2.3 30/11/20
# Corrected bug in Stereo svol calculation
# Save IAT and colocation plots at 10 minute intervals
# Create hybrid stereo-standard PSDs with a threshold of 300um to switch between them.
# Derive PSD flags from array element histograms.

#v2.4 6/1/21
# added more flights to data paths
# add diagnostics for the difference in y diameter between the two channels
# Dy filtering can be applied to stereo particle psds
# Flag changed to match amof data format
# Corrected a bug that meant edge particles weren't being rejected correctly.
# Version used for OShea et al AMT 2020 paper

#v2.5 11/2/21
# Probe settings passed to FindParticlesOnBothChannelsV2() functions via GetFlightInfo2DS()
# All files saved to Info2DS[FlightNumberStr, 'Path2DSsave']
# Added thresholds to IAT and colocation time histograms


import datetime
import numpy as np
import pandas as pd
#from netCDF4 import Dataset
#import math
#import bisect
import h5py
#import scipy.io as sio
import matplotlib.pyplot as plt  
import matplotlib.dates as mdates
import os
import netCDF4
#from scipy import stats
from scipy.optimize import curve_fit
from scipy import stats
from FlightInfo2DS import GetFlightInfo2DS

#_________________________________________________________________________________
#Loop through oasis .h5 files

def BatchFindStereo(Info2DS, FlightNumberStr):
    
    print(FlightNumberStr)
    #FlightNumberStr = 'C172'
    Path2DS = Info2DS[FlightNumberStr, 'Path2DS']
    PathSave = Info2DS[FlightNumberStr, 'Path2DSsave']
    FlightDate = Info2DS[FlightNumberStr,'FlightDate']
    
    if not os.path.exists(PathSave):
        os.makedirs(PathSave)
    
    for filena in os.listdir(Path2DS):
        if filena.endswith(".h5") and filena.startswith('base'):
            print(filena)
            FindStereo(filena,Info2DS,FlightNumberStr)
            Flag2DS(Path2DS,PathSave,'Export_'+filena,filena,  FlightDate)


#_________________________________________________________________________________
        
# Select particles with time seperation below ColocationThreshold on the other channel

# Particles ColocationDelta < ColocationThreshold are classed as stereo particles
# Additionally the difference in y diameter must be less than ThresholdDeltaDiameterY um

def FindStereo(filena, Info2DS,FlightNumberStr ):
    SaveFile =1
    
    Path2DS = Info2DS[FlightNumberStr, 'Path2DS']
    PathSave = Info2DS[FlightNumberStr, 'Path2DSsave']
    PixelSize = Info2DS['PixelSize']
    MeanXYFlag =Info2DS[FlightNumberStr,'MeanXYFlag'] # 1= mean xy, 0= max
    BiggestParticle = Info2DS[FlightNumberStr, 'BiggestParticle'] # 0 = BBox, 1 = largest particle
    Streak = Info2DS['Streak']
    MaxAspectRatio = Info2DS['MaxAspectRatio']
    IAT_threshold = Info2DS[FlightNumberStr, 'IAT_threshold']
    ColocationThreshold = Info2DS[FlightNumberStr, 'ColocationThreshold']
    ThresholdDeltaDiameterY = Info2DS[FlightNumberStr, 'ThresholdDeltaDiameterY']

    #load OASIS stats
    Data_h5 = h5py.File(Path2DS+filena, 'r')              
    HeaderMatrixWv=np.array(Data_h5['HeaderMatrixWv'])
    ParticleTimesWv=np.array(Data_h5['ParticleTimesWv'])
    ParticleStatsWv=np.array(Data_h5['ParticleStatsWv'])
    Data_h5.close()
    ImageID = ParticleTimesWv[:,3]
    BufferNumber = ParticleTimesWv[:,6]
    #Get particle buffer time in seconds from midnight. 3 decimal places
    ParticleBufferTimeS= np.zeros(len(BufferNumber))*np.nan
    for i in range(len(BufferNumber)):
        if BufferNumber[i] == 4294967295 :
            ParticleBufferTimeS[i] = np.nan
        else:   
            ParticleBufferTimeS[i] = 3600*HeaderMatrixWv[BufferNumber[i],3]+ 60*HeaderMatrixWv[BufferNumber[i],4] + HeaderMatrixWv[BufferNumber[i],5] + (HeaderMatrixWv[BufferNumber[i],6]/1000)
    ParticleTimeSeconds = ParticleTimesWv[:,0] +(ParticleTimesWv[:,1]/1E9)
    Channel = ParticleTimesWv[:,4] 
    
    if len(Channel[Channel == 0 ]) < 20 : return # if less than 100 particles in file don't calculate flag
         
    if BiggestParticle == 0 : # bounding box
        DiameterX=PixelSize+PixelSize*(ParticleStatsWv[:,4]-ParticleStatsWv[:,3]) #x diameter
        DiameterY = PixelSize+PixelSize*(ParticleStatsWv[:,6]-ParticleStatsWv[:,5])
        MIDx = (ParticleStatsWv[:,4] + ParticleStatsWv[:,3])/2
        #MaxDiameter= np.sqrt(DiameterX**2 + DiameterY**2)
        MaxDiameter=np.where(ParticleStatsWv[:,20] <1, PixelSize,PixelSize*ParticleStatsWv[:,20])
        MeanXYDiameter = (DiameterY+DiameterX) /2 
    if BiggestParticle == 1 :  # largest particle
        DiameterX=PixelSize+PixelSize*(ParticleStatsWv[:,10]-ParticleStatsWv[:,9]) #x diameter
        DiameterY = PixelSize+PixelSize*(ParticleStatsWv[:,12]-ParticleStatsWv[:,11])
        MIDx = (ParticleStatsWv[:,10] + ParticleStatsWv[:,9])/2
        #MaxDiameter= np.sqrt(DiameterX**2 + DiameterY**2)
        MaxDiameter=np.where(ParticleStatsWv[:,19] <1, PixelSize,PixelSize*ParticleStatsWv[:,19])
        MeanXYDiameter = (DiameterY+DiameterX) /2 

    SlicesY = PixelSize+PixelSize*(ParticleStatsWv[:,6]-ParticleStatsWv[:,5]) # Same diameterY when using BBox
    SlicesX = PixelSize+PixelSize*(ParticleStatsWv[:,4]-ParticleStatsWv[:,3]) # Same diameterY when using BBox
    StreakFlag =(np.where(np.logical_and(DiameterX == PixelSize, DiameterY >= Streak*PixelSize), 1, 0))  # 1 = streak, 0 = not 
    AspectRatio = DiameterY / DiameterX
    #AspectRatio = SlicesY / SlicesX
    
    # Edge = 0 not touching array edge
    # Edge = 1 touching array edge
    Edge = (np.where(np.logical_or(ParticleStatsWv[:,3]==0 , ParticleStatsWv[:,4]==127), 1, 0))   
    
    # Select channels
    DiameterX_Ch0 = DiameterX[Channel == 0]
    MeanXYDiameter_Ch0 = MeanXYDiameter[Channel == 0]
    Seconds_Ch0 = ParticleTimeSeconds[Channel == 0]
    Edge_Ch0 = Edge[Channel == 0]
    MaxDiameter_Ch0 = MaxDiameter[Channel == 0]
    Streak_Ch0 = StreakFlag[Channel == 0]
    ParticleBufferTimeS_Ch0 = ParticleBufferTimeS[Channel == 0]
    ImageID_Ch0 = ImageID[Channel == 0]
    MIDx_Ch0 = MIDx[Channel == 0]
    AspectRatio_Ch0 = AspectRatio[Channel == 0]
    DiameterY_Ch0 = DiameterY[Channel == 0]
    DiameterX_Ch0 = DiameterX[Channel == 0]
    SlicesY_Ch0 = SlicesY[Channel == 0]
    
    DiameterX_Ch1 = DiameterX[Channel == 1]
    MeanXYDiameter_Ch1 = MeanXYDiameter[Channel == 1]
    Seconds_Ch1 = ParticleTimeSeconds[Channel == 1]
    Edge_Ch1 = Edge[Channel == 1]
    MaxDiameter_Ch1 = MaxDiameter[Channel == 1]
    Streak_Ch1 = StreakFlag[Channel == 1]
    ParticleBufferTimeS_Ch1 = ParticleBufferTimeS[Channel == 1]
    ImageID_Ch1 = ImageID[Channel == 1]
    MIDx_Ch1 = MIDx[Channel == 1]
    AspectRatio_Ch1 = AspectRatio[Channel == 1]
    DiameterY_Ch1 = DiameterY[Channel == 1]
    DiameterX_Ch1 = DiameterX[Channel == 1]
    SlicesY_Ch1 = SlicesY[Channel == 1]

    print('File length '+ str(len(Channel[~np.isnan(Channel)])))
    print('Channel 0 = '+ str(len(Seconds_Ch0[~np.isnan(Seconds_Ch0)])))
    print('Channel 1 = '+ str(len(Seconds_Ch1[~np.isnan(Seconds_Ch1)])))

    #remove streaking
    Seconds_Ch0 = (np.where(np.logical_or(AspectRatio_Ch0 >= MaxAspectRatio, Streak_Ch0 == 1),np.nan,Seconds_Ch0))
    Seconds_Ch1 = (np.where(np.logical_or(AspectRatio_Ch1 >= MaxAspectRatio, Streak_Ch1 == 1),np.nan,Seconds_Ch1))

    DiameterX_Ch0 = DiameterX_Ch0[~np.isnan(Seconds_Ch0)]
    DiameterY_Ch0 = DiameterY_Ch0[~np.isnan(Seconds_Ch0)]
    MeanXYDiameter_Ch0 = MeanXYDiameter_Ch0[~np.isnan(Seconds_Ch0)]
    MaxDiameter_Ch0 = MaxDiameter_Ch0[~np.isnan(Seconds_Ch0)]
    Edge_Ch0 = Edge_Ch0[~np.isnan(Seconds_Ch0)]
    ImageID_Ch0 = ImageID_Ch0[~np.isnan(Seconds_Ch0)] 
    ParticleBufferTimeS_Ch0 = ParticleBufferTimeS_Ch0[~np.isnan(Seconds_Ch0)]
    MIDx_Ch0=MIDx_Ch0[~np.isnan(Seconds_Ch0)]
    SlicesY_Ch0 = SlicesY_Ch0[~np.isnan(Seconds_Ch0)]
    Seconds_Ch0= Seconds_Ch0[~np.isnan(Seconds_Ch0)]
    
    DiameterX_Ch1 = DiameterX_Ch1[~np.isnan(Seconds_Ch1)]
    DiameterY_Ch1 = DiameterY_Ch1[~np.isnan(Seconds_Ch1)]
    MeanXYDiameter_Ch1 = MeanXYDiameter_Ch1[~np.isnan(Seconds_Ch1)]
    MaxDiameter_Ch1 = MaxDiameter_Ch1[~np.isnan(Seconds_Ch1)]
    Edge_Ch1 = Edge_Ch1[~np.isnan(Seconds_Ch1)]
    ImageID_Ch1 = ImageID_Ch1[~np.isnan(Seconds_Ch1)] 
    ParticleBufferTimeS_Ch1 = ParticleBufferTimeS_Ch1[~np.isnan(Seconds_Ch1)]
    MIDx_Ch1=MIDx_Ch1[~np.isnan(Seconds_Ch1)]
    SlicesY_Ch1 = SlicesY_Ch1[~np.isnan(Seconds_Ch1)]
    Seconds_Ch1= Seconds_Ch1[~np.isnan(Seconds_Ch1)]
    
    #Remove shattering, Get IAT times for each channel seprately
    IAT_Ch0 = GetIAT_TimeInS_2DS(Seconds_Ch0)    
    IAT_Ch1 = GetIAT_TimeInS_2DS(Seconds_Ch1)
    
    # Plot IAT for both channels indpendently
    IATHist(PathSave,IAT_Ch0,IAT_Ch1,IAT_threshold,filena)
    # PLot IAT histograms at 10 min intervals
    IATHist_interval(PathSave,IAT_Ch0,IAT_Ch1,Seconds_Ch0, Seconds_Ch1 )
    
    # Remove particles with IAT below IAT_threshold
    Seconds_Ch0[IAT_Ch0<IAT_threshold] =np.nan
    Seconds_Ch1[IAT_Ch1<IAT_threshold] =np.nan
    
    #Remove nans again
    DiameterX_Ch0 = DiameterX_Ch0[~np.isnan(Seconds_Ch0)]
    DiameterY_Ch0 = DiameterY_Ch0[~np.isnan(Seconds_Ch0)]
    MeanXYDiameter_Ch0 = MeanXYDiameter_Ch0[~np.isnan(Seconds_Ch0)]
    MaxDiameter_Ch0 = MaxDiameter_Ch0[~np.isnan(Seconds_Ch0)]
    Edge_Ch0 = Edge_Ch0[~np.isnan(Seconds_Ch0)]
    ImageID_Ch0 = ImageID_Ch0[~np.isnan(Seconds_Ch0)] 
    ParticleBufferTimeS_Ch0 = ParticleBufferTimeS_Ch0[~np.isnan(Seconds_Ch0)]
    MIDx_Ch0=MIDx_Ch0[~np.isnan(Seconds_Ch0)]
    IAT_Ch0 =  IAT_Ch0[~np.isnan(Seconds_Ch0)]
    SlicesY_Ch0 = SlicesY_Ch0[~np.isnan(Seconds_Ch0)]
    Seconds_Ch0= Seconds_Ch0[~np.isnan(Seconds_Ch0)]
    
    DiameterX_Ch1 = DiameterX_Ch1[~np.isnan(Seconds_Ch1)]
    DiameterY_Ch1 = DiameterY_Ch1[~np.isnan(Seconds_Ch1)]
    MeanXYDiameter_Ch1 = MeanXYDiameter_Ch1[~np.isnan(Seconds_Ch1)]
    MaxDiameter_Ch1 = MaxDiameter_Ch1[~np.isnan(Seconds_Ch1)]
    Edge_Ch1 = Edge_Ch1[~np.isnan(Seconds_Ch1)]
    ImageID_Ch1 = ImageID_Ch1[~np.isnan(Seconds_Ch1)] 
    ParticleBufferTimeS_Ch1 = ParticleBufferTimeS_Ch1[~np.isnan(Seconds_Ch1)]
    MIDx_Ch1=MIDx_Ch1[~np.isnan(Seconds_Ch1)]
    IAT_Ch1 =  IAT_Ch1[~np.isnan(Seconds_Ch1)]
    SlicesY_Ch1 = SlicesY_Ch1[~np.isnan(Seconds_Ch1)]
    Seconds_Ch1= Seconds_Ch1[~np.isnan(Seconds_Ch1)]
    
    # Look for colocated particles. 
    ChTimeDelta_high= np.zeros(len(Seconds_Ch1))*np.nan
    ChIDX_high= np.zeros(len(Seconds_Ch1))*np.nan
    ChTimeDelta_low= np.zeros(len(Seconds_Ch1))*np.nan
    ChIDX_low= np.zeros(len(Seconds_Ch1))*np.nan
    ChTimeDelta= np.zeros(len(Seconds_Ch1))*np.nan
    ChIDX= np.zeros(len(Seconds_Ch1))*np.nan
    
    # Go element by element in channel1. looking for the closest time match in channel0
    ChIDX_high = np.searchsorted(Seconds_Ch0, Seconds_Ch1, side="left") # side="left" means	Seconds_CH0[i-1] < Seconds_CH1 <= Seconds_CH0[i]
    ChIDX_high[ChIDX_high>=len(Seconds_Ch0)] = len(Seconds_Ch0)-1 # if outside the array
    #check whether i or i-1 gives lowest TimeDelta
    ChTimeDelta_high=np.absolute(Seconds_Ch1 - Seconds_Ch0[ChIDX_high])
    ChIDX_low = ChIDX_high - 1
    ChIDX_low[ChIDX_low<0] = 0
    ChTimeDelta_low = np.absolute(Seconds_Ch1 - Seconds_Ch0[ChIDX_low])
    ChIDX = (np.where(ChTimeDelta_low < ChTimeDelta_high,ChIDX_low, ChIDX_high)) # Index with lowest TimeDelta
    ChTimeDelta = (np.where(ChTimeDelta_low < ChTimeDelta_high,ChTimeDelta_low, ChTimeDelta_high)) # Lowest TimeDelta

    print('Number stereo / not stereo')
    print(len(ChTimeDelta[ChTimeDelta < ColocationThreshold])/len(ChTimeDelta[ChTimeDelta > ColocationThreshold]))
     
    # Plot colocation time histogram
    ColocationTimeHist(ChTimeDelta,ColocationThreshold,PathSave,filena)
    #Colocation histograms for 10 minute intervals
    ColocationTimeHist_interval(ChTimeDelta,Seconds_Ch1,PathSave)

    #Select colocated particle stats
    ColocationIDX = (ChIDX[ChTimeDelta <= ColocationThreshold]).astype(int) # Indexes are for channel 0 (same length as channel 1)
    ColocationDelta = ChTimeDelta[ChTimeDelta < ColocationThreshold] 
    #ColocationParticleTime_CH1 = ParticleTime_CH1[ChTimeDelta < ColocationThreshold]
    ColocationDiameterY_Ch1 = DiameterY_Ch1[ChTimeDelta < ColocationThreshold] 
    ColocationDiameterX_Ch1 = DiameterX_Ch1[ChTimeDelta < ColocationThreshold]
    ColocationMeanXYDiameter_Ch1 = MeanXYDiameter_Ch1[ChTimeDelta < ColocationThreshold]
    ColocationMaxDiameter_Ch1 = MaxDiameter_Ch1[ChTimeDelta < ColocationThreshold]
    ColocationSecondsCh1 = Seconds_Ch1[ChTimeDelta < ColocationThreshold]
    ColocationEdgeCh1 = Edge_Ch1[ChTimeDelta < ColocationThreshold]
    ColocationImageID_Ch1 = ImageID_Ch1[ChTimeDelta < ColocationThreshold]
    ColocationParticleBufferTimeS_Ch1 = ParticleBufferTimeS_Ch1[ChTimeDelta < ColocationThreshold]  
    ColocationMIDx_Ch1=MIDx_Ch1[ChTimeDelta < ColocationThreshold]
    ColocationSlicesY_Ch1 = SlicesY_Ch1[ChTimeDelta < ColocationThreshold]
    
    ColocationDiameterY_Ch0 = DiameterY_Ch0[ColocationIDX]
    ColocationDiameterX_Ch0 = DiameterX_Ch0[ColocationIDX]
    ColocationSecondsCh0 = Seconds_Ch0[ColocationIDX]
    ColocationMeanXYDiameter_Ch0= MeanXYDiameter_Ch0[ColocationIDX]
    ColocationMaxDiameter_Ch0 = MaxDiameter_Ch0[ColocationIDX]
    ColocationEdgeCh0 = Edge_Ch0[ColocationIDX]
    ColocationImageID_Ch0 = ImageID_Ch0[ColocationIDX]
    ColocationParticleBufferTimeS_Ch0 = ParticleBufferTimeS_Ch0[ColocationIDX]
    ColocationMIDx_Ch0=MIDx_Ch0[ColocationIDX]
    ColocationSlicesY_Ch0 = SlicesY_Ch0[ColocationIDX]
      
    # select inter arrival time for only colocated particles 
    #IAT_Ch0_colocation = IAT_Ch0[ColocationIDX]
    #IAT_Ch1_colocation = IAT_Ch1[ChTimeDelta < ColocationThreshold]

    DeltaDiameterY_hist(ColocationDiameterY_Ch0, ColocationDiameterY_Ch1, ColocationDelta, ColocationThreshold, PathSave,filena)
  
    # Plot CH1 vs CH0 diameter for stereo particles
    ChannelDComparison(ColocationDiameterY_Ch0,ColocationDiameterY_Ch1,ColocationMeanXYDiameter_Ch0,ColocationMeanXYDiameter_Ch1,
                           IAT_threshold,ColocationThreshold,ColocationSlicesY_Ch0,ColocationSlicesY_Ch1,ColocationDelta,
                            PathSave,filena)
        
    # Plot Midx vs Midx for colocated particles 
    PlotArrayPositionSmallParticles(PathSave,filena, ColocationMeanXYDiameter_Ch0, ColocationMeanXYDiameter_Ch1, ColocationMIDx_Ch1, ColocationMIDx_Ch0)

    #  Particles that are below the colocation threshold with two particles
    RepeatIDs_CH0=[x for x, y in zip(ColocationImageID_Ch0, ColocationImageID_Ch0[1:]) if x>=y]
    RepeatIDs_CH1=[x for x, y in zip(ColocationImageID_Ch1, ColocationImageID_Ch1[1:]) if x>=y]
    print('Double co-location matches ch0 = '+ str(len(RepeatIDs_CH0)))
    print('Double co-location matches ch1 = '+ str(len(RepeatIDs_CH1)))
    
    # save matching particles
    if SaveFile == 1 : 
        SavePath= PathSave+'Colocate_'+filena
        file = h5py.File(SavePath, 'w')
        file.create_dataset('ColocationSecondsCh1', data=ColocationSecondsCh1)
        file.create_dataset('ColocationSecondsCh0', data=ColocationSecondsCh0)
        file.create_dataset('ColocationDelta', data=ColocationDelta)
        file.create_dataset('ColocationMeanXYDiameter_Ch1', data=ColocationMeanXYDiameter_Ch1)
        file.create_dataset('ColocationMeanXYDiameter_Ch0', data=ColocationMeanXYDiameter_Ch0)
        file.create_dataset('ColocationDiameterY_Ch0', data=ColocationDiameterY_Ch0)
        file.create_dataset('ColocationDiameterY_Ch1', data=ColocationDiameterY_Ch1)
        file.create_dataset('ColocationMaxDiameter_Ch1', data=ColocationMaxDiameter_Ch1)
        file.create_dataset('ColocationMaxDiameter_Ch0', data=ColocationMaxDiameter_Ch0)
        file.create_dataset('ColocationEdgeCh0', data=ColocationEdgeCh0)
        file.create_dataset('ColocationEdgeCh1', data=ColocationEdgeCh1)
        file.create_dataset('ColocationParticleBufferTimeS_Ch0', data=ColocationParticleBufferTimeS_Ch0)
        file.create_dataset('ColocationParticleBufferTimeS_Ch1', data=ColocationParticleBufferTimeS_Ch1)
        file.create_dataset('ColocationImageID_Ch0', data=ColocationImageID_Ch0)
        file.create_dataset('ColocationImageID_Ch1', data=ColocationImageID_Ch1)        
        file.create_dataset('ColocationSlicesY_Ch0', data=ColocationSlicesY_Ch0)
        file.create_dataset('ColocationSlicesY_Ch1', data=ColocationSlicesY_Ch1)
        file.create_dataset('IAT_threshold', data=IAT_threshold)
        file.create_dataset('ColocationThreshold', data=ColocationThreshold)
        file.close()
        
        #Create PSDs using all data and just colocated data
    
        #Diameter y filter for stereo particles 
        if ThresholdDeltaDiameterY == - 1: 
            Idx = np.where(ColocationSlicesY_Ch0==ColocationSlicesY_Ch0) # 
        else : 
            MinDiameterY = np.minimum(ColocationSlicesY_Ch0, ColocationSlicesY_Ch1)
            MaxDiameterY = np.maximum(ColocationSlicesY_Ch0, ColocationSlicesY_Ch1)
            Idx = np.where(MaxDiameterY<=(MinDiameterY+ThresholdDeltaDiameterY)) #Select colocation indexes that meet size criteria
            print('Fraction of stereo above Dy threshold ='+str(len(Idx)/len(MinDiameterY)))
        
        if len(ColocationSecondsCh1[Idx]) > 0 and len(ColocationSecondsCh0[Idx]) > 0 :
            if MeanXYFlag == 1: # mean xy
                PSD_Colocate_1hzV2(Info2DS,FlightNumberStr,'dNdD_L_',Seconds_Ch0, DiameterX_Ch0, MeanXYDiameter_Ch0, Edge_Ch0, 
                                   Seconds_Ch1, DiameterX_Ch1, MeanXYDiameter_Ch1, Edge_Ch1,SaveFile,filena,0)
                PSD_Colocate_1hzV2(Info2DS,FlightNumberStr,'dNdD_L_Colocate_',ColocationSecondsCh0[Idx], ColocationDiameterX_Ch0[Idx], ColocationMeanXYDiameter_Ch0[Idx], ColocationEdgeCh0[Idx], 
                                   ColocationSecondsCh1[Idx], ColocationDiameterX_Ch1[Idx], ColocationMeanXYDiameter_Ch1[Idx], ColocationEdgeCh1[Idx],SaveFile,filena,1)
            
            if MeanXYFlag == 0: # max
                PSD_Colocate_1hzV2(Info2DS,FlightNumberStr,'dNdD_L_',Seconds_Ch0, DiameterX_Ch0, MaxDiameter_Ch0, Edge_Ch0, 
                                   Seconds_Ch1, DiameterX_Ch1, MaxDiameter_Ch1, Edge_Ch1,SaveFile,filena,0)
                PSD_Colocate_1hzV2(Info2DS,FlightNumberStr,'dNdD_L_Colocate_',ColocationSecondsCh0[Idx], ColocationDiameterX_Ch0[Idx], ColocationMaxDiameter_Ch0[Idx], ColocationEdgeCh0[Idx], 
                                   ColocationSecondsCh1[Idx], ColocationDiameterX_Ch1[Idx], ColocationMaxDiameter_Ch1[Idx], ColocationEdgeCh1[Idx],SaveFile,filena,1)
                



#__________________________________________________________________________________
 

#Create arrays with dNdD_L_CH0 and dNdD_L_CH1. 
#ColocatedFlag =0  Lawson et al sample volume 
#ColocatedFlag =1  colocation sample volume 

def PSD_Colocate_1hzV2(Info2DS,FlightNumberStr,SavePrefix,SecondsCh0, DiameterXCh0, DiameterCh0, Edge_Ch0, SecondsCh1, DiameterXCh1, DiameterCh1,Edge_Ch1,SaveFile,filena,ColocatedFlag) :
    
    #Probe settings
    Lambda = Info2DS['Lambda_um'] # in um 
    PixelSize = Info2DS['PixelSize']
    ArrayElements = Info2DS['ArrayElements']
    ArmSep = Info2DS[FlightNumberStr,'ArmSep'] 
    c = Info2DS['c']
    TAS = Info2DS[FlightNumberStr,'TAS'] #m/s
    Path2DSsave = Info2DS[FlightNumberStr, 'Path2DSsave']
    
    #Set up size bins
    SizeBinsEdge=np.linspace((PixelSize*0.5),(PixelSize*ArrayElements+PixelSize*0.5),num=ArrayElements+1)
    SizeBinsMid=(SizeBinsEdge[1:]+SizeBinsEdge[0:-1]) / 2

    #Set up time bins
    #Startime=int(np.minimum(SecondsCh0[0], SecondsCh1[0]))
    #Endtime=int(np.maximum(SecondsCh0[-1], SecondsCh1[-1]))    
    
    Startime=int(np.minimum(min(SecondsCh0), min(SecondsCh1)))
    Endtime=int(np.maximum(max(SecondsCh0), max(SecondsCh1)))
    
    if Endtime > 172800:
        Endtime = 172800
    TimeBinsEdge=np.arange(Startime-0.5,Endtime+0.5,1)    
    TimeBinsMid=(TimeBinsEdge[1:]+TimeBinsEdge[0:-1]) / 2 
    
  
    # Calculate PSDs  
    if ColocatedFlag == 1 :# colocated particles psd
        
        #Channel0
        SecondsCh0 = SecondsCh0[Edge_Ch0 == 0]
        DiameterXCh0 = DiameterXCh0[Edge_Ch0 == 0] 
        DiameterCh0 = DiameterCh0[Edge_Ch0 == 0]
        SVolCh0_L_s = Colocated_sVol_arrayV2(DiameterXCh0,DiameterCh0,Lambda,ArmSep,ArrayElements,PixelSize,TAS,c)
        InverseSVol_L_s =1/SVolCh0_L_s # s L-1
        Counts_PSD_Ch0,tmp,tmp = np.histogram2d( SecondsCh0, DiameterCh0, bins=[TimeBinsEdge,SizeBinsEdge], weights=None)
        TotalCounts_Ch0 = np.nansum(Counts_PSD_Ch0, axis = 1 )
        dN_L_Ch0,tmp,tmp = np.histogram2d( SecondsCh0, DiameterCh0, bins=[TimeBinsEdge,SizeBinsEdge], weights=InverseSVol_L_s) 
        #TotalN_L_Ch0 = np.nansum(dN_L_Ch0)
        dNdD_L_Ch0 =  dN_L_Ch0/PixelSize # normalise using bin width
        
        
        #Channel1
        SecondsCh1 = SecondsCh1[Edge_Ch1 == 0]
        DiameterXCh1= DiameterXCh1[Edge_Ch1 == 0] 
        DiameterCh1= DiameterCh1[Edge_Ch1 == 0]
        SVolCh1_L_s = Colocated_sVol_arrayV2(DiameterXCh1,DiameterCh1,Lambda,ArmSep,ArrayElements,PixelSize,TAS,c)
        InverseSVol_L_s =1/SVolCh1_L_s # s L-1
        Counts_PSD_Ch1,tmp,tmp = np.histogram2d( SecondsCh1, DiameterCh1, bins=[TimeBinsEdge,SizeBinsEdge], weights=None)
        TotalCounts_Ch1 = np.nansum(Counts_PSD_Ch1, axis = 1 )
        dN_L_Ch1,tmp,tmp = np.histogram2d( SecondsCh1, DiameterCh1, bins=[TimeBinsEdge,SizeBinsEdge], weights=InverseSVol_L_s) 
        dNdD_L_Ch1 =  dN_L_Ch1/PixelSize # normalise using bin width

    else : # use lawson sample volume
    
        # Channel 0 
        #Edge rejection
        SecondsCh0 = SecondsCh0[Edge_Ch0 == 0]
        DiameterXCh0 = DiameterXCh0[Edge_Ch0 == 0] 
        DiameterCh0= DiameterCh0[Edge_Ch0 == 0] 
        # Sample volume calculation
        SVolCh0_L_s = Lawson_sVol_2DS(DiameterXCh0,DiameterCh0,Lambda,c,ArmSep,ArrayElements,PixelSize,TAS)    
        InverseSVol_L_s =1/SVolCh0_L_s # s L-1
        # Create 1Hz PSDs
        Counts_PSD_Ch0,tmp,tmp = np.histogram2d( SecondsCh0, DiameterCh0, bins=[TimeBinsEdge,SizeBinsEdge], weights=None)
        TotalCounts_Ch0 = np.nansum(Counts_PSD_Ch0, axis = 1 )
        dN_L_Ch0,tmp,tmp = np.histogram2d( SecondsCh0, DiameterCh0, bins=[TimeBinsEdge,SizeBinsEdge], weights=InverseSVol_L_s) 
        dNdD_L_Ch0 =  dN_L_Ch0/PixelSize # normalise using bin width

        # Channel 1 
        #Edge rejection
        SecondsCh1 = SecondsCh1[Edge_Ch1 == 0]
        DiameterXCh1 = DiameterXCh1[Edge_Ch1 == 0] 
        DiameterCh1 = DiameterCh1[Edge_Ch1 == 0] 
        # Sample volume calculation
        SVolCh1_L_s = Lawson_sVol_2DS(DiameterXCh1,DiameterCh1,Lambda,c,ArmSep,ArrayElements,PixelSize,TAS)    
        InverseSVol_L_s =1/SVolCh1_L_s # s L-1
        # Create 1Hz PSDs
        Counts_PSD_Ch1,tmp,tmp = np.histogram2d( SecondsCh1, DiameterCh1, bins=[TimeBinsEdge,SizeBinsEdge], weights=None)
        TotalCounts_Ch1 = np.nansum(Counts_PSD_Ch1, axis = 1 )
        dN_L_Ch1,tmp,tmp = np.histogram2d( SecondsCh1, DiameterCh1, bins=[TimeBinsEdge,SizeBinsEdge], weights=InverseSVol_L_s) 
        dNdD_L_Ch1=  dN_L_Ch1/PixelSize # normalise using bin width


    if SaveFile == 1 : 
        
        # if ColocatedFlag == 1 :
        #     SavePath= Path2DS+'dNdD_L_Colocate_'+filena
        # else:
        #     SavePath= Path2DS+'dNdD_L_'+filena
        SavePath= Path2DSsave+SavePrefix+filena
        file = h5py.File(SavePath, 'w')
        file.create_dataset('dNdD_L_Ch0', data=dNdD_L_Ch0)
        file.create_dataset('Counts_PSD_Ch0', data=Counts_PSD_Ch0)
        file.create_dataset('dNdD_L_Ch1', data=dNdD_L_Ch1)
        file.create_dataset('Counts_PSD_Ch1', data=Counts_PSD_Ch1)
        file.create_dataset('TimeBinsMid', data=TimeBinsMid)
        file.create_dataset('SizeBinsMid', data=SizeBinsMid)
                
        file.close()


    # if 1==2 : 
    #     test = np.datetime64('1900-01-01 00:00:00')
    #     dt= [test + np.timedelta64(np.int32(x),'s') for x in TimeBinsMid]
   
        
        
    #     #dt = datetime.timedelta(seconds = TimeBinsMid)
    #     fig=plt.figure(figsize=(12,8)) 
    #     date_format = mdates.DateFormatter('%H:%M:%S')
    #     plt.plot(dt,TotalCounts_Ch0, label = 'Ch0')
    #     plt.plot(dt,TotalCounts_Ch1, label = 'Ch1')
    #     plt.gca().xaxis.set_major_formatter(date_format)
    #     plt.yscale('log')
    #     plt.xlabel('Time')
    #     plt.ylabel('Total counts')
    #     plt.legend()
        
    #     if ColocatedFlag == 1 :
    #         plt.title('Colocated particles')
    #         plt.savefig(Path2DSsave+'CountsColocate'+filena[:-3]+'.png',dpi=200)
            
    #     else : 
    #         plt.title('All particles')
    #         plt.savefig(Path2DSsave+'Counts'+filena[:-3]+'.png',dpi=200)
            
    #     plt.close(fig)
        

##________________________________________________________________________________-

# Calculate sample volume L/s of co-located particles with edge rejection

def Colocated_sVol_arrayV2(DiameterX,Diameter,Lambda_um,ArmSep,ArrayElements,PixelSize,TAS,c):
    
    radius=Diameter/2 # um
    DoF=(2*c*radius**2)/Lambda_um ## 2* because there is a +/- in the equation
    DoF/=1000 # return in mm
    ArrayWidth_mm = ArrayElements*PixelSize*0.001 
    
    DoF = np.minimum(DoF,ArrayWidth_mm) # sVol Z dimension
    
    TAS*=1000 # mm/s # SVol Airflow dimension 
    
    EffectiveArrayWidth_mm = (((ArrayElements-1)*PixelSize) - (DiameterX))*0.001 # mm
    
    EffectiveArrayWidth_mm = np.minimum(EffectiveArrayWidth_mm,DoF) # Svol optical array dimension
    
    SVol_mm3_s = TAS * EffectiveArrayWidth_mm * DoF
    SVol_L_s = SVol_mm3_s / 1E6
    
    return SVol_L_s


#_______________________________________________________________________________________  
# gets IAT if times are in seconds. Rather than using ParticleTimesWv
    
def GetIAT_TimeInS_2DS(Seconds):
    IAT=np.zeros(len(Seconds))
    for i in range (len(IAT)-1):
        IAT[i]=IAT_TimeInS_2DS(Seconds,i)
    return IAT

#_______________________________________________________________________________________  


def IAT_TimeInS_2DS(Seconds,row):
    
    t1= -np.inf
    t2= Seconds[row]
    t3= np.inf
	
    if(row>0):
        t1=Seconds[row-1]
	
    if(row<(len(Seconds)-1)):
        t3=Seconds[row+1]
        
    IAT1=t2-t1
    IAT2=t3-t2
	
    return min(IAT1,IAT2)

#__________________________________________________________________________________

#Diameter um
#Lambda um
# ArmSep mm
# ProbeRes um
#TAS m/s

# Sample volume from Lawson et 2006

def Lawson_sVol_2DS(DiameterX,Diameter,Lambda_um,c,ArmSep,ArrayElements,ProbeRes,TAS):
    
    radius=Diameter/2 # um
    DoF=(2*c*(radius**2))/Lambda_um ## 2* because there is a +/- in the equation
    DoF/=1000 # return in mm
    
    DoF = np.minimum(DoF,ArmSep)
        
    ArrayWidth = (((ArrayElements-1)*ProbeRes) - (DiameterX))*0.001 # mm
    #ArrayWidth = (((ArrayElements)*ProbeRes) - (DiameterX))*0.001 # mm
    
    TAS*=1000 # mm/s
    
    sVol = DoF * ArrayWidth * TAS # mm3/s
   
    sVol /= 1E6 #L/s
    
    return sVol#, DoF



#_______________________________________________________________________________________  

# Generate data flags from diode histograms.
# looks at how often each pixel is on. compares across array and between channels

 # 0 = not used, 1= valid data, 2= reduced quality data

def Flag2DS(Path2DS,PathSave,ImageFileName,ParticleFileName,FlightDate):
        
    SaveFile = 1
    ArrayElements = 128 
    Nslices = 100000 # Number of slices to calculate stats over
    #Thresholds for flag
    ThresholdMeanDiff = 0.2 #fraction
    ThresholdMinMean = 0
    ThresholdMaxMean = 10
    ThresholdSdevMean = 0.5
    
    #Load image
    Image_h5 = h5py.File(Path2DS + ImageFileName, 'r')
    ImageTimes=np.array(Image_h5['ImageTimes'][:,0])
    LastTimeIDX = np.argmin(ImageTimes)-1 # last particle IDX. There is some filler at end of file
    ImageTimes = ImageTimes[0:LastTimeIDX+1]
    ImageSlices  =np.array(Image_h5['ImageTimes'][:,1])
    ImageSlices = ImageSlices[0:LastTimeIDX+1]
    
    #Find start position of each image within ImageData
    ImagePosition = np.cumsum(ImageSlices, axis = 0)
    ImagePosition = np.append([0], ImagePosition, axis= 0)
    
    #Position of last image in file
    FileSizeSlices= int(ImagePosition[LastTimeIDX])
    
    # Load particle stats file 
    Data_h5 = h5py.File(Path2DS+ ParticleFileName, 'r')             
    HeaderMatrixWv=np.array(Data_h5['HeaderMatrixWv'])
    ParticleTimesWv=np.array(Data_h5['ParticleTimesWv'])
    ParticleStatsWv=np.array(Data_h5['ParticleStatsWv'])
    Data_h5.close()
    
    Channel = np.array(ParticleTimesWv[:,4]).astype(np.float32)
    Channel[Channel>1] = np.nan
    Channel_Image= np.zeros(FileSizeSlices)*np.nan # the channel of each image slice
    
    # if less than 10 particles in file don't calculate flag
    if len(Channel[Channel == 0 ]) < 10 : return 
   
    # the channel of each image slice
    for i in range(len(ImagePosition)-1): 
        Channel_Image[int(ImagePosition[i]):int(ImagePosition[i+1])] = Channel[i]
    
    HistBinsEdge = np.linspace(-0.5,ArrayElements-0.5,num=ArrayElements+1, endpoint=True)
    HistBinsMid = (HistBinsEdge[:-1] + HistBinsEdge[1:]) / 2
    
    # Image idx to calculate stats over
    StatsSteps = np.append(np.arange(0,FileSizeSlices , Nslices),[FileSizeSlices], axis = 0) 

    MeanCh0 =np.zeros(len(StatsSteps)-1)*np.nan
    MaxCh0 =np.zeros(len(StatsSteps)-1)*np.nan
    MinCh0 = np.zeros(len(StatsSteps)-1)*np.nan
    SdevCh0 =np.zeros(len(StatsSteps)-1)*np.nan
    MedianCh0 =np.zeros(len(StatsSteps)-1)*np.nan
    Pc75Ch0 =np.zeros(len(StatsSteps)-1)*np.nan
    Pc25Ch0 =np.zeros(len(StatsSteps)-1)*np.nan
    HistogramsCh0 = np.zeros([len(StatsSteps)-1, ArrayElements ])*np.nan
    
    MeanCh1 =np.zeros(len(StatsSteps)-1)*np.nan
    MaxCh1 =np.zeros(len(StatsSteps)-1)*np.nan
    MinCh1 = np.zeros(len(StatsSteps)-1)*np.nan
    SdevCh1 =np.zeros(len(StatsSteps)-1)*np.nan
    MedianCh1 =np.zeros(len(StatsSteps)-1)*np.nan
    Pc75Ch1 =np.zeros(len(StatsSteps)-1)*np.nan
    Pc25Ch1 =np.zeros(len(StatsSteps)-1)*np.nan
    HistogramsCh1 = np.zeros([len(StatsSteps)-1, ArrayElements ])*np.nan

    SliceMidTime=np.zeros(len(StatsSteps)-1)*np.nan
    SliceStartTime = np.zeros(len(StatsSteps)-1)*np.nan
    SliceEndTime = np.zeros(len(StatsSteps)-1)*np.nan
    
    for i, x in enumerate(StatsSteps[:-1]):  
        #select image slices to calculate stats over
        Image = np.array(Image_h5['ImageData'][:,StatsSteps[i]:StatsSteps[i+1]])
        # whether each slice corresponds to channel or channel 1
        Channel_ImageSlices = Channel_Image[StatsSteps[i]:StatsSteps[i+1]]
        Image[Image == 0 ] = 1
        Image[Image == 255 ] = 0
        #Number of times each pixel is on
        Image_histCh0 = np.nansum(Image[:,Channel_ImageSlices==0], axis= 1)
        Image_histCh1 = np.nansum(Image[:,Channel_ImageSlices==1], axis= 1)
    #
        HistogramsCh0[i,:] = Image_histCh0
        MedianCh0[i] = np.nanmedian(Image_histCh0)
        MeanCh0[i]= np.nanmean(Image_histCh0)
        SdevCh0[i]=np.nanstd(Image_histCh0)#/np.sqrt(np.nansum(Image1_hist))  
        MaxCh0[i] = np.nanmax(Image_histCh0)
        MinCh0[i] = np.nanmin(Image_histCh0)
        Pc75Ch0[i] = np.nanpercentile(Image_histCh0, 75) 
        Pc25Ch0[i] = np.nanpercentile(Image_histCh0, 25)
        
        HistogramsCh1[i,:] = Image_histCh1
        MedianCh1[i] = np.nanmedian(Image_histCh1)
        MeanCh1[i]= np.nanmean(Image_histCh1)
        SdevCh1[i]=np.nanstd(Image_histCh1)#/np.sqrt(np.nansum(Image1_hist))  
        MaxCh1[i] = np.nanmax(Image_histCh1)
        MinCh1[i] = np.nanmin(Image_histCh1)
        Pc75Ch1[i] = np.nanpercentile(Image_histCh1, 75) 
        Pc25Ch1[i] = np.nanpercentile(Image_histCh1, 25)
        
        StartIDX = np.searchsorted(ImagePosition,StatsSteps[i])  # -1 since ImagePosition now has 1 more point than ImageTime
        EndIDX = np.searchsorted(ImagePosition,StatsSteps[i+1]-1)
        if EndIDX > LastTimeIDX:
            EndIDX = LastTimeIDX
        
        SliceMidTime[i] =np.nanmean(ImageTimes[StartIDX:EndIDX])
        SliceStartTime[i] =ImageTimes[StartIDX]
        SliceEndTime[i] =ImageTimes[EndIDX]
    Image_h5.close()
    
    SliceMidDateTime= [FlightDate + np.timedelta64(np.int32(x),'s') for x in SliceMidTime]

    Max_MedianCH0 =  MaxCh0 / MedianCh0
    Min_MedianCH0 =  MinCh0 / MedianCh0
    Max_MeanCH0 =  MaxCh0 / MeanCh0
    Min_MeanCH0 =  MinCh0 / MeanCh0
    IQRCh0 = Pc75Ch0 - Pc25Ch0
    Sdev_meanCH0 = SdevCh0/ MeanCh0
    
    Max_MedianCH1 =  MaxCh1 / MedianCh1 
    Min_MedianCH1 =  MinCh1 / MedianCh1
    Max_MeanCH1 =  MaxCh1 / MeanCh1
    Min_MeanCH1 =  MinCh1 / MeanCh1
    IQRCh1 = Pc75Ch1 - Pc25Ch1
    Sdev_meanCH1 = SdevCh1/ MeanCh1
    
    MeanDiff = (MeanCh1 - MeanCh0) / ((MeanCh1 + MeanCh0)/2)
    
    fig=plt.figure(figsize=(12,12)) 
    date_format = mdates.DateFormatter('%H:%M:%S')
    ax= plt.subplot(4,1,1)
    plt.plot(SliceMidDateTime, MeanCh0,color = 'r', label = 'Ch0')
    plt.plot(SliceMidDateTime, MeanCh1,color = 'b', label = 'Ch1')
    #plt.gca().xaxis.set_major_formatter(date_format)
    plt.ylabel('Mean')
    plt.legend()
    
    ax2 = ax.twinx() 
    color = 'g'
    #ax2.set_ylabel('Difference', color=color)  # we already handled the x-label with ax1
    plt.ylabel('Difference', color=color)  # we already handled the x-label with ax1
    plt.ylim([-1,1])
    plt.plot(SliceMidDateTime,MeanDiff, color=color)
    plt.axhline(-0.2,linestyle = ':')
    plt.axhline(0.2,linestyle = ':')
    plt.axhline(0,linestyle = ':')
    plt.tick_params(axis='y', labelcolor=color)
    plt.gca().xaxis.set_major_formatter(date_format)
    
    #plt.axhline(0,linestyle = ':')
    #ax2.tick_params(axis='y', labelcolor=color)
    #plt.xlim([np.nanmin(SliceMidDateTime), np.nanmax(SliceMidDateTime)])
    
    plt.subplot(4,1,2)
    plt.plot(SliceMidDateTime, Sdev_meanCH0,color = 'r', label = 'Ch0')
    plt.plot(SliceMidDateTime, Sdev_meanCH1,color = 'b', label = 'Ch1')
    plt.ylabel('Sdev / mean')
    #plt.yscale('log')
    plt.ylim([0.01,1])
    plt.axhline(0.5,linestyle = ':')
    plt.gca().xaxis.set_major_formatter(date_format)
    
    plt.subplot(4,1,3)
    plt.plot(SliceMidDateTime, Max_MeanCH0,color = 'r', label = 'Ch0')
    plt.plot(SliceMidDateTime, Max_MeanCH1,color = 'b', label = 'Ch1')
    plt.ylabel('Max / mean')
    #plt.yscale('log')
    plt.ylim([0,2])
    plt.axhline(1.5,linestyle = ':')
    plt.axhline(1,linestyle = ':')
    
    plt.gca().xaxis.set_major_formatter(date_format)
    
    plt.subplot(4,1,4)
    plt.plot(SliceMidDateTime, Min_MeanCH0,color = 'r', label = 'Ch0')
    plt.plot(SliceMidDateTime, Min_MeanCH1,color = 'b', label = 'Ch1')
    plt.ylabel('Min / Mean')
    #plt.yscale('log')
    plt.ylim([0,2])
    plt.axhline(0.5,linestyle = ':')
    plt.axhline(1,linestyle = ':')
    plt.gca().xaxis.set_major_formatter(date_format)
    
    #plt.show()
    
    # Data flags 
    
    # 0 = not used, 1= valid data, 2= reduced quality data
    FlagCh0 = np.where(np.logical_and.reduce((np.absolute(MeanDiff) < ThresholdMeanDiff, Sdev_meanCH0 < ThresholdSdevMean,Max_MeanCH0 < ThresholdMaxMean, Min_MeanCH0 > ThresholdMinMean)),1,2)
    FlagCh1 = np.where(np.logical_and.reduce((np.absolute(MeanDiff) < ThresholdMeanDiff, Sdev_meanCH1 < ThresholdSdevMean,Max_MeanCH1 < ThresholdMaxMean, Min_MeanCH1 > ThresholdMinMean)),1,2)
                    
    
    if SaveFile == 1 : 
        plt.savefig(PathSave+'DiodeStats'+ParticleFileName[:-3]+'.png',dpi=200)
        
        file = h5py.File(PathSave+'Flag'+ParticleFileName, 'w')
        file.create_dataset('FlagCh0', data=FlagCh0)
        file.create_dataset('FlagCh1', data=FlagCh1)
        file.create_dataset('SliceMidTime', data = SliceMidTime)
        file.create_dataset('SliceStartTime', data = SliceStartTime)
        file.create_dataset('SliceEndTime', data = SliceEndTime)
        file.close()

    #ParticleFileName
    
    plt.close(fig)
    
#_______________________________________________________________________________________  

def HybridStereoProcessing(Info2DS,FlightNumberStr, FillValue):
    
    #FlightNumberStr = 'C174'
    #Path2DS = Info2DS[FlightNumberStr, 'Path2DS']
    PathSave = Info2DS[FlightNumberStr, 'Path2DSsave']
    FlightDate = Info2DS[FlightNumberStr,'FlightDate']
    ThresholdSize =  Info2DS[FlightNumberStr,'ThresholdSize']
    
    files = [F for F in os.listdir(PathSave) if F.endswith(".h5") and F.startswith('Flagbase')]

    # go though each 2ds file and calculate stereo psds
    for  i, FlagFileName in enumerate(files):
        #if FlagFileName.endswith(".h5") and FlagFileName.startswith('Flagbase'):
        FileName = FlagFileName[4:]
        TmpFlagHybrid, TmpCounts_PSD_Hybrid, TmpdNdD_L_Hybrid, TmpTimeMidBins_s_colocation, PSD_SizeMid = HybridStereoPSD(PathSave,FileName,ThresholdSize)            

        if i ==0 : 
            FlightDateStr = FileName[4:10]
            FlagHybrid= TmpFlagHybrid
            TimeMidBins_s_colocation = TmpTimeMidBins_s_colocation
            Counts_PSD_Hybrid = TmpCounts_PSD_Hybrid
            dNdD_L_Hybrid = TmpdNdD_L_Hybrid
        else : 
            FlagHybrid = np.append(FlagHybrid,TmpFlagHybrid,axis = 0)
            TimeMidBins_s_colocation= np.append(TimeMidBins_s_colocation,TmpTimeMidBins_s_colocation,axis = 0)
            Counts_PSD_Hybrid= np.append(Counts_PSD_Hybrid,TmpCounts_PSD_Hybrid,axis = 0)
            dNdD_L_Hybrid=np.append(dNdD_L_Hybrid,TmpdNdD_L_Hybrid,axis = 0)
    
    #Fill gaps between files. 
    StartPoint = np.nanmin(TimeMidBins_s_colocation)
    EndPoint = np.nanmax(TimeMidBins_s_colocation)
    
    MaxNumberHours  = 72
    
    if EndPoint - StartPoint > 3600* MaxNumberHours : 
        print('File longer than max number of hours allowed ')
        print(StartPoint)
        print(EndPoint)
        return 
    
    TimeMid = np.arange(StartPoint,EndPoint+1, step = 1 )
    TimeIdx = np.searchsorted(TimeMid,TimeMidBins_s_colocation )
    
    dNdD_L_Output = np.zeros([len(TimeMid),128])*np.nan
    Counts_Output = np.zeros([len(TimeMid),128])*np.nan
    Flag_output = np.ones([len(TimeMid)])* 3 # missing data = 3
    
    dNdD_L_Output[TimeIdx,:] = dNdD_L_Hybrid
    Counts_Output[TimeIdx,:] = Counts_PSD_Hybrid
    Flag_output[TimeIdx] = FlagHybrid
    
    fig, ax1 =plt.subplots()
    ax1.plot(TimeMid, np.nanmean(dNdD_L_Output, axis = 1), label='Output')
    ax1.plot(TimeMidBins_s_colocation, np.nanmean(dNdD_L_Hybrid, axis = 1),'+', label='Raw')
    plt.legend()
    
    # Missing data rows to FillValue 
    Flag_output[np.isnan(dNdD_L_Output).any(axis=1)] = 3
    dNdD_L_Output[np.isnan(dNdD_L_Output).any(axis=1)] = FillValue
    Counts_Output[np.isnan(Counts_Output).any(axis=1)] = FillValue 
    
    ax2 = ax1.twinx()
    ax2.plot(TimeMid, Flag_output, label='Flag')
    
    
    return dNdD_L_Output, Counts_Output, Flag_output, TimeMid, PSD_SizeMid
    

    #return FlagHybrid, Counts_PSD_Hybrid, dNdD_L_Hybrid, TimeMidBins_s_colocation, PSD_SizeMid
#_______________________________________________________________________________________  

# Use ThresholdSize to switch between colocate and traditional 2ds data processing
#ThresholdSize = 300


def HybridStereoPSD(FilePath,FileName, ThresholdSize):
    
    #Load colocation PSDs
    Data_h5 = h5py.File(FilePath + 'dNdD_L_Colocate_'+FileName, 'r')
    dNdD_L_CH0_colocation=np.array(Data_h5['dNdD_L_Ch0'])
    dNdD_L_CH1_colocation=np.array(Data_h5['dNdD_L_Ch1'])
    Counts_PSD_CH0_colocation=np.array(Data_h5['Counts_PSD_Ch0'])
    Counts_PSD_CH1_colocation=np.array(Data_h5['Counts_PSD_Ch1'])
    TimeMidBins_s_colocation=np.array(Data_h5['TimeBinsMid'])
    PSD_SizeMid_colocation=np.array(Data_h5['SizeBinsMid'])
    Data_h5.close()
    #Load PSDs standard
    Data_h5 = h5py.File(FilePath + 'dNdD_L_'+FileName, 'r')
    dNdD_L_CH0=np.array(Data_h5['dNdD_L_Ch0'])
    dNdD_L_CH1=np.array(Data_h5['dNdD_L_Ch1'])
    Counts_PSD_CH0=np.array(Data_h5['Counts_PSD_Ch0'])
    Counts_PSD_CH1=np.array(Data_h5['Counts_PSD_Ch1'])
    TimeMidBins_s=np.array(Data_h5['TimeBinsMid'])
    PSD_SizeMid=np.array(Data_h5['SizeBinsMid'])
    Data_h5.close()
    #load flag
    Data_h5 = h5py.File(FilePath + 'Flag'+FileName, 'r')
    FlagCh0=np.array(Data_h5['FlagCh0'])
    FlagCh1=np.array(Data_h5['FlagCh1'])
    SliceMidTime=np.array(Data_h5['SliceMidTime'])
    SliceStartTime=np.array(Data_h5['SliceStartTime'])
    SliceEndTime=np.array(Data_h5['SliceEndTime'])
    Data_h5.close()
    
    ThresholdSizeIdx = np.searchsorted( PSD_SizeMid, ThresholdSize) # Size bin index of threshold size
    TimeStartIdx = np.searchsorted(TimeMidBins_s, TimeMidBins_s_colocation[0]) # traditional should always start before colocation
    TimeEndIdx = 1+np.searchsorted(TimeMidBins_s, TimeMidBins_s_colocation[-1]) # traditional should always end after colocation
        
    #Greater than threshold use average both channels using traditional procesisng
    dNdD_L_Hybrid = (dNdD_L_CH0[TimeStartIdx:TimeEndIdx,:]+dNdD_L_CH1[TimeStartIdx:TimeEndIdx,:])/2
    Counts_PSD_Hybrid = (Counts_PSD_CH0[TimeStartIdx:TimeEndIdx,:]+Counts_PSD_CH1[TimeStartIdx:TimeEndIdx,:])
    
    #Use channel 0 stereo for less than size threshold
    dNdD_L_Hybrid[:,0:ThresholdSizeIdx] = dNdD_L_CH0_colocation[:,0:ThresholdSizeIdx]
    Counts_PSD_Hybrid[:,0:ThresholdSizeIdx] = Counts_PSD_CH0_colocation[:,0:ThresholdSizeIdx]
    
    # TestColocation = np.nanmean(dNdD_L_CH0_colocation, axis = 0)
    # TestTrad = np.nanmean(dNdD_L_CH0, axis = 0)
    # TestHybrid = np.nanmean(dNdD_L_Hybrid, axis = 0) 
    
    # plt.plot(PSD_SizeMid,TestColocation )
    # plt.plot(PSD_SizeMid,TestTrad )
    # plt.plot(PSD_SizeMid,TestHybrid,'o')
    # plt.xscale('log')
    # plt.yscale('log')
    
    #Flag to 1hz timebase
    
    #Put flag on same timebase as PSDs
    FlagStartIdx = np.searchsorted( TimeMidBins_s_colocation, SliceStartTime)
    FlagEndIdx = np.searchsorted( TimeMidBins_s_colocation, SliceEndTime)
    
    FlagHybrid = np.ones(len(TimeMidBins_s_colocation))*2
    
    for x in range(len(FlagStartIdx)):
        FlagHybrid[FlagStartIdx[x]:FlagEndIdx[x]] = np.maximum(FlagCh0[x], FlagCh1[x]) 
           
    return FlagHybrid, Counts_PSD_Hybrid, dNdD_L_Hybrid, TimeMidBins_s_colocation, PSD_SizeMid
    


#_______________________________________________________________________________________  
#
# Plot colocation time histogram

def ColocationTimeHist(ChTimeDelta,ColocationThreshold,PathSave,filena) : 
    
    #fig=plt.figure(figsize=(7,7)) 
    #plt.rcParams.update({'font.size': 12})
    plt.subplot(2,1,2)
    
    #log spaced bins but also include 0 
    ColocationBinsEdge = np.logspace(-7,0,num=100)
    ColocationBinsEdge = np.append([0],ColocationBinsEdge)
    ColocationBinsMid=(ColocationBinsEdge[1:]+ColocationBinsEdge[0:-1]) /2
    ColocationBinsMid[0]=0 
    ColocationHist, ColocationBinsEdge =  np.histogram(ChTimeDelta, ColocationBinsEdge )
    ColocationHist = np.where(ColocationHist == 0, np.nan, ColocationHist )        
    #ColocationHist[ColocationHist==0] =np.nan
    plt.plot(ColocationBinsMid,ColocationHist,'o')
    #plt.hist(ChTimeDelta,HistBins)
    plt.vlines(ColocationThreshold, ymin =0 , ymax= np.nanmax(ColocationHist), color ='k', linestyle ='--')
    plt.xscale('symlog', linthreshx=1E-7)
    plt.xlabel('Co-location time, s')
    plt.ylabel('Counts')
    #plt.title('zeros = '+ str(len(ChTimeDelta[ChTimeDelta == 0]))+', IAT = '+str(IAT_treshold)+'s, Colocation = ' + str(ColocationThreshold)+ 's' )
    #plt.title('IAT = '+str(IAT_treshold)+'s Colocation = ' + str(ColocationThreshold)+ 's' )
    plt.savefig(PathSave+filena[:-3]+'_deltat.png',dpi=200)
    #plt.close(fig)





#_______________________________________________________________________________________  
#Calculate histograms of colocation times at set intervals


def ColocationTimeHist_interval(ChTimeDelta,Seconds_Ch1,PathSave):
    
    Interval = 60
    FitMinX = 5E-6
    ColocationBinsEdge = np.logspace(-7.9,0,num=100)
    ColocationBinsEdge = np.append([0],ColocationBinsEdge)
    
    ColocationBinsMid=(ColocationBinsEdge[1:]+ColocationBinsEdge[0:-1]) /2
    ColocationBinsMid[0]=0 
    
    TimeBinsEdge = np.arange(np.nanmin(Seconds_Ch1), np.nanmax(Seconds_Ch1),  Interval)
    TimeBinsMid=(TimeBinsEdge[1:]+TimeBinsEdge[0:-1]) /2
    if len(TimeBinsEdge) < 2 :
        return

    ColocationHist,tmp,tmp = np.histogram2d(Seconds_Ch1, ChTimeDelta, bins=[TimeBinsEdge,ColocationBinsEdge], weights=None)
    ColocationHist = np.where(ColocationHist == 0, np.nan, ColocationHist )  # remove zeros   
    
    test = datetime.datetime(1,1,1)
    dt = [test +datetime.timedelta(seconds = x) for x in TimeBinsEdge ]
        
    
    FullPathSave = PathSave + 'ColocationHistograms/' 
    if not os.path.exists(FullPathSave):
        os.makedirs(FullPathSave)
   
    for x in range(len(TimeBinsMid)) :
    
        if np.nansum(ColocationHist[x,:]) > 50 : # Only create histogram plot if greater than 50 counts
            fig=plt.figure(figsize=(8,8)) 
            plt.plot(ColocationBinsMid,ColocationHist[x,:],'o')
            #plt.plot(PlotFitx, PlotFity, label= 'Fit')
            #plt.hist(ChTimeDelta,HistBins)
            plt.xscale('symlog', linthreshx=1E-7)
            plt.xlim([0,1])
            plt.xlabel('Co-location time, s')
            plt.ylabel('Counts')
                
            # guassian fit
            # FitMinXIdx = np.searchsorted(ColocationBinsMid,FitMinX)
            # FitX = np.log10(ColocationBinsMid[FitMinXIdx:])
            # FitY = ColocationHist[x,FitMinXIdx:]
            # H, A, x0, sigma = gauss_fit(FitX, FitY)
            # PlotFitx = ColocationBinsMid[1:] # skip colocationtime = 0 
            # PlotFity = gauss(np.log10(PlotFitx), H, A, x0, sigma)
            # plt.plot(PlotFitx, PlotFity, label= 'Log-normal')
            
            plt.title(dt[x].strftime("%H:%M:%S")+' to '+dt[x+1].strftime("%H:%M:%S"))    
            plt.savefig(FullPathSave+'ColocationTime_'+dt[x].strftime("%H-%M-%S")+'.png',dpi=200)
            plt.close(fig)
        
    # ColocationMode = ColocationBinsMid[np.nanargmax(ColocationHist, axis = 1)]   
    # FirstTimeStr = dt[0].strftime("%H-%M-%S")
    # test = np.datetime64('1900-01-01 00:00:00')
    # dt= [test + np.timedelta64(np.int32(x),'s') for x in TimeBinsMid]
    
    # fig=plt.figure(figsize=(12,8)) 
    # date_format = mdates.DateFormatter('%H:%M:%S')
    # plt.plot(dt,ColocationMode,'o')
    # plt.gca().xaxis.set_major_formatter(date_format)
    # plt.ylim([1e-7,1])
    # plt.yscale('log')
    # plt.xlabel('Time')
    # plt.ylabel('Colocation mode, s')
    # plt.savefig(FullPathSave+'ColocationTimeMode_'+FirstTimeStr+'.png',dpi=200)
       
  
#_______________________________________________________________________________________  
 
#Plot and save IAT histogram for both channel   
  
def IATHist(PathSave,IAT_Ch0,IAT_Ch1,IAT_threshold,filena ):
    fig=plt.figure(figsize=(8,12)) 
    plt.rcParams.update({'font.size': 12})
    plt.subplot(2,1,1)
    HistBins = np.logspace(-8,0,num=100)
    IAT_Ch0_hist, tmp= np.histogram(IAT_Ch0,bins=HistBins)
    HistBinsMid = (HistBins[:-1] + HistBins[1:]) / 2
    plt.plot(HistBinsMid,IAT_Ch0_hist,marker='o', linewidth=0, color = 'r',label='Channel 0')       
    IAT_Ch1_hist, tmp= np.histogram(IAT_Ch1,bins=HistBins)
    plt.plot(HistBinsMid,IAT_Ch1_hist,marker='o', linewidth=0, color = 'c',label='Channel 1')
    plt.vlines(IAT_threshold, ymin =0 , ymax= np.nanmax(IAT_Ch1_hist),color ='k', linestyle ='--')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Inter-arrival time, s')
    plt.ylabel('Counts')
    plt.legend()
    #plt.title('IAT = '+str(IAT_treshold)+'s, Colocation = ' + str(ColocationThreshold)+ 's' )
    plt.savefig(PathSave+filena[:-3]+'_IAT.png',dpi=200)
    #plt.close(fig)
    
    
#_______________________________________________________________________________________  

# Plot and save IAT histograms for a given time interval


def IATHist_interval(PathSave,IAT_Ch0,IAT_Ch1,Seconds_Ch0, Seconds_Ch1 ):

    Interval = 60
    IATBinsEdge = np.logspace(-7,0,num=100)
    #IATBinsEdge = np.append([0],IATBinsEdge)
    
    IATBinsMid=(IATBinsEdge[1:]+IATBinsEdge[0:-1]) /2
    #ColocationBinsMid[0]=0 
    
    TimeBinsEdge = np.arange(np.nanmin(Seconds_Ch1), np.nanmax(Seconds_Ch1),  Interval)
    TimeBinsMid=(TimeBinsEdge[1:]+TimeBinsEdge[0:-1]) /2
    if len(TimeBinsEdge) < 2 :
        return

    IATCh0Hist,tmp,tmp = np.histogram2d(Seconds_Ch0, IAT_Ch0, bins=[TimeBinsEdge,IATBinsEdge], weights=None)
    IATCh1Hist,tmp,tmp = np.histogram2d(Seconds_Ch1, IAT_Ch1, bins=[TimeBinsEdge,IATBinsEdge], weights=None)
    
    
    
    #test = np.datetime64('1900-01-01 00:00:00')
    #dt= [test + np.timedelta64(np.int32(x),'s') for x in TimeBinsMid]
   
    test = datetime.datetime(1,1,1)
    dt = [test +datetime.timedelta(seconds = x) for x in TimeBinsEdge ]
        
    
    FullPathSave = PathSave + 'IATHistograms/' 
    if not os.path.exists(FullPathSave):
        os.makedirs(FullPathSave)
   
    
   
    for x in range(len(TimeBinsMid)) :
        #dt = datetime.timedelta(seconds = TimeBinsMid)
        fig=plt.figure(figsize=(8,8)) 
        plt.plot(IATBinsMid,IATCh0Hist[x,:], label='Ch0')
        plt.plot(IATBinsMid,IATCh1Hist[x,:], label='Ch1')
        #plt.hist(ChTimeDelta,HistBins)
        #plt.xscale('symlog', linthreshx=1E-7)
        plt.xscale('log')
        plt.xlabel('IAT, s')
        plt.ylabel('Counts')
        plt.title(dt[x].strftime("%H:%M:%S")+' to '+dt[x+1].strftime("%H:%M:%S"))    
        plt.savefig(FullPathSave+'IAT_'+dt[x].strftime("%H-%M-%S")+'.png',dpi=200)
        plt.close(fig)
    
    
#____________________________________________________________________________________


def gauss(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
#____________________________________________________________________________________

def gauss_fit(x, y):
    valid = ~(np.isnan(x) | np.isnan(y))
    x = x[valid]
    y= y[valid]
    
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    popt, pcov = curve_fit(gauss, x, y, p0=[min(y), max(y), mean, sigma])
    return popt

#____________________________________________________________________________________

def DeltaDiameterY_hist(ColocationDiameterY_Ch0, ColocationDiameterY_Ch1, ColocationDelta, ColocationThreshold, PathSave,filena):

    ColocationBinsEdge = np.arange((-1E-7)/2,1E-5,step=1E-7)
    #ColocationBinsEdge = np.append([0],ColocationBinsEdge)
    ColocationBinsMid=(ColocationBinsEdge[1:]+ColocationBinsEdge[0:-1]) /2
    ColocationBinsMid[0]=0 
    
    # Histogram y size difference between the two channels
    SizeBinEdge = np.arange(-5,1285,step = 10)
    SizeBinMid=(SizeBinEdge[1:]+SizeBinEdge[0:-1]) /2
    DeltaDiameterY = np.absolute(ColocationDiameterY_Ch0 - ColocationDiameterY_Ch1)
    DeltaTimeSizeHist,tmp,tmp = np.histogram2d(ColocationDelta,DeltaDiameterY, bins=[ColocationBinsEdge,SizeBinEdge], weights=None)
    RowSum = np.nansum(DeltaTimeSizeHist, axis = 1)
    
    DeltaSizeHist = np.nansum(DeltaTimeSizeHist, axis = 0) 
    
    DeltaTimeSizeHist /= RowSum[:,None] # normalise DeltaDiameterY histograms
    
    fig=plt.figure(figsize=(12,8))
    plt.rcParams.update({'font.size': 14})
    ax = plt.subplot(2,1,1)
    DeltaSizePdf = DeltaSizeHist/np.nansum(DeltaSizeHist)
    plt.plot(SizeBinMid, DeltaSizePdf)
    #plt.xscale('log')
    plt.xlim([0,300])
    plt.ylim([0,1])
    plt.ylabel('Normalised frequency')
    plt.xlabel('Diameter Y Ch0 - Ch1, μm')
    
    ax2 = ax.twinx()
    color = 'r'
    DeltaSizeCumulative = np.cumsum(DeltaSizePdf)
    plt.plot(SizeBinMid, DeltaSizeCumulative, color=color)
    plt.tick_params(axis='y', labelcolor=color)
    plt.ylabel('Cumulative frequency', color = color)
    plt.ylim([0,1])

    
    plt.subplot(2,1,2)
    plt.pcolormesh(SizeBinEdge,ColocationBinsEdge,DeltaTimeSizeHist, cmap='gist_earth',vmin=0,vmax=1) #The grid orientation follows the standard matrix convention: An array C with shape (nrows, ncolumns) is plotted with the column number as X and the row number as Y.
    plt.xlabel('Diameter Y Ch0 - Ch1, μm')
    plt.xlim([0,300])
    #plt.xscale('log')
    plt.ylabel('Co-location time, s')
    plt.ylim([0,ColocationThreshold])

    #plt.ylabel('D measured, μm')
    cbar=plt.colorbar(orientation='vertical')
    cbar.set_label('Normalised frequency', rotation=270, labelpad=20)
    plt.savefig(PathSave+filena[:-3]+'_DeltaTimeSize.png',dpi=200)    
    plt.close(fig)

#____________________________________________________________________________________

# Great stereo psds applying different particle filters

def Stereo_Dy_filter(Path2DS,filena,StartTime, EndTime,suffix):
            
    #Load particles with ColocationDelta < ColocationThreshold 

    Data_h5 = h5py.File(Path2DS+'Colocate_'+filena, 'r')
    ColocationSecondsCh1=np.array(Data_h5['ColocationSecondsCh1'])
    ColocationSecondsCh0=np.array(Data_h5['ColocationSecondsCh0'])
    ColocationDelta=np.array(Data_h5['ColocationDelta'])
    ColocationMeanXYDiameter_Ch1=np.array(Data_h5['ColocationMeanXYDiameter_Ch1'])
    ColocationMeanXYDiameter_Ch0=np.array(Data_h5['ColocationMeanXYDiameter_Ch0'])
    ColocationDiameterY_Ch0=np.array(Data_h5['ColocationDiameterY_Ch0'])
    ColocationDiameterY_Ch1=np.array(Data_h5['ColocationDiameterY_Ch1'])
    ColocationMaxDiameter_Ch1=np.array(Data_h5['ColocationMaxDiameter_Ch1'])
    ColocationMaxDiameter_Ch0=np.array(Data_h5['ColocationMaxDiameter_Ch0'])
    ColocationEdgeCh0=np.array(Data_h5['ColocationEdgeCh0'])
    ColocationEdgeCh1=np.array(Data_h5['ColocationEdgeCh1'])
    ColocationSlicesY_Ch0 = np.array(Data_h5['ColocationSlicesY_Ch0'])
    ColocationSlicesY_Ch1 = np.array(Data_h5['ColocationSlicesY_Ch1'])
    ColocationParticleBufferTimeS_Ch0 =np.array(Data_h5['ColocationParticleBufferTimeS_Ch0'])
    ColocationParticleBufferTimeS_Ch1 =np.array(Data_h5['ColocationParticleBufferTimeS_Ch1'])
    ColocationImageID_Ch0=np.array(Data_h5['ColocationImageID_Ch0'])
    ColocationImageID_Ch1=np.array(Data_h5['ColocationImageID_Ch1'])
    IAT_threshold=np.array(Data_h5['IAT_threshold'])
    ColocationThreshold=np.array(Data_h5['ColocationThreshold'])
    
    Data_h5.close()
    Path2DS += 'PSDFiltered/'
    
    #Make directory to save 
    if not os.path.exists(Path2DS):
        os.makedirs(Path2DS)
    
    # Dx not in colocation file. Use MeanXY and Dy. ColocationSlicesY maybe different to Dy if using largest image fragment. 
    ColocationDiameterX_Ch0 =  2*ColocationMeanXYDiameter_Ch0 - ColocationDiameterY_Ch0
    ColocationDiameterX_Ch1 =  2*ColocationMeanXYDiameter_Ch1 - ColocationDiameterY_Ch1 
    
    MinDiameterY = np.minimum(ColocationSlicesY_Ch0, ColocationSlicesY_Ch1)
    MaxDiameterY = np.maximum(ColocationSlicesY_Ch0, ColocationSlicesY_Ch1)
    
    # No Dy filter
    PSD_Colocate_1hzV2('dNdD_L_Colocate_NoDy',ColocationSecondsCh0, ColocationDiameterX_Ch0, ColocationMeanXYDiameter_Ch0, ColocationEdgeCh0, ColocationSecondsCh1, ColocationDiameterX_Ch1, ColocationMeanXYDiameter_Ch1, ColocationEdgeCh1,1,Path2DS,filena,1)

    # Flag images using maxD > minD**1.1 + 10
    ThresholdDiameterYExponent = 1.1
    Idx = np.where(MaxDiameterY<(MinDiameterY**ThresholdDiameterYExponent + 10)) #Select colocation indexes that meet size criteria
    PSD_Colocate_1hzV2('dNdD_L_Colocate_Dy^1_1',ColocationSecondsCh0[Idx], ColocationDiameterX_Ch0[Idx], ColocationMeanXYDiameter_Ch0[Idx], ColocationEdgeCh0[Idx], ColocationSecondsCh1[Idx], ColocationDiameterX_Ch1[Idx], ColocationMeanXYDiameter_Ch1[Idx], ColocationEdgeCh1[Idx],1,Path2DS,filena,1)
  
    # maxD <= minD+20
    ThresholdDeltaDiameterY = 20
    Idx = np.where(MaxDiameterY<=(MinDiameterY+ThresholdDeltaDiameterY)) #Select colocation indexes that meet size criteria
    PSD_Colocate_1hzV2('dNdD_L_Colocate_Dy_plus20',ColocationSecondsCh0[Idx], ColocationDiameterX_Ch0[Idx], ColocationMeanXYDiameter_Ch0[Idx], ColocationEdgeCh0[Idx], ColocationSecondsCh1[Idx], ColocationDiameterX_Ch1[Idx], ColocationMeanXYDiameter_Ch1[Idx], ColocationEdgeCh1[Idx],1,Path2DS,filena,1)
  
    # maxD <= minD+40
    ThresholdDeltaDiameterY = 40
    Idx = np.where(MaxDiameterY<=(MinDiameterY+ThresholdDeltaDiameterY)) #Select colocation indexes that meet size criteria
    PSD_Colocate_1hzV2('dNdD_L_Colocate_Dy_plus40',ColocationSecondsCh0[Idx], ColocationDiameterX_Ch0[Idx], ColocationMeanXYDiameter_Ch0[Idx], ColocationEdgeCh0[Idx], ColocationSecondsCh1[Idx], ColocationDiameterX_Ch1[Idx], ColocationMeanXYDiameter_Ch1[Idx], ColocationEdgeCh1[Idx],1,Path2DS,filena,1)
  
    # maxD <= minD+10
    #ThresholdDeltaDiameterY = 10
    #Idx = np.where(MaxDiameterY<=(MinDiameterY+ThresholdDeltaDiameterY)) #Select colocation indexes that meet size criteria
    #PSD_Colocate_1hzV2('dNdD_L_Colocate_Dy_plus10',ColocationSecondsCh0[Idx], ColocationDiameterX_Ch0[Idx], ColocationMeanXYDiameter_Ch0[Idx], ColocationEdgeCh0[Idx], ColocationSecondsCh1[Idx], ColocationDiameterX_Ch1[Idx], ColocationMeanXYDiameter_Ch1[Idx], ColocationEdgeCh1[Idx],1,Path2DS,filena,1)
  
  
    #PLot avg PSDs

    plt.figure(figsize=(8,8))
    FilterNames = ['NoDy','Dy^1_1','Dy_plus20','Dy_plus40']
    
    for Names in FilterNames : 

        Data_h5 = h5py.File(Path2DS + 'dNdD_L_Colocate_'+Names+filena, 'r')
        dNdD_L_CH0_colocation=np.array(Data_h5['dNdD_L_Ch0'])
        dNdD_L_CH1_colocation=np.array(Data_h5['dNdD_L_Ch1'])
        Counts_PSD_CH0_colocation=np.array(Data_h5['Counts_PSD_Ch0'])
        Counts_PSD_CH1_colocation=np.array(Data_h5['Counts_PSD_Ch1'])
        TimeMidBins_s_colocation=np.array(Data_h5['TimeBinsMid'])
        PSD_SizeMid_colocation=np.array(Data_h5['SizeBinsMid'])
        Data_h5.close() 

        if (StartTime == 0) and (EndTime == 0): #Avg over file
            dNdD_L_CH0_avg = np.nanmean(dNdD_L_CH0_colocation, axis = 0)
        else : # avg between StartTime and EndTime 
            StartIdx = np.searchsorted(TimeMidBins_s_colocation,StartTime)
            EndIdx = np.searchsorted(TimeMidBins_s_colocation,EndTime) 
            dNdD_L_CH0_avg = np.nanmean(dNdD_L_CH0_colocation[StartIdx:EndIdx,:], axis = 0)
        
        plt.plot(PSD_SizeMid_colocation, dNdD_L_CH0_avg, label = Names)
        
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel('dN dDp$^{-1}$, L$^{-1}$ μm$^{-1}$')
    plt.xlabel('Diameter, μm')
    plt.legend()
    plt.savefig(Path2DS+filena[:-3]+'_PSD'+suffix+'.png',dpi=200)
    

#____________________________________________________________________________________
# 2DS sample volume used in OASIS version 1.49962

def OASIS_svol(Diameter,ArmSep, ArrayElements,ProbeRes,TAS ):
    
    TAS*=1000
    
    ArrayWidth = ((ArrayElements*ProbeRes) - Diameter)*0.001 #mm
    
    
    z=5130									# 5.13 um-1 is 5130 mm-1
    Diameter/=1000							# Diameter in mm (diamter is passed in microns)
    DoF = np.minimum(z*Diameter*Diameter,ArmSep) #mm
    
    sVol = DoF * ArrayWidth * TAS # mm 
    return sVol

#____________________________________________________________________________________
# Scatter plots of diameter from CH0 and CH1 for particles classified as stereo

def ChannelDComparison(ColocationDiameterY_Ch0,ColocationDiameterY_Ch1,ColocationMeanXYDiameter_Ch0,ColocationMeanXYDiameter_Ch1,
                       IAT_threshold,ColocationThreshold,ColocationSlicesY_Ch0,ColocationSlicesY_Ch1,ColocationDelta,
                        PathSave,filena) : 
    # Plot CH1 vs CH0 diameter
    if (len(ColocationDiameterY_Ch0) > 10):
        fig=plt.figure(figsize=(8,12)) 
        plt.rcParams.update({'font.size': 12})
        plt.subplot(2,1,1)
        plt.plot(ColocationMeanXYDiameter_Ch0,ColocationMeanXYDiameter_Ch1,'o',color='tab:gray',markersize=1)        
        #plt.scatter(ColocationMeanXYDiameter_Ch0,ColocationMeanXYDiameter_Ch1,c= ColocationDelta,markersize=1)        
        
        plt.xlabel('Channel 0 mean, μm')
        plt.ylabel('Channel 1 mean, μm')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim([5,1280])
        plt.ylim([5,1280])
        plt.title('IAT = '+str(IAT_threshold)+'s, Colocation = ' + str(ColocationThreshold)+ 's' )

        plt.subplot(2,1,2)
        #plt.plot(ColocationDiameterY_Ch0,ColocationDiameterY_Ch1,'o',color='tab:gray',markersize=1)        
        plt.scatter(ColocationSlicesY_Ch0,ColocationSlicesY_Ch1,c= ColocationDelta, s=1)        
        cbar=plt.colorbar(orientation='vertical')
        cbar.set_label('Co-location time, s', rotation=270, labelpad=20)
        
        # linear fit
        slope, intercept, r_value, p_value, std_err = stats.linregress(ColocationDiameterY_Ch0,ColocationDiameterY_Ch1)
        plt.annotate('R = '+str(np.round(r_value,3)), xy=(0.05, 0.95), xycoords='axes fraction')
        
        SizeBins = np.arange(10,1280, step = 10 )
        UpperBound = SizeBins**1.1 + 10
        plt.plot(SizeBins,UpperBound, color ='r', linewidth=2 )
        plt.plot(UpperBound, SizeBins, color ='r', linewidth=2 )
        
        plt.xlabel('Channel 0 Y diameter (bounding box), μm')
        plt.ylabel('Channel 1 Y diameter (bounding box), μm')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim([5,1280])
        plt.ylim([5,1280])
        plt.savefig(PathSave+filena[:-3]+'_colocationCH0vsCH1.png',dpi=200)
        plt.close(fig)
        
#____________________________________________________________________________________
#
# plot array position for stereo particles below certain size

def PlotArrayPositionSmallParticles(PathSave,filena, ColocationMeanXYDiameter_Ch0, ColocationMeanXYDiameter_Ch1, ColocationMIDx_Ch1, ColocationMIDx_Ch0): 
    
    ColocationMIDx_Ch1_lessthan = (np.where(np.logical_and(ColocationMeanXYDiameter_Ch0<20, ColocationMeanXYDiameter_Ch1<20),ColocationMIDx_Ch1,np.nan))
    ColocationMIDx_Ch0_lessthan = (np.where(np.logical_and(ColocationMeanXYDiameter_Ch0<20, ColocationMeanXYDiameter_Ch1<20),ColocationMIDx_Ch0,np.nan)) 
    
    fig=plt.figure(figsize=(7,7)) 
    plt.rcParams.update({'font.size': 12})
    
    HistBins = np.linspace(0,128,num=129, endpoint=True)
    HistBinsMid = (HistBins[:-1] + HistBins[1:]) / 2
    ColocationMIDx_Ch1_lessthan_hist, tmp= np.histogram(ColocationMIDx_Ch1_lessthan,bins=HistBins)
    ColocationMIDx_Ch0_lessthan_hist, tmp= np.histogram(ColocationMIDx_Ch0_lessthan,bins=HistBins)
    plt.plot(HistBinsMid, ColocationMIDx_Ch1_lessthan_hist, label ='Channel 1')
    plt.plot(HistBinsMid, ColocationMIDx_Ch0_lessthan_hist, label ='Channel 0')       
    plt.xlabel('Pixel number')
    plt.ylabel('Counts')
    plt.legend()
    plt.savefig(PathSave+filena[:-3]+'_ArrayPosition.png',dpi=200)
    plt.close(fig)
    
#____________________________________________________________________________________
#
