# -*- coding: utf-8 -*-

#Save stereo images as individual .pngs
#FileName are '2DS raw file name + index number in Colocate_base.... + channel number

#v1 20/04/2021
#original



import numpy as np 
import h5py
import matplotlib.pyplot as plt  
import matplotlib.dates as mdates
import matplotlib.colors as colors
#import matplotlib.colors.Colormap
import datetime
import os
from FlightInfo2DS import GetFlightInfo2DS
import imageio

#FlightNumberStr = 'MAC_233'
Info2DS = GetFlightInfo2DS()
SizeThreshold  = 50 # min size particle to save 

# Flights = ['MAC_235', 'MAC_233', 'MAC_232','MAC_231', 'MAC_230','MAC_228', 'MAC_227',
#            'MAC_226', 'MAC_225','MAC_224','MAC_223','MAC_222','MAC_221',
#            'MAC_220','MAC_219','MAC_218',
#            'C174_dataPC', 'C172_dataPC', 'C171_dataPC', 'C170_dataPC', 'C169_dataPC',
#            'C098_dataPC', 'C097_dataPC'] 

# Flights = ['MAC_231',
#            'C174_dataPC', 'C172_dataPC', 'C171_dataPC', 'C170_dataPC', 'C169_dataPC',
#            'C098_dataPC', 'C097_dataPC']

Flights = ['B895_dataPC']
 
#missed 231

for FlightNumberStr in Flights : 
    print(FlightNumberStr)
    Path2DS = Info2DS[FlightNumberStr,'Path2DS']
    Path2DSsave = Info2DS[FlightNumberStr,'Path2DSsave']
    ThresholdDeltaDimaterY = Info2DS[FlightNumberStr,'ThresholdDeltaDiameterY']
    tmp = [F for F in os.listdir(Path2DSsave) if F.endswith(".h5") and F.startswith('dNdD_L_Colocate_')]
    files = [x.replace('dNdD_L_Colocate_', '') for x in tmp] 
        
    for filena in files :
        print(filena)
        #Load colocation pbp statistics
        SaveFullPath = Path2DSsave + filena[:-3]+'_stereoimages_i/'
        if not os.path.exists(SaveFullPath):
            os.makedirs(SaveFullPath)
        Data_h5 = h5py.File(Path2DSsave + 'Colocate_'+filena, 'r')              
        ColocationParticleBufferTimeS_Ch0=np.array(Data_h5['ColocationParticleBufferTimeS_Ch0'])
        ColocationParticleBufferTimeS_Ch1=np.array(Data_h5['ColocationParticleBufferTimeS_Ch1'])
        ColocationImageID_Ch0=np.array(Data_h5['ColocationImageID_Ch0'])
        ColocationImageID_Ch1=np.array(Data_h5['ColocationImageID_Ch1'])
        ColocationEdgeCH0=np.array(Data_h5['ColocationEdgeCh0'])
        ColocationEdgeCH1=np.array(Data_h5['ColocationEdgeCh1'])
        ColocationMeanXYDiameter_Ch1 =np.array(Data_h5['ColocationMeanXYDiameter_Ch1'])
        ColocationMeanXYDiameter_Ch0 =np.array(Data_h5['ColocationMeanXYDiameter_Ch0'])
        ColocationSlicesY_Ch0 = np.array(Data_h5['ColocationSlicesY_Ch0'])
        ColocationSlicesY_Ch1 = np.array(Data_h5['ColocationSlicesY_Ch1'])
        ColocationSecondsCh0 = np.array(Data_h5['ColocationSecondsCh0'])
        ColocationSecondsCh1 = np.array(Data_h5['ColocationSecondsCh1'])
        Data_h5.close()
        
        #Load images
        Data_h5 = h5py.File(Path2DS+ 'Export_'+filena, 'r')
        ImageTimes=np.array(Data_h5['ImageTimes'][:,0])
        ImageSlices  =np.array(Data_h5['ImageTimes'][:,1])
        ImageID_Ch0 =np.array(Data_h5['ImageTimes'][:,2])
        ImageID_Ch1 =np.array(Data_h5['ImageTimes'][:,2])
        ImageSlices[ImageSlices<0] = np.nan
        #Find start position of image within ImageData
        ImagePosition = np.cumsum(ImageSlices, axis = 0)
        ImagePosition = np.append(0, ImagePosition)
            
        #Indexes to save images
        Stereo_Idxs = np.nonzero(((ColocationMeanXYDiameter_Ch1 > SizeThreshold) | (ColocationMeanXYDiameter_Ch0 > SizeThreshold))
                       & ((ColocationEdgeCH0 == 0) & (ColocationEdgeCH1 == 0))
                       & ((ThresholdDeltaDimaterY == -1) | (np.absolute(ColocationSlicesY_Ch0 - ColocationSlicesY_Ch1) < ThresholdDeltaDimaterY)))
        Stereo_Idxs=Stereo_Idxs[0]
        
        # Number of slices per stereo images
        OutputSlicesY_Ch0 = ColocationSlicesY_Ch0[Stereo_Idxs] /10 # needs to be in pixels not size
        OutputSlicesY_Ch1 = ColocationSlicesY_Ch1[Stereo_Idxs] /10 
        #Position within output image array
        OutputImagePositionCh0 = np.cumsum(OutputSlicesY_Ch0, axis = 0)
        OutputImagePositionCh0 = np.append(0, OutputImagePositionCh0)
        OutputImagePositionCh1 = np.cumsum(OutputSlicesY_Ch1, axis = 0)
        OutputImagePositionCh1 = np.append(0, OutputImagePositionCh1)
        
        # Set up output array 
        OutputImageCh0 = np.ones([128,int(np.nansum(OutputSlicesY_Ch0))], dtype=np.uint8)*255
        OutputImageCh1 = np.ones([128,int(np.nansum(OutputSlicesY_Ch1))], dtype=np.uint8)*255
        
        # select particles and put images in OutputImage
        for j, Idx in enumerate(Stereo_Idxs) : 
            # find each image in array 
            #channel 0
            Ch0i = np.nonzero((ImageTimes == ColocationParticleBufferTimeS_Ch0[Idx]) & (ImageID_Ch0 == ColocationImageID_Ch0[Idx]))
            i = Ch0i[0]
            if (len(i)==0): 
                print('Missing =' + str(i)) # if can't find particle 
            else : 
                if (len(i) > 1): 
                    print('Multiple particles with same ID=' + str(i)) #repeat particle
                    i=i[0]
                if (ImagePosition[i+1]-ImagePosition[i] != (ColocationSlicesY_Ch0[Idx]/10) ): 
                    print('zero slice =' + str(i)) # 0 slice particle   
                else:
                    ImageCH0 = np.array(Data_h5['ImageData'][:,int(ImagePosition[i]):int(ImagePosition[i+1])]) 
                    imageio.imwrite(SaveFullPath+ filena.replace('.h5', '')+'_'+str(Idx)+'_Ch0.png', ImageCH0)
            
            #Channel 1
            Ch1i = np.nonzero((ImageTimes == ColocationParticleBufferTimeS_Ch1[Idx]) & (ImageID_Ch1 == ColocationImageID_Ch1[Idx]))
            i = Ch1i[0]
            if (len(i)==0): 
                print('Missing =' + str(i)) # if can't find particle 
            else : 
                if (len(i) > 1): 
                    print('Multiple particles with same ID=' + str(i)) #repeat particle
                    i=i[0]
                if (ImagePosition[i+1]-ImagePosition[i] != (ColocationSlicesY_Ch1[Idx]/10) ): 
                    print('zero slice =' + str(i)) # 0 slice particle   
                else:
                    ImageCH1 = np.array(Data_h5['ImageData'][:,int(ImagePosition[i]):int(ImagePosition[i+1])]) 
                    imageio.imwrite(SaveFullPath+ filena.replace('.h5', '')+'_'+str(Idx)+'_Ch1.png', ImageCH1)
        
        Data_h5.close()