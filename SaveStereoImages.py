# -*- coding: utf-8 -*-


# Save stereo images .h5 file
#v1 17/3/21
#Original

import numpy as np 
import h5py
import matplotlib.pyplot as plt  
import matplotlib.dates as mdates
import matplotlib.colors as colors
#import matplotlib.colors.Colormap
import datetime
import os
from Process2DS_v2_5 import GetFlightInfo2DS



    
#__________________________________________________________________________________
#
    #ImageCH0[ImageCH0 == 0 ] = 1
    #ImageCH0[ImageCH0 == 255 ] = 0



def SaveStereoImagesh5(Info2DS,FlightNumberStr,SizeThreshold):
    
    # Path2DSsave = 'C:/Users/mbexjso2/OneDrive - The University of Manchester/Documents/Work/Picasso/C174/'
    # Path2DS = 'C:/Users/mbexjso2/OneDrive - The University of Manchester/Documents/Work/Picasso/C174/'
    # ThresholdDeltaDimaterY = 40
    # filena = 'base190523063310.h5'
    
    Path2DS = Info2DS[FlightNumberStr,'Path2DS']
    Path2DSsave = Info2DS[FlightNumberStr,'Path2DSsave']
    ThresholdDeltaDimaterY = Info2DS[FlightNumberStr,'ThresholdDeltaDiameterY']
    tmp = [F for F in os.listdir(Path2DSsave) if F.endswith(".h5") and F.startswith('dNdD_L_Colocate_')]
    files = [x.replace('dNdD_L_Colocate_', '') for x in tmp]
    filena = files[0] # select file index to plot. 
    
    
    #Load colocation pbp statistics
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
                   & ((ColocationEdgeCH0 == 0) | (ColocationEdgeCH1 == 0))
                   & (np.absolute(ColocationSlicesY_Ch0 - ColocationSlicesY_Ch1) < ThresholdDeltaDimaterY))
    Stereo_Idxs=Stereo_Idxs[0]
    
    # Number of slice of stereo images
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
    
    # Look through select particles and put images in OutputImage
    for j, Idx in enumerate(Stereo_Idxs) : 
        # find each image in array 
    
        Ch0i = np.nonzero((ImageTimes == ColocationParticleBufferTimeS_Ch0[Idx]) & (ImageID_Ch0 == ColocationImageID_Ch0[Idx]))
        i = Ch0i[0]
        if len(i) != 1: print('Cant find particle =' + str(i)) # if can't find particle  
            #i=i[0]
        
        ImageCH0 = np.array(Data_h5['ImageData'][:,int(ImagePosition[i]):int(ImagePosition[i+1])]) 
        #Add to output array
        OutputImageCh0[:,int(OutputImagePositionCh0[j]):int(OutputImagePositionCh0[j+1])] = ImageCH0
        
        #Channel 1
        Ch1i = np.nonzero((ImageTimes == ColocationParticleBufferTimeS_Ch1[Idx]) & (ImageID_Ch1 == ColocationImageID_Ch1[Idx]))
        i = Ch1i[0]
        if len(i) != 1: print('Cant find particle =' + str(i)) # if can't find particle 
            #i=i[0]
        ImageCH1 = np.array(Data_h5['ImageData'][:,int(ImagePosition[i]):int(ImagePosition[i+1])]) 
        #Add to output array
        OutputImageCh1[:,int(OutputImagePositionCh1[j]):int(OutputImagePositionCh1[j+1])] = ImageCH1
       
    
    Data_h5.close()
    
    # Save the images as .h5
    h5f = h5py.File(Path2DSsave+'StereoImages_'+filena, 'w')
    h5f.create_dataset('ImageCh0', data=OutputImageCh0)
    h5f.create_dataset('ImagePositionCh0', data=OutputImagePositionCh0)
    h5f.create_dataset('SecondsCh0', data=ColocationSecondsCh0)
    h5f.create_dataset('SlicesY_Ch0', data=OutputSlicesY_Ch0)
    h5f.create_dataset('ImageCh1', data=OutputImageCh1)
    h5f.create_dataset('ImagePositionCh1', data=OutputImagePositionCh1)
    h5f.create_dataset('SecondsCh1', data=ColocationSecondsCh1)
    h5f.create_dataset('SlicesY_Ch1', data=OutputSlicesY_Ch1)
    h5f.close()
    
    # plt.figure()
    # plt.pcolormesh(OutputImageCh0)
    # plt.xlim([1000,1500])
    
    # plt.figure()
    # plt.pcolormesh(OutputImageCh1)
    # plt.xlim([1000,500])


Info2DS = GetFlightInfo2DS()
FlightNumberStr = 'C174_dataPC'
SizeThreshold  = 50 # min size particle to save 
#SaveStereoImagesh5(Info2DS,FlightNumberStr,SizeThreshold)

#Combines files
# files = [F for F in os.listdir(Path2DSsave) if F.endswith(".h5") and F.startswith('StereoImages_')]
# for F in files : 
#     Data_h5 = h5py.File(Path2DSsave + F, 'r')              
#     ImageCh0=np.array(Data_h5['ImageCh0'])
#     h5f.create_dataset('ImagePositionCh0', data=ImagePositionCh0)
#     h5f.create_dataset('SecondsCh0', data=ColocationSecondsCh0)
#     h5f.create_dataset('SecondsCh0', data=SecondsCh0)

#     h5f.close()




