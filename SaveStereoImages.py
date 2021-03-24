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

# for each base.h5 file create 1 .h5 with stereo images

def SaveStereoImagesh5(Info2DS,FlightNumberStr,SizeThreshold):
    
    Path2DS = Info2DS[FlightNumberStr,'Path2DS']
    Path2DSsave = Info2DS[FlightNumberStr,'Path2DSsave']
    ThresholdDeltaDimaterY = Info2DS[FlightNumberStr,'ThresholdDeltaDiameterY']
    tmp = [F for F in os.listdir(Path2DSsave) if F.endswith(".h5") and F.startswith('dNdD_L_Colocate_')]
    files = [x.replace('dNdD_L_Colocate_', '') for x in tmp] 
    
    for filena in files : 
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
        
        # Look through select particles and put images in OutputImage
        for j, Idx in enumerate(Stereo_Idxs) : 
            # find each image in array 
            Ch0i = np.nonzero((ImageTimes == ColocationParticleBufferTimeS_Ch0[Idx]) & (ImageID_Ch0 == ColocationImageID_Ch0[Idx]))
            i = Ch0i[0]
            if len(i) != 1: print('Cant find particle =' + str(i)) # if can't find particle  
            else:
                ImageCH0 = np.array(Data_h5['ImageData'][:,int(ImagePosition[i]):int(ImagePosition[i+1])]) 
                #Add to output array
                OutputImageCh0[:,int(OutputImagePositionCh0[j]):int(OutputImagePositionCh0[j+1])] = ImageCH0
            
            #Channel 1
            Ch1i = np.nonzero((ImageTimes == ColocationParticleBufferTimeS_Ch1[Idx]) & (ImageID_Ch1 == ColocationImageID_Ch1[Idx]))
            i = Ch1i[0]
            if len(i) != 1: print('Cant find particle =' + str(i)) # if can't find particle 
            else :
                ImageCH1 = np.array(Data_h5['ImageData'][:,int(ImagePosition[i]):int(ImagePosition[i+1])]) 
                #Add to output array
                OutputImageCh1[:,int(OutputImagePositionCh1[j]):int(OutputImagePositionCh1[j+1])] = ImageCH1
           
        Data_h5.close()
        
        # Save the images as .h5
        h5f = h5py.File(Path2DSsave+'StereoImages_'+filena, 'w')
        h5f.create_dataset('ImageCh0', data=OutputImageCh0)
        h5f.create_dataset('ImagePositionCh0', data=OutputImagePositionCh0)
        h5f.create_dataset('SecondsCh0', data=ColocationSecondsCh0[Stereo_Idxs])
        h5f.create_dataset('SlicesY_Ch0', data=OutputSlicesY_Ch0)
        h5f.create_dataset('ImageCh1', data=OutputImageCh1)
        h5f.create_dataset('ImagePositionCh1', data=OutputImagePositionCh1)
        h5f.create_dataset('SecondsCh1', data=ColocationSecondsCh1[Stereo_Idxs])
        h5f.create_dataset('SlicesY_Ch1', data=OutputSlicesY_Ch1)
        h5f.close()
        
#__________________________________________________________________________________
#
# create one stereo image file per flight

Info2DS = GetFlightInfo2DS()
FlightNumberStr = 'C172_dataPC'
SizeThreshold  = 0 # min size particle to save 
SaveStereoImagesh5(Info2DS,FlightNumberStr,SizeThreshold)

vnumber = 0
rnumber = 0

Path2DSsave = Info2DS[FlightNumberStr,'Path2DSsave']
FlightNumber = Info2DS[FlightNumberStr,'FlightNumber'] 
FlightDate = (Info2DS[FlightNumberStr,'FlightDate']).astype(datetime.datetime)


#Combine files
files = [F for F in os.listdir(Path2DSsave) if F.endswith(".h5") and F.startswith('StereoImages_')]
for i,F in enumerate(files) : 
    
    Data_h5 = h5py.File(Path2DSsave + F, 'r')              
    tmpImageCh0=np.array(Data_h5['ImageCh0'])
    tmpImagePositionCh0=np.array(Data_h5['ImagePositionCh0'])
    tmpSecondsCh0=np.array(Data_h5['SecondsCh0'])
    tmpSlicesY_Ch0=np.array(Data_h5['SlicesY_Ch0'])
    tmpImageCh1=np.array(Data_h5['ImageCh1'])
    tmpImagePositionCh1=np.array(Data_h5['ImagePositionCh1'])
    tmpSecondsCh1=np.array(Data_h5['SecondsCh1'])
    tmpSlicesY_Ch1=np.array(Data_h5['SlicesY_Ch1'])
    Data_h5.close()
    
    print(len(tmpSecondsCh1))
    print(len(tmpImagePositionCh1))
    
    if i == 0 : 
        ImageCh0 = tmpImageCh0
        ImagePositionCh0= tmpImagePositionCh0
        SecondsCh0=tmpSecondsCh0
        SlicesY_Ch0=tmpSlicesY_Ch0
        ImageCh1=tmpImageCh1
        ImagePositionCh1=tmpImagePositionCh1
        SecondsCh1=tmpSecondsCh1
        SlicesY_Ch1=tmpSlicesY_Ch1
    else:  
        ImageCh0= np.append(ImageCh0,tmpImageCh0, axis=1)
        ImagePositionCh0= np.append(ImagePositionCh0,tmpImagePositionCh0[1:]+ImagePositionCh0[-1], axis=0)
        SecondsCh0= np.append(SecondsCh0, tmpSecondsCh0, axis=0)
        SlicesY_Ch0= np.append(SlicesY_Ch0,tmpSlicesY_Ch0, axis=0)
        ImageCh1= np.append(ImageCh1,tmpImageCh1, axis=1)
        ImagePositionCh1= np.append(ImagePositionCh1,tmpImagePositionCh1[1:]+ImagePositionCh1[-1], axis=0)
        SecondsCh1= np.append(SecondsCh1, tmpSecondsCh1, axis=0)
        SlicesY_Ch1= np.append(SlicesY_Ch1, tmpSlicesY_Ch1, axis=0)
   

# plt.pcolormesh(ImageCh1)
# plt.vlines(ImagePositionCh1, ymin=0,ymax=128)
# plt.xlim([2000,2500])

Mergedfilename='uman-2ds_faam_'+FlightDate.strftime("%Y%m%d")+'_v'+str(vnumber)+'_r'+str(rnumber)+'_'+FlightNumber+'_stereo_images.h5'
#save to file
h5f = h5py.File(Path2DSsave+Mergedfilename, 'w')
h5f.create_dataset('ImageCh0', data=ImageCh0)
h5f.create_dataset('ImagePositionCh0', data=ImagePositionCh0)
h5f.create_dataset('SecondsCh0', data=SecondsCh0)
h5f.create_dataset('ImageWidth_Ch0', data=SlicesY_Ch0)
h5f.create_dataset('ImageCh1', data=ImageCh1)
h5f.create_dataset('ImagePositionCh1', data=ImagePositionCh1)
h5f.create_dataset('SecondsCh1', data=SecondsCh1)
h5f.create_dataset('ImageWidthCh1', data=SlicesY_Ch1)
h5f.close()

