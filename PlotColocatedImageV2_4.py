# -*- coding: utf-8 -*-
 
# Procedures to plot stereo pairs of images 

# 1) Extract stereo particles using Process2DS_vXXX.py
# 2) PlotAllImages(Info2DS,FlightNumberStr), plot all stereo pairs
#


#v1.0  21/05/2020
#original 

#v2.0 21/08/2020
# Stack panels of 400 slices

#v2.1 28/08/2020
# Colour images grey if size of image above threshold 

#v2.2 06/01/2020
# Use D**x+10 function to colour image grey
# Use number of slices in bbox to flag non stereo pairs

#v2.3 12/2/2021
#Get paths from Info2DS
#Included option not to flag images.
# Added minimum size threshold to plot image

#v2.4 28/3/2021
# Simplified method to search for individual image within 'ImageData'. 
# Only selects files to plot images where colocation .h5 file already exists
# Catches error if a stereo image can't be found within 'ImageData'.


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
# Loop through .h5 image files in folder


def PlotAllImages(Info2DS,FlightNumberStr):  

    Path2DSsave = Info2DS[FlightNumberStr,'Path2DSsave']
    tmp = [F for F in os.listdir(Path2DSsave) if F.endswith(".h5") and F.startswith('Colocate_')]
    files = [x.replace('Colocate_', '') for x in tmp]
    print(files)
    
    
    for filena in files :
        filena.replace('Colocate_','')
        BatchPlotImages_2Channels(Info2DS,FlightNumberStr,filena)
        
    #filena = files[0] # select file index to plot. 
    #BatchPlotImages_2Channels(Info2DS,FlightNumberStr,filena)

    

#__________________________________________________________________________________
# Plot all colocated images in .h5 image file

def BatchPlotImages_2Channels(Info2DS,FlightNumberStr,filena):
    
    Path2DS = Info2DS[FlightNumberStr,'Path2DS']
    Path2DSsave = Info2DS[FlightNumberStr,'Path2DSsave']
    ThresholdDeltaDimaterY = Info2DS[FlightNumberStr,'ThresholdDeltaDiameterY']
    
    print(filena)
    
    SavePath = Path2DSsave + filena[:-3]+'/'
    if not os.path.exists(SavePath):
        os.makedirs(SavePath)
    
    #ThresholdDeltaDimaterY = 40 # threshold size difference between y dimension of stereo images. Images above threshold are coloured grey.
    #ThresholdDiameterYExponent = 1.1

    #Nfraction2Plot = 1 # fraction of particles to plot 
    MinSizeThreshold  = 50 # min size particle to plot
    MaxSizeThreshold = 2000
    Nslices = 1600 # Total number of slices per plot
    Npanels =4 # number of panels per plot
    Data_h5 = h5py.File(Path2DSsave + 'Colocate_'+filena, 'r')              
    ColocationParticleBufferTimeS_Ch0=np.array(Data_h5['ColocationParticleBufferTimeS_Ch0'])
    ColocationParticleBufferTimeS_Ch1=np.array(Data_h5['ColocationParticleBufferTimeS_Ch1'])
    ColocationImageID_Ch0=np.array(Data_h5['ColocationImageID_Ch0'])
    ColocationImageID_Ch1=np.array(Data_h5['ColocationImageID_Ch1'])
    ColocationEdgeCH0=np.array(Data_h5['ColocationEdgeCh0'])
    ColocationEdgeCH1=np.array(Data_h5['ColocationEdgeCh1'])
    ColocationMeanXYDiameter_Ch1 =np.array(Data_h5['ColocationMeanXYDiameter_Ch1'])
    ColocationMeanXYDiameter_Ch0 =np.array(Data_h5['ColocationMeanXYDiameter_Ch0'])
    ColocationDelta = np.array(Data_h5['ColocationDelta'])
    ColocationSlicesY_Ch0 = np.array(Data_h5['ColocationSlicesY_Ch0'])
    ColocationSlicesY_Ch1 = np.array(Data_h5['ColocationSlicesY_Ch1'])

    Data_h5.close()

    
    AllColocationImages=np.zeros([128,Nslices])
    
    ZerosThreesZeros =np.append(np.zeros([128,1]),3*np.ones([128,1]), axis = 1)
    ZerosThreesZeros =np.append(ZerosThreesZeros,np.zeros([128,1]), axis = 1)   
    
    if 1 == 1:
        FileName = 'Export_'+filena
        x=0 #Image idx in colocation file
        while x < len(ColocationImageID_Ch0) :
            
            # new plot
            AllColocationImages[:,:]=0 
            TotalSize = 0
            FirstTime = ColocationParticleBufferTimeS_Ch0[x]
            
            while  TotalSize < Nslices and x < len(ColocationImageID_Ch0):
                #print(x)
                #if (((ColocationMeanXYDiameter_Ch1[x] > SizeThreshold) or (ColocationMeanXYDiameter_Ch0[x] > SizeThreshold) ) and ColocationEdgeCH0[x] == 0 and ColocationEdgeCH1[x] == 0): 
                if (((ColocationMeanXYDiameter_Ch1[x] > MinSizeThreshold) or (ColocationMeanXYDiameter_Ch0[x] > MinSizeThreshold)) 
                    and ((ColocationMeanXYDiameter_Ch1[x] < MaxSizeThreshold) or (ColocationMeanXYDiameter_Ch0[x] > MaxSizeThreshold))
                    and (ColocationEdgeCH0[x] == 0 and ColocationEdgeCH1[x] == 0)): 
                     
                    # Flag images using maxD > minD**1.1 + 10
                    #MinDiameterY = min(ColocationSlicesY_Ch0[x], ColocationSlicesY_Ch1[x])
                    #MaxDiameterY = max(ColocationSlicesY_Ch0[x], ColocationSlicesY_Ch1[x])
                    #PairFlag = np.where(MaxDiameterY>(MinDiameterY**ThresholdDiameterYExponent + 10), 0,1)
                    
                    #Flag images using maxD-minD
                    DeltaDiameterY = np.absolute(ColocationSlicesY_Ch0[x] - ColocationSlicesY_Ch1[x])
                    if ThresholdDeltaDimaterY == -1 :
                        PairFlag = 1 #No filtering for DeltaDimaterY
                    else :
                        PairFlag = np.where(DeltaDiameterY>ThresholdDeltaDimaterY, 0,1)
                    
                    ImagePair =  CombineImage2Channels(Path2DS,FileName,ColocationParticleBufferTimeS_Ch0[x],ColocationImageID_Ch0[x],ColocationParticleBufferTimeS_Ch1[x],ColocationImageID_Ch1[x],PairFlag)                        
                    ImagePair = np.append(ImagePair,ZerosThreesZeros, axis = 1)
                    ImageSize = np.size(ImagePair,axis=1)
                    IDXstart = TotalSize
                    IDXend = TotalSize+ImageSize

                    if TotalSize + ImageSize > Nslices -1 : # check that the new image will fit in the plot
                        TotalSize = Nslices # Image too large to fit
                    else: # add image to plot
                        TotalSize += ImageSize
                        AllColocationImages[:,IDXstart:IDXend] =ImagePair
                        LastTime = ColocationParticleBufferTimeS_Ch0[x]
                        x+=1
                else:
                    LastTime = ColocationParticleBufferTimeS_Ch0[x]
                    x+=1
                
            
            PlotAllColocationImages(AllColocationImages, FirstTime, LastTime,Npanels,Nslices,filena,SavePath)
            #x+=1 
            

#__________________________________________________________________________________

        
def CombineImage2Channels(ImagePath,FileName,ParticleBufferTime_Ch0,ParticleID_Ch0,ParticleBufferTime_Ch1,ParticleID_Ch1, PairFlag):
    
    Data_h5 = h5py.File(ImagePath + FileName, 'r')
    ImageTimes=np.array(Data_h5['ImageTimes'][:,0])
    ImageSlices  =np.array(Data_h5['ImageTimes'][:,1])
    ImageID_Ch0 =np.array(Data_h5['ImageTimes'][:,2])
    ImageID_Ch1 =np.array(Data_h5['ImageTimes'][:,2])
    ImageSlices[ImageSlices<0] = np.nan
    #Find start position of image within ImageData
    ImagePosition = np.cumsum(ImageSlices, axis = 0)
    ImagePosition = np.append(0, ImagePosition)
    
    # Channel 0 
    #Search for particle image
    idx = np.nonzero((ImageTimes == ParticleBufferTime_Ch0) & (ImageID_Ch0 == ParticleID_Ch0))
    i = idx[0]
    
    if (len(i)==0): 
        print('Missing =' + str(i)) # if can't find particle 
        ImageCH0= np.ones([128,1])*255 #return blank image
    else : 
        if len(i) > 1:
            print('Multiple particles with same ID=' + str(i))
            i=i[0]
        ImageCH0 = np.array(Data_h5['ImageData'][:,int(ImagePosition[i]):int(ImagePosition[i+1])]) 
        ImageCH0[ImageCH0 == 0 ] = 1
        ImageCH0[ImageCH0 == 255 ] = 0  

    # Channel 1
    #Search for particle image
    idx = np.nonzero((ImageTimes == ParticleBufferTime_Ch1) & (ImageID_Ch1 == ParticleID_Ch1))
    i = idx[0]
    if (len(i)==0): 
        print('Missing =' + str(i)) # if can't find particle 
        ImageCH1= np.ones([128,1])*255 #return blank image
    else : 
        if len(i) > 1:
            print('Multiple particles with same ID=' + str(i))
            i=i[0]
        ImageCH1 = np.array(Data_h5['ImageData'][:,int(ImagePosition[i]):int(ImagePosition[i+1])]) 
        ImageCH1[ImageCH1 == 0 ] = 2
        ImageCH1[ImageCH1 == 255 ] = 0

    ImageCombine = np.append(ImageCH0, ImageCH1, axis = 1)

    if PairFlag == 0 : 
        ImageCombine = np.where(ImageCombine >0, 4, ImageCombine) # colour differently if PairFlag == 0
    
    Data_h5.close()  
    
    return ImageCombine#, ImageSize

#__________________________________________________________________________________
   
# plot continuous stream of colocated images

def PlotAllColocationImages(AllColocationImages, FirstTime, LastTime,Npanels,Nslices,filena,SavePath):

    #PixelSize= 10
    ArrayWidth = 128
    Nslices /=Npanels

    fig=plt.figure(figsize=(10, 10*((ArrayWidth*Npanels)/Nslices)))


    plt.suptitle(str(datetime.timedelta(seconds=FirstTime))+' to '+str(datetime.timedelta(seconds=LastTime)))

    
    for i in range(1,Npanels+1,1) :
        plt.subplot(Npanels,1,i)
        cmap = colors.ListedColormap(["w", "darkkhaki", "royalblue", "k",'silver'])
        plt.pcolormesh(AllColocationImages,  vmin=0, vmax=4, cmap=cmap)
        plt.axis('off')
        plt.axhline(y=128, color='k', linestyle=':')
        plt.axhline(y=-1, color='k', linestyle=':')
        plt.xlim([Nslices*(i-1),Nslices*i])
        
        
    plt.savefig(SavePath+'ColocateImages_'+filena[:-3]+'_'+str(FirstTime)+'.png',dpi=200)
    
    plt.close(fig)
   
 
#__________________________________________________________________________________

