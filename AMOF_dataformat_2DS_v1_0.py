# -*- coding: utf-8 -*-

import numpy as np
import h5py
import os
import datetime
import netCDF4
import matplotlib.pyplot as plt 
import pandas as pd
from Process2DS_v2_5 import BatchBothChannels,HybridStereoProcessing
from FlightInfo2DS import GetFlightInfo2DS

# Generate an .nc file for the 2DS that meets CEDA criteria

# v1.0 21/01/2021
# Original
# Generate file format for upload to CEDA.
# Correction to PSD for using constant TAS in data processing
# Change "_" to "-" as requested by CEDA


# #_________________________________________________________________
#load .h5 that has been generated using OASIS containing PSDs and data flags

def LoadOAPH5(DataPath, DataFileName):
    
    Data_h5 = h5py.File(DataPath+DataFileName, 'r')              
    DataFlag = np.array(Data_h5['DataFlag']).astype(np.uint32)
    Psd= np.array(Data_h5['PSD_All_accept_submit'])
    Counts= np.array(Data_h5['PSD_Counts_All_accept_submit'])
    Size_mid= np.array(Data_h5['Size_mid_submit'])
    Size_edge= np.array(Data_h5['Size_submit'])
    Time_edge= np.array(Data_h5['Time_edge_submit'])
    Time_mid= np.array(Data_h5['Time_mid_submit'])
    Data_h5.close()
    
    
    return DataFlag, Psd, Counts, Size_mid, Size_edge, Time_mid, Time_edge 

#_________________________________________________________________
# Change timebase of FAAM core variable to OAP timebase
# Also change data fill values to FillValueNew

def CovertTimeBase2NC(FAAMCore, VariableStr, TimeIndex, TimeDiff, FillValueOld,FillValueNew):
    tmp = np.array(FAAMCore[VariableStr][:])
    tmp_nc = tmp[TimeIndex]
    tmp_nc[TimeDiff>=1] = FillValueNew
    tmp_nc[tmp_nc == FillValueOld]= FillValueNew
    return tmp_nc

#_________________________________________________________________

# Load faam core file and change time base and format ready for AMOF data submisison

def LoadCoreVariables4AMOF(CorePath,CoreFileName,Time_mid):

    #Time_mid is the time base of returned array
    
    FAAMCore = netCDF4.Dataset(CorePath+CoreFileName)
    
    TimeCore=np.array(FAAMCore['Time'][:])
    FlightDate= datetime.datetime(int(CoreFileName[10:14]), int(CoreFileName[14:16]), int(CoreFileName[16:18]), 0, 0, 0)
    # index to change core timebase to Oap timebase
    TimeIndex = np.searchsorted(TimeCore,Time_mid)
    TimeDiff = np.absolute(Time_mid-TimeCore[TimeIndex])

    FillValueOld = -9999
    FillValueNew = -1E20
    
    LAT_GIN_nc = CovertTimeBase2NC(FAAMCore, 'LAT_GIN', TimeIndex, TimeDiff,FillValueOld,FillValueNew)
    LON_GIN_nc = CovertTimeBase2NC(FAAMCore, 'LON_GIN', TimeIndex, TimeDiff,FillValueOld,FillValueNew)
    ALT_GIN_nc = CovertTimeBase2NC(FAAMCore, 'ALT_GIN', TimeIndex, TimeDiff,FillValueOld,FillValueNew)
    TAS_RVSM_nc = CovertTimeBase2NC(FAAMCore, 'TAS_RVSM', TimeIndex, TimeDiff,FillValueOld,FillValueNew)
    TRCK_GIN_nc = CovertTimeBase2NC(FAAMCore, 'TRCK_GIN', TimeIndex, TimeDiff,FillValueOld,FillValueNew)
    GSPD_GIN_nc = CovertTimeBase2NC(FAAMCore, 'GSPD_GIN', TimeIndex, TimeDiff,FillValueOld,FillValueNew)
    PTCH_GIN_nc = CovertTimeBase2NC(FAAMCore, 'PTCH_GIN', TimeIndex, TimeDiff,FillValueOld,FillValueNew)
    PITR_GIN_nc = CovertTimeBase2NC(FAAMCore, 'PITR_GIN', TimeIndex, TimeDiff,FillValueOld,FillValueNew)
    ROLL_GIN_nc = CovertTimeBase2NC(FAAMCore, 'ROLL_GIN', TimeIndex, TimeDiff,FillValueOld,FillValueNew)
    ROLR_GIN_nc = CovertTimeBase2NC(FAAMCore, 'ROLR_GIN', TimeIndex, TimeDiff,FillValueOld,FillValueNew)
    HDG_GIN_nc = CovertTimeBase2NC(FAAMCore, 'HDG_GIN', TimeIndex, TimeDiff,FillValueOld,FillValueNew)
    HDGR_GIN_nc = CovertTimeBase2NC(FAAMCore, 'HDGR_GIN', TimeIndex, TimeDiff,FillValueOld,FillValueNew)
    
    FAAMCore.close()
    
    return LAT_GIN_nc, LON_GIN_nc, ALT_GIN_nc, TAS_RVSM_nc, TRCK_GIN_nc, GSPD_GIN_nc, PTCH_GIN_nc, PITR_GIN_nc, ROLL_GIN_nc, ROLR_GIN_nc, HDG_GIN_nc, HDGR_GIN_nc

#_________________________________________________________________


FlightNumberStr = 'C168'

FillValue = -1E20 # Fill value of output data file
Info2DS = GetFlightInfo2DS() # list of paths and variables
BatchBothChannels(Info2DS,FlightNumberStr) # locate stereo particles and save psd files for each 2DS file 
Psd, Counts, DataFlag, Time_mid, Size_mid = HybridStereoProcessing(Info2DS,FlightNumberStr, FillValue) # Merge files, create hybrid (stereo, traditional) and create flag
Psd = np.float32(Psd)
Counts = np.float32(Counts)
Size_edge = np.append([5],Size_mid+5, axis=0)

# Get paths and data processing settings for file.
FlightNumber = Info2DS[FlightNumberStr,'FlightNumber'] 
FlightDate = (Info2DS[FlightNumberStr,'FlightDate']).astype(datetime.datetime)
DataPath= Info2DS[FlightNumberStr,'Path2DS']
IAT_threshold = Info2DS[FlightNumberStr, 'IAT_threshold']
ColocationThreshold = Info2DS[FlightNumberStr, 'ColocationThreshold']
ThresholdDeltaDiameterY = Info2DS[FlightNumberStr, 'ThresholdDeltaDiameterY']
ThresholdSize =  Info2DS[FlightNumberStr,'ThresholdSize']

CorePath = Info2DS[FlightNumberStr,'CorePath']
CoreFileName = Info2DS[FlightNumberStr,'CoreFileName']
rnumber = 0 # revision number
vnumber = '001'# Don't change relates to versions of processing software (currently not defined)
NormFlag = 1 # Change if data normalised, 0 = dN, 1 = dN/dD
UnitConversion =1E-3 # Factor to convert input PSD to cm-3 
TAS_correctionFlag = 1 # 1 = apply correction, 0 = dont. If fixed TAS was used in Svol calc


#Load PSDs file
#DataFlag, Psd, Counts, Size_mid, Size_edge, Time_mid, Time_edge = LoadOAPH5(DataPath, DataFileName)
#DataFlag += 1 # Conversion factor to flag data so that '0 = not used, 1 = valid data, 2 = reduced quality data. 3 = missing data'
DataFlag_2D = np.array([DataFlag,]*len(Size_mid)).transpose() # convert flag to a 2d array

# Load variables from core files
LAT_GIN_nc, LON_GIN_nc, ALT_GIN_nc, TAS_RVSM_nc, TRCK_GIN_nc, GSPD_GIN_nc, PTCH_GIN_nc, PITR_GIN_nc, ROLL_GIN_nc, ROLR_GIN_nc, HDG_GIN_nc, HDGR_GIN_nc = LoadCoreVariables4AMOF(CorePath,CoreFileName,Time_mid)

#Convert time to seconds since 1970
DateTime_df = pd.to_datetime(FlightDate) + pd.to_timedelta(Time_mid, unit='s')

TimeSecondsSince = np.array((DateTime_df-datetime.datetime(1970,1,1)) / np.timedelta64(1,'s'))
year = np.array((DateTime_df.year)).astype(np.int32) 
month = np.array((DateTime_df.month)).astype(np.int32) 
day = np.array((DateTime_df.day)).astype(np.int32) 
hour = np.array((DateTime_df.hour)).astype(np.int32) 
minute = np.array((DateTime_df.minute)).astype(np.int32) 
second = np.array((DateTime_df.second))
day_of_year = np.array((DateTime_df.dayofyear)) + Time_mid/86400 # day of year with fraction of day

#Convert to cm-3
Psd *= UnitConversion

# correct psd if constant TAS was used in data processing
if TAS_correctionFlag == 1:    
    TAS_constant =  Info2DS[FlightNumberStr,'TAS']
    TAS_correction = (np.where(TAS_RVSM_nc == FillValue, 1, TAS_constant/TAS_RVSM_nc))
    Psd *= TAS_correction[:,None]

#Use chosen FillValue                    
Psd[Psd < 0] = FillValue
Counts[Counts < 0] = FillValue 
Psd[np.isnan(Psd)] = FillValue
Counts[np.isnan(Counts)] = FillValue 



#######################
#Save data as netcdf
Nc_filename='uman-2ds_faam_'+FlightDate.strftime("%Y%m%d")+'_v'+vnumber+'_r'+str(rnumber)+'_'+FlightNumber+'_particle-size-distribution.nc'
print(Nc_filename)
NcFile= netCDF4.Dataset(DataPath+Nc_filename,'w', format='NETCDF4') #'w' stands for write

# Global attributes
NcFile.conventions ='CF-1.6'
NcFile.source = 'University of Manchester 2D-S'
NcFile.instrument_manufacturer = 'Spec Inc., USA'
NcFile.instrument_model = '2DS'
NcFile.instrument_serial_number = 'Not available'
NcFile.instrument_software = 'spec2d.exe'
NcFile.instrument_software_version ='Not available'
NcFile.creator_name = 'Sebastian OShea'
NcFile.creator_email = 'sebastian.oshea@manchester.ac.uk'
NcFile.url = 'https://orcid.org/0000-0002-0489-1723'
NcFile.institution = 'University of Manchester'
NcFile.processing_software_url = 'Not available'
NcFile.processing_software_version = 'OASIS_routines_v1.49962, Process2DS_v2_4'
NcFile.calibration_sensitivity = 'The instrument counting uncertainty for each time and size bin is given by (100% / N^(1/2)), where N is the number of counts in the bin.'
NcFile.calibration_certification_date ='not applicable'
NcFile.calibration_certification_url ='not applicable'
NcFile.sampling_interval='1 second'
NcFile.averaging_interval='1 second'
NcFile.product_version='v'+str(rnumber+1)+'.0'
NcFile.processing_level= int(1)
NcFile.last_revised_date= datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
#NcFile.project='Parameterizing Ice Clouds using Airborne obServationS and triple-frequency dOppler radar (PICASSO)'
#NcFile.project_principal_investigator='Jonathan Crosier'
#NcFile.project_principal_investigator_email='jonathan.crosier@manchester.ac.uk'
#NcFile.project_principal_investigator_url='https://orcid.org/0000-0002-3086-4729'
NcFile.project='Met Office PIKNMIX campaigns: cloud pyhsics and radiation events (PIKNMIX-F) project'
#NcFile.project_principal_investigator='Jonathan Crosier'
#NcFile.project_principal_investigator_email='jonathan.crosier@manchester.ac.uk'
#NcFile.project_principal_investigator_url='https://orcid.org/0000-0002-3086-4729'
NcFile.licence='Data usage licence - UK Government Open Licence agreement: http://www.nationalarchives.gov.uk/doc/open-government-licence'
NcFile.acknowledgement='To discuss use and interest in this data please contact Sebastian OShea.'
NcFile.platform='faam'
NcFile.platform_type='mobile'
NcFile.deployment_mode='air'
NcFile.title='In situ measurement of particle size distributions from the FAAM BAe-146 research aircraft.' 
NcFile.featureType='timeSeries'
NcFile.time_coverage_start=FlightDate.strftime("%Y-%m-%dT%H:%M:%S")
NcFile.time_coverage_end=(FlightDate+datetime.timedelta(seconds=int(Time_mid[-1]))).strftime("%Y-%m-%dT%H:%M:%S")
NcFile.geospatial_bounds= str(np.around(np.min(LAT_GIN_nc[LAT_GIN_nc != FillValue]),2))+'N '+str(np.around(np.min(LON_GIN_nc[LON_GIN_nc != FillValue]),2))+'E, '+str(np.around(np.max(LAT_GIN_nc[LAT_GIN_nc != FillValue]),2))+'N '+str(np.around(np.max(LON_GIN_nc[LON_GIN_nc != FillValue]),2))+'E'
NcFile.platform_altitude='Not applicable'
NcFile.location_keywords = 'Stornoway, UK'
#NcFile.amf_vocabularies_release = 'https://github.com/ncasuk/AMF_CVs/releases/tag/v2.0.0'
NcFile.history = FlightDate.strftime("%Y-%m-%dT%H:%M:%SZ")+' - Data collected.'
NcFile.comment = 'Data processing was performed to derive size segregated number concentrations as in OShea et al. (https://doi.org/10.5194/amt-2020-265). Particle size has been calculated as the mean particle extent parallel and perpendicular to the probe optical array. Particles with inter-arrival times <'+str(IAT_threshold)+'s have been rejected. The PSD <'+str(ThresholdSize)+' um has been calculated using stereo pairs. Pairs of images were identified as stereo if the time difference between them is <'+str(ColocationThreshold)+'s and diameter y difference is <='+str(ThresholdDeltaDiameterY)+'. Particles in contact with the edge of the sample array have been rejected. Other aircraft variables have been extracted from the file '+CoreFileName+'. The instrument counting uncertainty for each time and size bin is given by (100% / N^(1/2)), where N is the number of counts in the bin.'
NcFile.measurement_technique = 'optical'

#Dimensions
time = NcFile.createDimension('time', len(TimeSecondsSince))
latitude = NcFile.createDimension('latitude', len(TimeSecondsSince))
longitude = NcFile.createDimension('longitude', len(TimeSecondsSince))
index = NcFile.createDimension('index', len(Size_mid))

# Probe variables
time = NcFile.createVariable(varname='time', dimensions=('time',),datatype='float64')
time[:] = TimeSecondsSince
time.units = 'seconds since 1970-01-01 00:00:00'
time.long_name = 'Time (seconds since 1970-01-01)'
time.calendar = 'standard'
time.axis = 'T'
time.valid_min= np.min(TimeSecondsSince)
time.valid_max= np.max(TimeSecondsSince)

ambient_particle_diameter= NcFile.createVariable(varname='ambient_particle_diameter', dimensions=('index',),datatype='float32')
ambient_particle_diameter[:] = Size_mid
ambient_particle_diameter.long_name = 'Ambient Particle Diameter'
ambient_particle_diameter.units = 'um'
ambient_particle_diameter.valid_min= np.min(Size_mid[Size_mid!=FillValue])
ambient_particle_diameter.valid_max= np.max(Size_mid[Size_mid!=FillValue])
ambient_particle_diameter.coordinates = 'latitude longitude'

measurement_channel_lower_limit = NcFile.createVariable(varname='measurement_channel_lower_limit', dimensions=('index',),datatype='float32')
measurement_channel_lower_limit[:] = Size_edge[:-1]
measurement_channel_lower_limit.long_name = 'Lower Limit of Spectrometer Measurement Channel'
measurement_channel_lower_limit.units = 'um'
measurement_channel_lower_limit.valid_min= np.min(Size_edge[:-1])
measurement_channel_lower_limit.valid_max= np.max(Size_edge[:-1])
measurement_channel_lower_limit.coordinates = 'latitude longitude'

measurement_channel_upper_limit = NcFile.createVariable(varname='measurement_channel_upper_limit', dimensions=('index',),datatype='float32')
measurement_channel_upper_limit[:] = Size_edge[1:]
measurement_channel_upper_limit.long_name = 'Upper Limit of Spectrometer Measurement Channel'
measurement_channel_upper_limit.units = 'um'
measurement_channel_upper_limit.valid_min= np.min(Size_edge[1:])
measurement_channel_upper_limit.valid_max= np.max(Size_edge[1:])
measurement_channel_upper_limit.coordinates = 'latitude longitude'

if NormFlag == 0 : 
    ambient_particle_number_per_channel = NcFile.createVariable(varname='ambient_particle_number_per_channel', dimensions=('time','index'),datatype='float32', fill_value = FillValue) 
    ambient_particle_number_per_channel[:,:] = Psd
    ambient_particle_number_per_channel.long_name = 'Ambient Particle Number per Channel (dN)'
    ambient_particle_number_per_channel.units = 'cm-3'
    ambient_particle_number_per_channel.valid_min= np.min(Psd[Psd!=FillValue])
    ambient_particle_number_per_channel.valid_max= np.max(Psd[Psd!=FillValue])
    ambient_particle_number_per_channel.coordinates = 'latitude longitude'
    ambient_particle_number_per_channel.cell_methods = 'time: mean'
if NormFlag == 1:
    ambient_particle_size_distribution = NcFile.createVariable(varname='ambient_particle_size_distribution', dimensions=('time','index'),datatype='float32', fill_value = FillValue) 
    ambient_particle_size_distribution[:,:] = Psd
    ambient_particle_size_distribution.long_name = 'Ambient Particle Size Distribution (dN\dD)'
    ambient_particle_size_distribution.units = 'cm-3 um-1'
    ambient_particle_size_distribution.valid_min= np.min(Psd[Psd!=FillValue])
    ambient_particle_size_distribution.valid_max= np.max(Psd[Psd!=FillValue])
    ambient_particle_size_distribution.coordinates = 'latitude longitude'
    ambient_particle_size_distribution.cell_methods = 'time: mean'  

number_of_instrument_counts_per_channel = NcFile.createVariable(varname='number_of_instrument_counts_per_channel', dimensions=('time','index'),datatype='float32', fill_value = FillValue)  
number_of_instrument_counts_per_channel[:,:] = Counts
number_of_instrument_counts_per_channel.long_name = 'Number of Instrument Counts per Channel'
number_of_instrument_counts_per_channel.units = '1'
number_of_instrument_counts_per_channel.valid_min= np.min(number_of_instrument_counts_per_channel[number_of_instrument_counts_per_channel!=FillValue])
number_of_instrument_counts_per_channel.valid_max= np.max(number_of_instrument_counts_per_channel[number_of_instrument_counts_per_channel!=FillValue])
number_of_instrument_counts_per_channel.coordinates = 'latitude longitude'
number_of_instrument_counts_per_channel.cell_methods = 'time: point'

qc_flag_ambient_particle_number_per_channel = NcFile.createVariable(varname='qc_flag_ambient_particle_number_per_channel', dimensions=('time','index'),datatype='u4')  
qc_flag_ambient_particle_number_per_channel[:,:] = DataFlag_2D
qc_flag_ambient_particle_number_per_channel.long_name = 'Data Quality flag: Ambient Particle Number per Channel'
qc_flag_ambient_particle_number_per_channel.flag_values = '0b, 1b, 2b, 3b'
qc_flag_ambient_particle_number_per_channel.flag_meanings = 'not_used \n\rvalid_data \n\rreduced_quality_data \n\rmissing_data'

# Core variables 
altitude = NcFile.createVariable(varname='altitude', dimensions=('time',),datatype='float32', fill_value = FillValue) 
altitude[:] = ALT_GIN_nc
altitude.long_name = 'Geometric height above geoid (WGS84)'
altitude.units = 'm'
altitude.axis = 'Z'
altitude.valid_min = np.min(ALT_GIN_nc[ALT_GIN_nc!=FillValue])
altitude.valid_max = np.max(ALT_GIN_nc[ALT_GIN_nc!=FillValue])
altitude.cell_methods = 'time:point'

latitude = NcFile.createVariable(varname='latitude', dimensions=('latitude',),datatype='float32', fill_value = FillValue) 
latitude[:] = LAT_GIN_nc
latitude.long_name = 'Latitude'
latitude.units = 'degrees_north'
latitude.axis = 'Y'
latitude.valid_min = np.min(LAT_GIN_nc[LAT_GIN_nc!=FillValue])
latitude.valid_max =  np.max(LAT_GIN_nc[LAT_GIN_nc!=FillValue])
latitude.cell_methods= 'time:point'

longitude = NcFile.createVariable(varname='longitude', dimensions=('longitude',),datatype='float32', fill_value = FillValue) 
longitude[:] = LON_GIN_nc
longitude.long_name = 'Longitude'
longitude.units = 'degrees_east'
longitude.axis= 'X'
longitude.valid_min=  np.min(LON_GIN_nc[LON_GIN_nc!=FillValue])
longitude.valid_max= np.min(LAT_GIN_nc[LAT_GIN_nc!=FillValue])
longitude.cell_methods= 'time:point'

platform_speed_wrt_air = NcFile.createVariable(varname='platform_speed_wrt_air', dimensions=('time',),datatype='float32', fill_value = FillValue)  
platform_speed_wrt_air[:] = TAS_RVSM_nc
platform_speed_wrt_air.long_name = 'Platform speed with respect to air (air speed)' 
platform_speed_wrt_air.units = 'm s-1'
platform_speed_wrt_air.valid_min= np.min(TAS_RVSM_nc[TAS_RVSM_nc!=FillValue])
platform_speed_wrt_air.valid_max= np.max(TAS_RVSM_nc[TAS_RVSM_nc!=FillValue])
platform_speed_wrt_air.cell_methods = 'time: mean'
platform_speed_wrt_air.coordinates = 'latitude longitude'

platform_course = NcFile.createVariable(varname='platform_course', dimensions=('time',),datatype='float32', fill_value = FillValue) 
platform_course[:] = TRCK_GIN_nc
platform_course.long_name = 'Direction in which the platform is travelling'
platform_course.units = 'degree'
platform_course.valid_min= np.min(TRCK_GIN_nc[TRCK_GIN_nc!=FillValue])
platform_course.valid_max= np.max(TRCK_GIN_nc[TRCK_GIN_nc!=FillValue])
platform_course.cell_methods = 'time: mean'
platform_course.coordinates = 'latitude longitude'

platform_speed_wrt_ground = NcFile.createVariable(varname='platform_speed_wrt_ground', dimensions=('time',),datatype='float32', fill_value = FillValue) 
platform_speed_wrt_ground[:] = GSPD_GIN_nc
platform_speed_wrt_ground.long_name = 'Platform speed with respect to ground (ground speed)'
platform_speed_wrt_ground.units = 'm s-1'
platform_speed_wrt_ground.valid_min= np.min(GSPD_GIN_nc[GSPD_GIN_nc!=FillValue])
platform_speed_wrt_ground.valid_max= np.max(GSPD_GIN_nc[GSPD_GIN_nc!=FillValue])
platform_speed_wrt_ground.cell_methods = 'time: mean'
platform_speed_wrt_ground.coordinates = 'latitude longitude'

platform_pitch_angle = NcFile.createVariable(varname='platform_pitch_angle', dimensions=('time',),datatype='float32', fill_value = FillValue) 
platform_pitch_angle[:] = PTCH_GIN_nc
platform_pitch_angle.long_name = 'Pitch Angle'
platform_pitch_angle.units = 'degree'
platform_pitch_angle.valid_min= np.min(PTCH_GIN_nc[PTCH_GIN_nc!=FillValue])
platform_pitch_angle.valid_max= np.max(PTCH_GIN_nc[PTCH_GIN_nc!=FillValue])
platform_pitch_angle.cell_methods = 'time: mean'
platform_pitch_angle.coordinates = 'latitude longitude'

platform_pitch_rate = NcFile.createVariable(varname='platform_pitch_rate', dimensions=('time',),datatype='float32', fill_value = FillValue) 
platform_pitch_rate[:] = PITR_GIN_nc
platform_pitch_rate.long_name = 'Pitch Angle Rate of Change'
platform_pitch_rate.units = 'degree s-1'
platform_pitch_rate.valid_min= np.min(PITR_GIN_nc[PITR_GIN_nc!=FillValue])
platform_pitch_rate.valid_max= np.max(PITR_GIN_nc[PITR_GIN_nc!=FillValue])
platform_pitch_rate.cell_methods = 'time: mean'
platform_pitch_rate.coordinates = 'latitude longitude'

platform_yaw_angle = NcFile.createVariable(varname='platform_yaw_angle', dimensions=('time',),datatype='float32', fill_value = FillValue) 
platform_yaw_angle[:] = HDG_GIN_nc
platform_yaw_angle.long_name = 'Yaw Angle'
platform_yaw_angle.units = 'degree'
platform_yaw_angle.valid_min= np.min(HDG_GIN_nc[HDG_GIN_nc!=FillValue])
platform_yaw_angle.valid_max= np.max(HDG_GIN_nc[HDG_GIN_nc!=FillValue])
platform_yaw_angle.cell_methods = 'time: mean'
platform_yaw_angle.coordinates = 'latitude longitude'

platform_yaw_rate = NcFile.createVariable(varname='platform_yaw_rate', dimensions=('time',),datatype='float32', fill_value = FillValue) 
platform_yaw_rate[:] = HDGR_GIN_nc
platform_yaw_rate.long_name = 'Yaw Angle Rate of Change'
platform_yaw_rate.units = 'degree s-1'
platform_yaw_rate.valid_min= np.min(HDGR_GIN_nc[HDGR_GIN_nc!=FillValue])
platform_yaw_rate.valid_max= np.max(HDGR_GIN_nc[HDGR_GIN_nc!=FillValue])
platform_yaw_rate.cell_methods = 'time: mean'
platform_yaw_rate.coordinates = 'latitude longitude'

DayOfYear = NcFile.createVariable(varname='day_of_year', dimensions=('time',),datatype='float32') 
DayOfYear[:] = day_of_year
DayOfYear.long_name = 'Day of Year'
DayOfYear.units = '1'
DayOfYear.valid_min = np.float32(np.nanmin(day_of_year))
DayOfYear.valid_max = np.float32(np.nanmax(day_of_year))

Year = NcFile.createVariable(varname='year', dimensions=('time',),datatype='u4') 
Year[:] = year
Year.long_name = 'Year'
Year.units = '1'
Year.valid_min = np.nanmin(year)
Year.valid_max = np.nanmax(year)

Month = NcFile.createVariable(varname='month', dimensions=('time',),datatype='u4') 
Month[:] = month
Month.long_name = 'Month'
Month.units = '1'
Month.valid_min = np.nanmin(month)
Month.valid_max = np.nanmax(month)

Day = NcFile.createVariable(varname='day', dimensions=('time',),datatype='u4') 
Day[:] = day
Day.long_name = 'Day'
Day.units = '1'
Day.valid_min = np.nanmin(day)
Day.valid_max = np.nanmax(day)

Hour = NcFile.createVariable(varname='hour', dimensions=('time',),datatype='u4') 
Hour[:] = hour
Hour.long_name = 'Hour'
Hour.units = '1'
Hour.valid_min = np.nanmin(hour)
Hour.valid_max = np.nanmax(hour)

Minute = NcFile.createVariable(varname='minute', dimensions=('time',),datatype='u4') 
Minute[:] = minute
Minute.long_name = 'Minute'
Minute.units = '1'
Minute.valid_min = np.nanmin(minute)
Minute.valid_max = np.nanmax(minute)

Second = NcFile.createVariable(varname='second', dimensions=('time',),datatype='float32') 
Second[:] = second
Second.long_name = 'Second'
Second.units = '1'
Second.valid_min = np.nanmin(second)
Second.valid_max = np.nanmax(second)

# LAT_GIN  # latitude
# LON_GIN  # longitude
# ALT_GIN # altitude
# TAS_RVSM # platform_speed_wrt_air
# TRCK_GIN # platform_course
# # platform_orientation
# GSPD_GIN # platform_speed_wrt_ground
# PTCH_GIN # platform_pitch_angle
# PITR_GIN # platform_pitch_rate
# ROLL_GIN # platform_roll_angle
# ROLR_GIN # platform_roll_rate
# HDG_GIN  # platform_yaw_angle
# HDGR_GIN # platform_yaw_rate


NcFile.close()






