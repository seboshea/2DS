# 2DS

Procedures to identify stereo imaged particles and calculate particle size distributions.


Instructions:
#1) OASIS igor software - raw file processing (roi extraction) and Image Output as .h5 file. 
#2) Add data paths and probe settings to GetFlightInfo2DS().
#3) Info2DS = GetFlightInfo2DS()
#4) BatchBothChannels(Info2DS,FlightNumberStr), find stereo particles and create PSDs using stereo and traditional methods.
#5) HybridStereoProcessing(Info2DS,FlightNumberStr), combine traditional and stereo psds for all files in folder.

Reference:
https://amt.copernicus.org/articles/14/1917/2021/amt-14-1917-2021.html
