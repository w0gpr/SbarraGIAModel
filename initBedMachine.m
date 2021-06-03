%% Initialize and scale BedMachine
% This code reads in BedMachine V3 data for Greenland (stored locally) and scales it down to 1 km grid from 150 m grid of the original data. For the sake of data size, the scaled output is saved to a .MAT file which takes the full 2.2 GB .nc file with errors at 150 m resolution and creates a 9.9 MB .MAT file without erros at 1 km resolution. For future testing and stream lining code, this should make it more portable and speed things up a little.
% {'Morlighem M. et al., (2017), BedMachine v3: Complete bed topography and ocean bathymetry mapping of Greenland from multi-beam echo sounding combined with mass conservation, Geophys. Res. Lett., 44, doi:10.1002/2017GL074954. (http://onlinelibrary.wiley.com/doi/10.1002/2017GL074954/full)'}
ncid =netcdf.open('BedMachineGreenland-2017-09-20.nc');
% gridSize = netcdf.getAtt(ncid,netcdf.inqAttID(ncid,'spacing')); % This is
% the wrong way to use this call. The grid size is 150 as directly read
% from the file...
gridSize = 150; % meters
newGridSize = 1000; % meters
sc = gridSize/newGridSize; % This is the scaling factor to increase grid size to 1 km
surface = imresize(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'surface')),sc);
bed = imresize(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'bed')),sc);
x = imresize(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'x')),sc);
y = imresize(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'y')),sc);
netcdf.close(ncid)

save('scaledBMV3.mat','bed','surface','x','y','newGridSize','sc')
