%% Initialize and scale BedMachine
% This code reads in BedMachine V3 data for Greenland (stored locally) and scales it down to 1 km grid from 150 m grid of the original data. For the sake of data size, the scaled output is saved to a .MAT file which takes the full 2.2 GB .nc file with errors at 150 m resolution and creates a 9.9 MB .MAT file without erros at 1 km resolution. For future testing and stream lining code, this should make it more portable and speed things up a little.

ncid =netcdf.open('BedMachineGreenland-2017-09-20.nc');
sc = 0.15; % This is the scaling factor to increase grid size to 1 km
surface = imresize(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'surface')),sc);
bed = imresize(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'bed')),sc);
x = imresize(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'x')),sc);
y = imresize(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'y')),sc);
netcdf.close(ncid)
save('scaledBMV3.mat','bed','surface','x','y')
