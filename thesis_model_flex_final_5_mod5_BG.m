%% Example 1: Simplest case

% This commented out section was separated to initBedMachine.m
% The simplest way to load BedMachine data is to just specify which variable
% you want to load, like this:
% ncid = netcdf.open('BedMachineGreenland-2017-09-20.nc');
% surface = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'surface')); 
% x = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'x')); 
% y = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'y')); 
% bed = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'bed')); 
% netcdf.close(ncid);
% sc = 150/1000; % scale for resizing the grid from 150 m to 1 km grid size
% Small_surf_x = imresize(x,sc);
% Small_surf_y = imresize(y,sc);
% Small_surf =imresize(surface,sc);
% x = imresize(x,sc); % Should use this for referencing point inputs
% y = imresize(y,sc);
% Small_bed =imresize(bed,sc);

clear
BM = open('scaledBMV3.mat');
Small_surf = BM.surface';
x = BM.x; % Should use this for referencing point inputs
y = BM.y;
Small_bed = BM.bed';
dx = BM.newGridSize;
clear BM

% %% Pick the points the line should pass through
% This section should be for user input of points of interest and/or for
% automatic point generation. The user should pick 2 points (either by
% entering lat/long in decimal degrees, or by picking on a map. The code
% should then translate those points to: Polar Stereographic North via
% ll2psn, to the indices in x and y that coorespond to those points on Bed
% Machine. A conditional should be added to check that those points
% actually fall within the domain.

% The original program had a lot of direct indices referenced in the code.

% These are to become user supplied points to do the comparison along. At
% the moment it is a static entry to allow for testing.
point1 = [66.953476, -50.873059];
point2 = [66.851245, -44.752894];

points = [point1;point2];
[c,d] = ll2psn(points(:,1),points(:,2));
c = int32(c); d = int32(d); % This forces the values to be not floating nums

% This is a check to make sure the points are not out of bounds This is
% especially helpful for a direct entry format rather than for future plans
% with points that are chosen based on clicking on a map.

for i = 1:length(c)
    if c(i)<x(1)||c(i)>x(end)
        error(['Longitude of point ',num2str(i),' out of bound'])
    end
    if d(i)>y(1)||d(i)<y(end)
        error(['Latitude of point ',num2str(i),' out of bound'])
    end
end

xI = zeros(2,1);
yI = zeros(2,1);
for i = 1:2
    [~,xI(i)] = min(abs(x-c(i)));
    [~,yI(i)] = min(abs(y-d(i)));
end


% This set of if statements tests to make sure the line is not vertical. If
% it is not vertical, then it solves for the x and y boundaries of the line
% projected through the 2 chosen points. If the ends of the line fall on
% the North and South edges of the frame, then the nested if is activated
% and solves for the points that way. If the line is vertical, then the
% least squares would fail due to a singular matric, and the bounds are 
% identified directly.

if diff(xI) ~= 0
    B = [ones(size(xI)) xI]\yI;
    xbound = [1 length(x)];
    ybound = xbound.*B(2)+B(1);
    ybound = round(ybound);
    if ybound(1)<0 || ybound(2)<0
        ybound = [1 length(y)];
        xbound = (ybound - B(1))./B(2);
        xbound = round(xbound);
    end
else
%     disp('derp')
    xbound = xI;
    ybound = [1 length(y)];
end
    

% %% Chris's original setup for drawing the line hardcoded
% % Study 66.986881,-53.733691 (1864 and 1869) Sis  66.951285,-53.725157 (1869)
% %kang  66.953476, -50.873059 (1885) Disco  69.451704,-52.581850 (1589)
% %aasiat  68.713605,-52.801423 (1679)
% [a,b] = ll2psn(point1(1),point1(2));
% Lat = find(y < b);% x and a           y and b
% lon = [1864 1869];  % These are hard coded locations and should instead be input and converted
% lat = [269 300];
% slope = ((lon(2)-lon(1))/((lat(2)-lat(1))));
% lon2 =[(lon(1)-(1533*slope*(lat(1)/1533))),(lon(2)+(1533*slope*((1533-lat(2))/1533)))];
% lat2=[1 1533];
% x1=lat2(1);x2=lat2(2);y1=round(lon2(1));y2=round(lon2(2));



figure(1)
clf
hold on
imagesc(x,y,Small_bed)
title('Bedmachine Bedrock topography under the Current Greenland Ice Sheet');
% set(gca,'YDir','reverse')
xlabel('distance (Km)');ylabel('distance(Km)');
line(x(xbound),y(ybound))
plot(c,d,'or')
% line(x([x1,x2]),y([y1,y2]))
hold off

% Sis_bed=improfile(Small_bed,lat2,lon2).';% if doing a complex profile
% Sis_surf = improfile(Small_surf,lat2,lon2).';

Sis_bed=improfile(Small_bed,xbound,ybound).';% if doing a complex profile
Sis_surf = improfile(Small_surf,xbound,ybound).';

% This next part is to separate the ice sheet in half for the ice divide.
% Ideally this would use a vector or raster to pick the intersection of the
% line of interest (determined earlier) and the ice divide. It could also
% be done using the slope of the ice sheet. I don't think max can be used
% automatically in the event that a mountain on the east side is higher
% than the ice divide. For now however, I'll use the max peak.
[maxIceElevation, iceDivide] = max(Sis_surf);

% AAA = (Sis_surf(1,435:725)-Sis_bed(1,435:725));
%Sis_bed = Small_bed(1869,:); %horizontal profile
TKM = iceDivide; % total km of the study path
SKM = 495; %total km for small ice path

figure(2)
clf
plot((x(1:TKM)/1000),(Sis_bed(1:TKM)))
hold on
plot((x(1:TKM)/1000),(Sis_surf(1:TKM)))
hold off

figure(3)
clf
plot(x/1000,[Sis_bed;Sis_surf])
hold on
plot(x(iceDivide)/1000,maxIceElevation,'ok')
hold off


% west_x = (0:1000:724000); 
west_x = (0:dx:(iceDivide-1)*dx); % west half of ice sheet

% Use this to zoom in to pick the LGM max location. This value would be
% better suited as an input or "automatically" picked based on the line
% generated above.
figure(4)
plot(1:iceDivide,Sis_bed(1:iceDivide))
%% max flex from 0

% dx = 1000;    % [m] choose some suitably small delta x
contShelf = 152; % this is the continental shelf edge 'eyballed' from fig4
iceRetreat = 607;

% Lmax could be set earlier as a parameter input
Lmax  = (iceDivide-contShelf)*1000*2; % [m] full width of the ice sheet at maximum (edge to edge)

% horizontal space
xk = -5*Lmax : dx : 5*Lmax;


% Lmin = (725-455)*1000*2; % ]m] full width of the ice sheet at minimum (min margin)
% Lmarg = Lmax:-dx:Lmin;
% up = zeros(607,11461);
% surf = zeros(1,11461);
% Sis_surf_flip= importdata('Sis_surf_flip.csv')';
% surf(5441:6021) =Sis_surf_flip; 

% constants
constants

% These are defined in the file constants.m
% icetau = 110e3;                     % [Pa], yield strength of ice 100 kPa
% rho_ice = 917;                      % [kg/m3] ice density
% rho_mantle = 3200;                  % [kg/m3] mantle density
% gk = 9.8;                           % [m/s2] gravity
% hk =40e3;                        	% [m] Lithosphere elastic thickness: Te ~ h ~ 100 km for cratons
% Ek = 70e9;                          % [Pa] Young's Modulus, Turcotte & Schubert p. 152 Problem 3.19
% nuk = 0.25;                         % [] Poisson's Ratio, p. 152 Problem 3.19
% Dk = Ek*hk^3/12/(1-nuk^2);          % [Pa m^3] Flexural rigidity parameter for lithoplate
% alpha = (4*Dk/rho_mantle/gk)^0.25;  % [m] length scale of bulging

% This for loop builds the ice sheet profile for each km step from the
% continental divide to the ice minima

up = zeros(iceRetreat,Lmax/100+1);
for i = 1:iceRetreat   
    Lk = (Lmax:-dx:Lmax-((i-1)*dx))';
%     Lk = Lk';

    % parabolic ice sheet profile
    Hk = real(sqrt(2*icetau / rho_ice / gk * (Lk(i)/2-xk)));
    Hk(xk<0) = fliplr(Hk(xk>0));    % This makes it symmetrical about the ice divide
    up(i,:) = Hk(:);
end

Hk = up(1:2:end,:); % This decimates the profile to simplify the choice of profiles
%%
% here you put in you for loop for creating the animation
% first initialize your output matrices to save speed
% then run the for loop

% load V(x)
nx = length(xk);
out = zeros(304, nx);
for i = 1:304
Vk = Hk(i,:)*dx;  % [Pa/m] spatially varying line load V(x)

% Flexure from a generic load
gensol = plate1bvp_v2(false);  % false will not make any plots
% gensol.x is x in units of alpha
% gensol.y is w in units of (V_0 alpha^3) / (8 D)
xscaled = gensol.x * alpha;  % now this is in meters

% Step through all the line loads V and scale the flexure solution to them.
% They can then be superimposed to get the full flexure from the ice sheet.
wk = zeros(1,nx);  % w is the total vertical displacement in plate
% Step through x space to get all the loads
% figure(12); clf; hold on
for ii=1:nx
    Vk = rho_ice * gk * Hk(i,ii) * dx;
    if Vk>0  % only do if there's a load here
        % find the displacement from the load centered here
        yscaled = -gensol.y(1,:) * Vk * alpha^3 / (8*Dk);
        % shift the line load to its location in x
        xshifted = xk - xk(ii);
        % the vertical displacement caused by this individual load: 
        wii = interp1([-fliplr(xscaled) xscaled(2:end)],[fliplr(yscaled) yscaled(2:end)],xshifted);
        % add it in to the total displacement w: 
        wk = wk + wii;
    end
end
out(i,:) = wk(:);
end

out1 = out;
out1 = out1-(out1(284,:));
valk = 1:304;
flex = zeros(725,725);
for i = 1:145 %(Big ice)
    flex(i,1:725) = out1(1,5007:5731);
end
% ice retrat over time and distance 
flex(146:261,1:725) = out1(1:116,5007:5731);%LGM to 12ka
flex(262:278,1:725) =out1(116:2:149,5007:5731);% 12ka to 11.6 ka
rep1 = 150:1:(150+23);
flex(279:326,1:725) = out1(repelem(rep1,2),5007:5731);% 11.6 to 10.4
flex(327:376,1:725) = out1(174:223,5007:5731);% 10.4 to 9.1 ka
flex(377:403,1:725) = out1(224:250,5007:5731);% 9.1 to 8.3 ka
rep2 = 251:1:(251+16);
flex(404:437,1:725) = out1(repelem(rep2,2),5007:5731);% 8.3 to 7.5 ka
flex(438:445,1:725) = out1(268:275,5007:5731); %7.5 to 7.3 ka
rep3 = 276:1:(276+27);
flex(446:557,1:725) = out1(repelem(rep3,4),5007:5731);%7.3 to 4.3 ka
for i = 558:649
    flex(i,1:725) = out1(304,5007:5731); %4.3 to 2 ka
end 
rep4 = 303:-1:285;
flex(650:725,1:725) = out1(repelem(rep4,4),5007:5731);%2ka to current

for i = 1:725
    flex_real(i, :) = flex(i,:) + Sis_bed(1,1:725); % sets area to bed topo
end
% small ice 
flexs = zeros(725,725);
for i = 1:146 %(Big ice)
    flexs(i,1:725) = out1(79,5007:5731);
end
% ice retrat over time and distance 
rep0 = 80:1:109;
flexs(147:236,1:725) = out1(repelem(rep0,3),5007:5731);%LGM to 12ka
flexs(237,1:725) = flexs(236,1:725);
rep00 = 110:1:115;
flexs(238:261,1:725) = out1(repelem(rep00,4),5007:5731);
flexs(262:278,1:725) =out1(116:2:149,5007:5731);% 12ka to 11.6 ka
rep1 = 150:1:(150+23);
flexs(279:326,1:725) = out1(repelem(rep1,2),5007:5731);% 11.6 to 10.4
flexs(327:376,1:725) = out1(174:223,5007:5731);% 10.4 to 9.1 ka
flexs(377:403,1:725) = out1(224:250,5007:5731);% 9.1 to 8.3 ka
rep2 = 251:1:(251+16);
flexs(404:437,1:725) = out1(repelem(rep2,2),5007:5731);% 8.3 to 7.5 ka
flexs(438:445,1:725) = out1(268:275,5007:5731); %7.5 to 7.3 ka
rep3 = 276:1:(276+27);
flexs(446:557,1:725) = out1(repelem(rep3,4),5007:5731);%7.3 to 4.3 ka
for i = 558:649
    flexs(i,1:725) = out1(304,5007:5731); %4.3 to 2 ka
end 
rep4 = 303:-1:285;
flexs(650:725,1:725) = out1(repelem(rep4,4),5007:5731);%2ka to current

for i = 1:725
    flexs_real(i, :) = flexs(i,:) + Sis_bed(1,1:725); % sets area to bed topo
end
%
iceflex = zeros(725,725);
for i = 1:146
    iceflex(i,1:725) = Hk(1,5007:5731);
end
% iceflex(152:455,1:725) = Hk(1:304,5007:5731);
iceflex(147:261,1:725) = Hk(1:115,5007:5731);%LGM to 12ka
iceflex(262:278,1:725) =Hk(116:2:149,5007:5731);% 12ka to 11.6 ka
rep1 = 150:1:(150+23);
iceflex(279:326,1:725) = Hk(repelem(rep1,2),5007:5731);% 11.6 to 10.4
iceflex(327:376,1:725) = Hk(174:223,5007:5731);% 10.4 to 9.1 ka
iceflex(377:403,1:725) = Hk(224:250,5007:5731);% 9.1 to 8.3 ka
rep2 = 251:1:(251+16);
iceflex(404:437,1:725) = Hk(repelem(rep2,2),5007:5731);% 8.3 to 7.5 ka
iceflex(438:445,1:725) = Hk(268:275,5007:5731); %7.5 to 7.3 ka
rep3 = 276:1:(276+27);
iceflex(446:557,1:725) = Hk(repelem(rep3,4),5007:5731);%7.3 to 4.3 ka
for i = 558:649
    iceflex(i,1:725) = Hk(304,5007:5731); %4.3 to 2 ka
end 
rep4 = 303:-1:285;
iceflex(650:725,1:725) = Hk(repelem(rep4,4),5007:5731);%2ka to current
%% Plotting flexure

position =[100 200 300 400 500 600];
%iceflex(1:position,1:725) = flex_real(1:position,1:725);

figure(2)
clf
hold on
% plot(west_x/(1e3),iceflex(position,:)+flex(position,:))
plot(west_x/1e3,flex(position,:))
plot(west_x/1e3, iceflex(position,:))
% plot(west_x/(1e3),Sis_bed(1,1:TKM))

% xL = get(gca, 'XLim')
% plot(xL, [0 0], 'linewidth',2)
xlabel('kilometers')
ylabel('meters')

% the foredeep is 150km from the front of the ice margin
figure(3)
clf
hold on
time = 18850:-26:26;
plot(time/1000,flex(:,268),'b')
plot(time/1000,flexs(:,268),'r--')
xlabel('time (ka)')


%% big ice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dist_X = (0:1000:1532000);% west to east distance used for plotting the entire ice sheet
west_x = (0:1000:724000); %west half of ice sheet

figure(4)
clf
hold on
plot(west_x,Sis_bed(1,1:TKM))
plot(west_x(435:TKM),Sis_surf(1,435:TKM))
% xlabel('Distance (m)');ylabel('elevation(m)');
title('Current ice profile from Bedmachine surface data');
%%
SEA = zeros(1,725);
rho_ice = 910; % kg/m^3 density ice
gk = 9.81; % gravity
F = 1;
% Ty1 = zeros(270,1)+30000;
% Ty2 =zeros(1263,1)+55000;
% Ty = cat(1,Ty1,Ty2);
edge = 152;     % edge of the continental shelf
END = 435;
flat = zeros(1,TKM);
for rrr = edge:TKM
    flat(1,rrr) = (4.7*(west_x(1,rrr-(edge-1)).^.5));
end

for rr=edge:TKM
Tyy(rr,1) = 55000;%1/2 of basal shear
end
Tyy(TKM-50:TKM,1) = 55000:-(55000/50):0;
fh = zeros(length(west_x),1); % ice surface (m)
fh(1:edge) = Sis_bed(1,1:edge);
fH = zeros(length(west_x),1); % ice thickness (m)
fH(1:edge) = 0;
fB = zeros(length(west_x),1); % b term of quadractic equation
fB(1:edge) = 0;
fC = zeros(length(west_x),1); % c term of quadratic equation
fC(1:edge) = 0;
Tyyy = zeros(TKM,TKM);
for mm = 1:edge
Tyyy (:,mm) = Tyy(:,1)  ;
end
for mmm = edge:TKM
Tyyy(mmm+1:end,mmm+1) = Tyyy(mmm:TKM-1,mmm);
Tyyy(TKM-(50):end,mmm+1) =55000:-(55000/50):0;
end
for j =(2:TKM-(edge-1))
   fB(j+(edge-1)) = -(Sis_bed(1,j+(edge-1))+Sis_bed(1,j+(edge-2)));
   fC(j+(edge-1)) = (fh(j+(edge-2))*(Sis_bed(1,j+(edge-1))-(fH(j+(edge-2)))))-((2*(west_x(1,j+(edge-1))-west_x(1,j+(edge-2)))*(Tyy(j+(edge-1))+Tyy(j+(edge-2))/2))/(rho_ice*gk));
   fh(j+(edge-1)) =(-fB(j+(edge-1))+abs((fB(j+(edge-1))^2)-(4*fC(j+(edge-1))))^0.5)/2;
   fH (j+(edge-1)) = (fh(j+(edge-1))-Sis_bed(j+(edge-1)));
end 
figure(5)
hold on
plot(west_x,fh)
plot(west_x,-140,'.k')
plot(west_x,Sis_bed(1,1:TKM))

% Sis_bed_abs = Sis_bed;
% LGM_SL = find(Sis_bed_abs<-140);
% Sis_bed_abs(LGM_SL) = -140;

Water_change(1,(1:(edge-2))) = Sis_bed(1,(1:(edge-2)))+44.9312;
Water_change(1,(edge-1):1148) = Sis_bed(1,(edge-1):1148); 

figure(6);
hold on
plot(dist_X,Sis_bed)
plot(dist_X(1,1:TKM),Water_change(1,1:TKM))
plot(dist_X,0,'.')
%% polygon
rockcolor = [0.629 0.494 0.125]; % browny
icecolor =  [0.678 0.921 1];     % blueish

polyx  = [west_x fliplr(west_x)];         % x locations of my polygon
polyys = [fh(1:TKM,1).' fliplr(Sis_bed(1,1:TKM))]; % y locations of ice polygon
polyyb = [Sis_bed(1,1:TKM) -1200*ones(size(Sis_bed(1,1:TKM)))];        % y locations of my bed polygon (rock)

    figure(4); clf; hold on
    pp= fill(polyx/1,polyys,icecolor);  % 2003 ice polygon
    p= fill(polyx/1,polyyb,rockcolor);     % bedrock polygon
    p.LineStyle='none';pp.LineStyle='none';
    text(6e5,-500,'Bedrock','fontsize',24)%,'rotation',22.5)
    text(6e5,1000,'Ice Sheet','fontsize',24)%,'rotation',22.5)
    title('Ice sheet and bed')
    xlabel('km from terminus'); ylabel('meters a.s.l.')
    

%% 
h = zeros(length(west_x),1); % ice surface (m)
h(1:edge) = Sis_bed(1,1:edge);
H = zeros(length(west_x),1); % ice thickness (m)
H(1:edge) = 0;
B = zeros(length(west_x),1); % b term of quadractic eqaution
B(1:edge) = 0;
C = zeros(length(west_x),1); % c term of quadratic equation
C(1:edge) = 0;

for j =(2:TKM-edge+1)
   B(j+(edge-1)) = -(Sis_bed(1,j+(edge-1))+Sis_bed(1,j+(edge-2)));
   C(j+(edge-1)) = (h(j+(edge-2))*(Sis_bed(1,j+(edge-1))-(H(j+(edge-2)))))-((2*(west_x(1,j+(edge-1))-west_x(1,j+(edge-2)))*(Tyyy(j+(edge-1),edge)+Tyyy(j+(edge-2),edge)/2))/(rho_ice*gk));
   h(j+(edge-1)) =(-B(j+(edge-1))+abs((B(j+(edge-1))^2)-(4*C(j+(edge-1))))^0.5)/2;
   H (j+(edge-1)) = (h(j+(edge-1))-Sis_bed(j+(edge-1)));
end 
%% Small ice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h2 = zeros(length(west_x),1); % ice surface (m)
h2(1:TKM-SKM) = Sis_bed(1,1:TKM-SKM);
H2 = zeros(length(west_x),1); % ice thickness (m)
H2(1:TKM-SKM) = 0;
B2 = zeros(length(west_x),1); % b term of quadractic eqaution
B2(1:TKM-SKM) = 0;
C2 = zeros(length(west_x),1); % c term of quadratic equation
C2(1:TKM-SKM) = 0;

flats = zeros(1,TKM);
for rrr = TKM-SKM:(TKM)
    flats(1,rrr) = (4.7*(west_x(1,rrr-(TKM-SKM-1)).^.5));
end
Tyys =zeros(TKM,1);
for rr=(TKM-SKM):TKM
Tyys(rr,1) = 55000;%1/2 of basal shear
end
Tyys(TKM-50:TKM,1) = 55000:-(55000/50):0;
Tyyys = zeros(TKM,TKM);
for mm = 1:(TKM-SKM)
Tyyys (:,mm) = Tyys(:,1)  ;
end
for mm = (TKM-SKM):TKM
Tyyys(mm+1:end,mm+1) = Tyyys(mm:TKM-1,mm);
Tyyys(TKM-(50):end,mm+1) =55000:-(55000/50):0;
end
for j =(2:SKM+1)
   B2(j+TKM-SKM-1) = -(Sis_bed(1,j+TKM-SKM-1)+Sis_bed(1,j+TKM-SKM-2));
   C2(j+TKM-SKM-1) = (h2(j+TKM-SKM-2)*(Sis_bed(1,j+TKM-SKM-1)-(H2(j+TKM-SKM-2))))-((2*(west_x(1,j+TKM-SKM-1)-west_x(1,j+TKM-SKM-2))*(Tyyys(j+(TKM-SKM-1),TKM-SKM)+Tyyys(j+(TKM-SKM-2),TKM-SKM)/2))/(rho_ice*gk));
   h2(j+TKM-SKM-1) =(-B2(j+TKM-SKM-1)+abs((B2(j+TKM-SKM-1)^2)-(4*C2(j+TKM-SKM-1)))^0.5)/2;
   H2 (j+TKM-SKM-1) = (h2(j+TKM-SKM-1)-Sis_bed(j+TKM-SKM-1));
end
figure
hold on
plot(west_x,h2)
plot(west_x,-140,'.')
plot(west_x,Sis_bed(1,1:TKM))
plot(west_x,-140,'.')

%% current ice sheet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h3 = zeros(length(west_x),1); % ice surface (m)
h3(1:END) = Sis_bed(1,1:END);
H3 = zeros(length(west_x),1); % ice thickness (m)
H3(1:END) = 0;
B3 = zeros(length(west_x),1); % b term of quadractic eqaution
B3(1:END) = 0;
C3 = zeros(length(west_x),1); % c term of quadratic equation
C3(1:END) = 0;
flatc = zeros(1,TKM);
for rrr = END:(TKM)
    flatc(1,rrr) = (4.7*(west_x(1,rrr-(END-1)).^.5));
end
Tyyc =zeros(TKM,1);
num1 = 55000;
for rr=(END):TKM
Tyyc(rr,1) =num1;
end
Tyyc(TKM-50:TKM,1) = num1:-(num1/50):0;
Tyyyc = zeros(TKM,TKM);
for mm = 1:(END)
Tyyyc (:,mm) = Tyyc(:,1)  ;
end
for j =(2:(TKM-(END-1)))
   B3(j+(END-1)) = -(Sis_bed(1,j+(END-1))+Sis_bed(1,j+(END-2)));
   C3(j+(END-1)) = (h3(j+(END-2))*(Sis_bed(1,j+(END-1))-(H3(j+(END-2)))))-((2*(west_x(1,j+(END-1))-west_x(1,j+(END-2)))*(Tyyyc(j+(END-1),END)+Tyyyc(j+(END-2),END)/2))/(rho_ice*gk));
   h3(j+(END-1)) =(-B3(j+(END-1))+abs((B3(j+(END-1))^2)-((4)*C3(j+(END-1))))^0.5)/2;
   H3 (j+(END-1)) = (h3(j+(END-1))-Sis_bed(j+(END-1)));
end 
H3 = zeros(TKM,1);
H3(END:TKM,1) = Sis_surf(1,435:TKM);
bb=1;
polyx  = [west_x fliplr(west_x)];         % x locations of my polygon
cPolyys(1:1450,1) = [h3(1:TKM,1).' fliplr(Sis_bed(1,1:TKM))]; % y locations of ice polygon
cPolyyb(1:1450,1) = [Sis_bed(1,1:TKM) -1200*ones(size(Sis_bed(1,1:TKM)))];        % y locations of my bed polygon (rock)
cpolywater(bb,1:2*END) = [1:END fliplr(1:END)];
cPolyyw(bb,1:2*END) = [SEA(1,1:END) min(fliplr(Sis_bed(1,1:END)),SEA(1,END))];
figure
hold on    
fill(polyx,cPolyys(:,1),icecolor,'edgecolor','k');hold on;
plot(west_x(1,435:725),Sis_surf(1,435:725),'r','linewidth',2)
hold on; fill(polyx,cPolyyb(:,1),rockcolor);hold on;
p=fill(cpolywater(1,:)*1000,cPolyyw(1,:),b);
p.LineStyle='none';
text(5.5e5,-500,'Bedrock','fontsize',12);
text(5.5e5,1000,'Ice Sheet','fontsize',12);
text(5e4,-500,'Ocean','fontsize',8);
title('Modeled Current Ice Sheet', 'FontSize', 12);
xlabel('Distance (km)');ylabel('Elevation (m)');
ax = gca;ax.XLim = [0 TKM*1000];ax.YLim = [-1200 3500];
xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/1000) 
hLegend = legend({ 'Modeled output of current ice sheet','Elevation data of current ice sheet'}, 'Location', 'northwest','fontsize',12);

polyx  = [west_x fliplr(west_x)];         % x locations of my polygon
ccPolyys(1:1450,1) = [Sis_surf(1,1:TKM) fliplr(Sis_bed(1,1:TKM))]; % y locations of ice polygon
ccPolyyb(1:1450,1) = [Sis_bed(1,1:TKM) -1200*ones(size(Sis_bed(1,1:TKM)))];        % y locations of my bed polygon (rock)
ccpolywater(bb,1:2*END) = [1:END fliplr(1:END)];
ccPolyyw(bb,1:2*END) = [SEA(1,1:END) min(fliplr(Sis_bed(1,1:END)),SEA(1,END))];

figure
hold on    
fill(polyx,ccPolyys(:,1),icecolor,'edgecolor','k');hold on;
hold on; fill(polyx,ccPolyyb(:,1),rockcolor);hold on;
p=fill(ccpolywater(1,:)*1000,ccPolyyw(1,:),b);
p.LineStyle='none';
text(5.5e5,-500,'Bedrock','fontsize',12);
text(5.5e5,1000,'Ice Sheet','fontsize',12);
text(5e4,-500,'Ocean','fontsize',8);
title('True Current Ice Sheet', 'FontSize', 12);
xlabel('Distance (km)');ylabel('Elevation (m)');
ax = gca;ax.XLim = [0 TKM*1000];ax.YLim = [-1200 3500];
xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/1000) 

[hval,pval,ci,stats] = ttest2(Sis_surf(1,435:725),h3(435:725,1))
%% big ice with depression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B_LGM_Dep = (H-(H3)).*(.91/3.2);
%Big_LGM_bed = Sis_bed(1,1:TKM)-B_LGM_Dep.';
Big_LGM_bed = flex_real(152,:);
%Big_LGM_bed = B_LGM_Dep.'-Sis_bed(1,1:766);
%Big_LGM_bed(1,1:152) = Water_change(1,1:152);
figure
hold on
plot(west_x,Water_change(1,1:TKM))
plot(west_x,Sis_bed(1,1:TKM))
plot(west_x,Big_LGM_bed)

hb = zeros(length(west_x),1); % ice surface (m)
hb(1:edge) = Big_LGM_bed(1,1:edge);
Hb = zeros(length(west_x),1); % ice thickness (m)
Hb(1:edge) = 0;
Bb = zeros(length(west_x),1); % b term of quadractic eqaution
Bb(1:edge) = 0;
Cb = zeros(length(west_x),1); % c term of quadratic equation
Cb(1:edge) = 0;

for j =(2:(TKM-edge+1))
   Bb(j+(edge-1)) = -(Big_LGM_bed(1,j+(edge-1))+Big_LGM_bed(1,j+(edge-2)));
   Cb(j+(edge-1)) = (hb(j+(edge-2))*(Big_LGM_bed(1,j+(edge-1))-(Hb(j+(edge-2)))))-((2*(west_x(1,j+(edge-1))-west_x(1,j+(edge-2)))*(Tyyy(j+(edge-1),edge)+Tyyy(j+(edge-2),edge)/2))/(rho_ice*gk));
   hb(j+(edge-1)) =(-Bb(j+(edge-1))+abs((Bb(j+(edge-1))^2)-(4*Cb(j+(edge-1))))^0.5)/2;
   Hb (j+(edge-1)) = (hb(j+(edge-1))-Big_LGM_bed(j+(edge-1)));
end 

figure
hold on
plot(west_x,hb)
plot(west_x,-140,'.')
plot(west_x,Big_LGM_bed(1,1:TKM))
plot(west_x,-140,'.')
%% depression iteration
% B_LGM_Dep1 = (Hb-(H3)).*(.91/3.2);
% Big_LGM_bed1 = Sis_bed(1,1:TKM)-B_LGM_Dep1.';
%Big_LGM_bed = B_LGM_Dep.'-Sis_bed(1,1:766);
%Big_LGM_bed1(1,1:152) = Water_change(1,1:152);
Big_LGM_bed1 = flex_real(152,:);
figure
hold on
plot(west_x,Water_change(1,1:TKM))
plot(west_x,Sis_bed(1,1:TKM))
plot(west_x,Big_LGM_bed1)

hb1 = zeros(length(west_x),1); % ice surface (m)
hb1(1:edge) = Big_LGM_bed1(1,1:edge);
Hb1 = zeros(length(west_x),1); % ice thickness (m)
Hb1(1:edge) = 0;
Bb1 = zeros(length(west_x),1); % b term of quadractic eqaution
Bb1(1:edge) = 0;
Cb1 = zeros(length(west_x),1); % c term of quadratic equation
Cb1(1:edge) = 0;

for j =(2:TKM-edge+1)
   Bb1(j+(edge-1)) = -(Big_LGM_bed1(1,j+(edge-1))+Big_LGM_bed1(1,j+(edge-2)));
   Cb1(j+(edge-1)) = (hb1(j+(edge-2))*(Big_LGM_bed1(1,j+(edge-1))-(Hb1(j+(edge-2)))))-((2*(west_x(1,j+(edge-1))-west_x(1,j+(edge-2)))*(Tyyy(j+(edge-1),edge)+Tyyy(j+(edge-2),edge)/2))/(rho_ice*gk));
   hb1(j+(edge-1)) =(-Bb1(j+(edge-1))+abs((Bb1(j+(edge-1))^2)-(4*Cb1(j+(edge-1))))^0.5)/2;
   Hb1 (j+(edge-1)) = (hb1(j+(edge-1))-Big_LGM_bed1(j+(edge-1)));
end 
figure
hold on
plot(west_x,hb1)
plot(west_x,-140,'.')
plot(west_x,Big_LGM_bed1(1,1:TKM))
plot(west_x,-140,'.')
%% Small ice with depression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S_LGM_Dep = (H2-H3).*(.91/3.2);
% Small_LGM_bed = Sis_bed(1,1:TKM)-S_LGM_Dep.';
%Small_LGM_bed(1,1:152) = Water_change(1,1:152);
Small_LGM_bed = flex_real(230, :);
figure
hold on
plot(west_x,Small_LGM_bed)

hs = zeros(length(west_x),1); % ice surface (m)
hs(1:TKM-SKM) = flex_real(230,1:TKM-SKM);
Hs = zeros(length(west_x),1); % ice thickness (m)
Hs(1:TKM-SKM) = 0;
Bs = zeros(length(west_x),1); % b term of quadractic eqaution
Bs(1:TKM-SKM) = 0;
Cs = zeros(length(west_x),1); % c term of quadratic equation
Cs(1:TKM-SKM) = 0;

for j =(2:SKM+1)
   Bs(j+TKM-SKM-1) = -(Small_LGM_bed(1,j+TKM-SKM-1)+Small_LGM_bed(1,j+TKM-SKM-2));
   Cs(j+TKM-SKM-1) = (hs(j+TKM-SKM-2)*(Small_LGM_bed(1,j+TKM-SKM-1)-(Hs(j+TKM-SKM-2))))-((2*(west_x(1,j+TKM-SKM-1)-west_x(1,j+TKM-SKM-2))*(Tyyys(j+(TKM-SKM-1),TKM-SKM)+Tyyys(j+(TKM-SKM-2),TKM-SKM)/2))/(rho_ice*gk));
   hs(j+TKM-SKM-1) =(-Bs(j+TKM-SKM-1)+abs((Bs(j+TKM-SKM-1)^2)-(4*Cs(j+TKM-SKM-1)))^0.5)/2;
   Hs (j+TKM-SKM-1) = (hs(j+TKM-SKM-1)-Small_LGM_bed(j+TKM-SKM-1));
end
figure
hold on
plot(west_x,hs)
plot(west_x,-140,'.')
plot(west_x,Small_LGM_bed(1,1:TKM))
plot(west_x,-140,'.')
%% Small depression iteration
% S_LGM_Dep1 = (Hs-(H3)).*(.91/3.2);
% Small_LGM_bed1 = Sis_bed(1,1:TKM)-S_LGM_Dep1.';
Small_LGM_bed1 = flex_real(230,:);
figure
hold on
plot(west_x,Water_change(1,1:TKM))
plot(west_x,Sis_bed(1,1:TKM))
plot(west_x,Small_LGM_bed1)

hs1 = zeros(length(west_x),1); % ice surface (m)
hs1(1:TKM-SKM) = flex_real(230,TKM-SKM);
Hs1 = zeros(length(west_x),1); % ice thickness (m)
Hs1(1:TKM-SKM) = 0;
Bs1 = zeros(length(west_x),1); % b term of quadractic eqaution
Bs1(1:TKM-SKM) = 0;
Cs1 = zeros(length(west_x),1); % c term of quadratic equation
Cs1(1:TKM-SKM) = 0;

for j =(2:SKM+1)
   Bs1(j+TKM-SKM-1) = -(Small_LGM_bed1(1,j+TKM-SKM-1)+Small_LGM_bed1(1,j+TKM-SKM-2));
   Cs1(j+TKM-SKM-1) = (hs1(j+TKM-SKM-2)*(Small_LGM_bed1(1,j+TKM-SKM-1)-(Hs1(j+TKM-SKM-2))))-((2*(west_x(1,j+TKM-SKM-1)-west_x(1,j+TKM-SKM-2))*(Tyyys(j+(TKM-SKM-1),TKM-SKM)+Tyyys(j+(TKM-SKM-2),TKM-SKM)/2))/(rho_ice*gk));
   hs1(j+TKM-SKM-1) =(-Bs1(j+TKM-SKM-1)+abs((Bs1(j+TKM-SKM-1)^2)-(4*Cs1(j+TKM-SKM-1)))^0.5)/2;
   Hs1 (j+TKM-SKM-1) = (hs1(j+TKM-SKM-1)-Small_LGM_bed1(j+TKM-SKM-1));
end
figure
hold on
plot(west_x,hs1)
% plot(west_x,-140,'.')
% plot(west_x,Small_LGM_bed1(1,1:TKM))
% plot(west_x,-140,'.')
title('Funder et al 2011 LGM Ice sheet');

% figure
% hold on
% plot(west_x,hb1)
% plot(west_x,-140,'.')
% plot(west_x,Big_LGM_bed1(1,1:TKM))
% plot(west_x,-140,'.')

%% Retreating ice (change shear stress in big ice and hit run) 
hb2 = zeros(length(west_x),length(west_x)); % ice surface (m)
Hb2 = zeros(length(west_x),length(west_x)); % ice thickness (m)
Hb2(1:TKM,1) = 0;
Bb2 = zeros(length(west_x),length(west_x)); % b term of quadractic eqaution
Cb2 = zeros(length(west_x),length(west_x)); % c term of quadratic equation
Cb2 (1:TKM,1) = 0;



for u=(1:TKM-1)  
hb2(u,u) = (Big_LGM_bed1(1,u)); 
for j =(u+1:TKM)
   Bb2(j,u) = -(Big_LGM_bed1(1,j)+Big_LGM_bed1(1,j-1));
   Cb2(j,u) = (hb2(j-1,u)*(Big_LGM_bed1(1,j)-(Hb2(j-1,u))))-((2*(west_x(1,j)-west_x(1,j-1))*(Tyyy(j,u)+Tyyy(j-1,u)/2))/(rho_ice*gk));
   hb2(j,u) =(-Bb2(j,u)+abs((Bb2(j,u)^2)-(4*Cb2(j,u)))^0.5)/2;
   Hb2 (j,u) = (hb2(j,u)-Big_LGM_bed1(1,j));
       
end
end

for LL = (1:TKM-1)
hb2(1:LL,LL+1) = Big_LGM_bed1(1,1:LL);  
end
for cvt = 1:edge
hb2(:,cvt) = hb2(:,edge);
end

fixH = zeros(TKM:edge);
for rr = 1:edge
    fixH(:,rr) = Hb;
end
Hb2(1:edge-1,1:edge-1) = 0;

figure
hold on
plot(west_x,hb2)
plot(west_x,Big_LGM_bed1)
plot(west_x,-140,'.')

% polyx  = [west_x fliplr(west_x)];         % x locations of my polygon
% for ab = 1:TKM;
% Polyys(1:1450,ab) = [hb2(1:TKM,ab).' fliplr(Big_LGM_bed1(1,1:TKM))]; % y locations of ice polygon
% end
% polyyb = [Sis_bed(1,1:TKM) -1200*ones(size(Sis_bed(1,1:TKM)))];        % y locations of my bed polygon (rock)
% 
%     figure(4); clf; hold on
%     pp= fill(polyx/1,Polyys,icecolor);  % 2003 ice polygon
%     p= fill(polyx/1,polyyb,rockcolor);     % bedrock polygon
%     p.LineStyle='none';pp.LineStyle='none';
%     text(6e5,-500,'Bedrock','fontsize',24)%,'rotation',22.5)
%     text(6e5,1000,'Ice Sheet','fontsize',24)%,'rotation',22.5)
%     title('Ice sheet and bed')
%     xlabel('km from terminus'); ylabel('meters a.s.l.')
%% small ice retreat 
hs2 = zeros(length(west_x),length(west_x)); % ice surface (m)
Hs2 = zeros(length(west_x),length(west_x)); % ice thickness (m)
Hs2(1:TKM,1) = 0;
Bs2 = zeros(length(west_x),length(west_x)); % b term of quadractic eqaution
Cs2 = zeros(length(west_x),length(west_x)); % c term of quadratic equation
Cs2 (1:TKM,1) = 0;
%  for rr = 1:TKM
%  hs2(1:TKM-SKM,rr) = Sis_bed(1,1:TKM-SKM);
%  Hs2(1:TKM-SKM,rr) = 0;
%  Bs2(1:TKM-SKM,rr) = 0;
%  Cs2(1:TKM-SKM,rr) = 0;
%  end

srt = (TKM-SKM)-edge; land = (268); smll = TKM-SKM;
for u=(TKM-SKM:TKM)  
% hs2(u,u) = (Small_LGM_bed1(1,u));
for j =(u+1:TKM)
   Bs2(j,u-srt) = -(Small_LGM_bed1(1,j)+Small_LGM_bed1(1,j-1));
   Cs2(j,u-srt) = (hs2(j-1,u-srt)*(Small_LGM_bed1(1,j)-(Hs2(j-1,u-srt))))-((2*(west_x(1,j)-west_x(1,j-1))*(Tyyys(j,u-srt)+Tyyys(j-1,u-srt)/2))/(rho_ice*gk));
   hs2(j,u-srt) =(-Bs2(j,u-srt)+abs((Bs2(j,u-srt)^2)-(4*Cs2(j,u-srt)))^0.5)/2;
   Hs2 (j,u-srt) = (hs2(j,u-srt)-Small_LGM_bed1(1,j));
end
end
 for c = 1:3
     Small_LGM_bed1(c,:) = Small_LGM_bed1(1,:);
 end
szhs2 = size(hs2(:,152:190)); 
nn = 3;
hs2s=reshape(repmat(hs2(:,152:190),nn,1),szhs2(1),nn*szhs2(2));
hs2a(:,152:268) = hs2s(:,1:117);

% for c = 1:26
%     x(23:c+22,(3*c-2:3*c)+17 ) = f(1:3,23:c+22).';%emerging land
% end
for LL = (1:land)
hs2a(smll:LL+smll-1,(3*LL-2:3*LL)+edge-1) = Small_LGM_bed1(1:3,smll:LL+smll-1).';  
end

for LL=smll:TKM
    hs2a(LL,1:edge-1) = hs2a(LL,edge);
end
for LL=1:TKM
    hs2a(1:smll,LL) = Small_LGM_bed1(1,1:smll).';
end
for LL = 1:TKM
hs2(LL:TKM,:) = hs2a(LL:TKM,1:TKM);
end
szHs2 = size(Hs2(:,152:190)); 
nn = 3;
Hs2s=reshape(repmat(Hs2(:,152:190),nn,1),szHs2(1),nn*szHs2(2));
Hs2(:,152:268) = Hs2s(:,1:117);

for u=(land:TKM)  
 hs2(u,u) = (Small_LGM_bed1(1,u));
for j =(u+1:TKM)
   Bs2(j,u) = -(Small_LGM_bed1(1,j)+Small_LGM_bed1(1,j-1));
   Cs2(j,u) = (hs2(j-1,u)*(Small_LGM_bed1(1,j)-(Hs2(j-1,u))))-((2*(west_x(1,j)-west_x(1,j-1))*(Tyyys(j,u)+Tyyys(j-1,u)/2))/(rho_ice*gk));
   hs2(j,u) =(-Bs2(j,u)+abs((Bs2(j,u)^2)-(4*Cs2(j,u)))^0.5)/2;
   Hs2 (j,u) = (hs2(j,u)-Small_LGM_bed1(1,j));
end
end
for LL = (1:TKM)
hs2(1:LL,LL) = Small_LGM_bed1(1,1:LL);  
end
% for cvt = 1:(TKM-SKM)
% hs2(:,cvt) = hs2(:,(TKM-SKM));
% end
% for ff = 1:(TKM-SKM)
% Hs2(:,ff) = Hs2(:,TKM-SKM);    
% end
figure
hold on
plot(west_x,hs2(:,:))
plot(west_x,Small_LGM_bed1)
%plot(west_x,-140,'.')
  
%% uplift if ice retreats (retreated 20 km more than today)
%make a matrix of the crust Hieght based on the location of the ice sheet
%so row one in big_LGM_bed, row 2 is one km of reateat uplift and so on. i
%will aslo need to replace the Big_LGM_bed function in the second loop
%with the new function and use u,j insead of 1,j

% water retreat 152-262  land reatreat 263-456
%Hb2 each coloum is ice sheet retreat one km
DIST_XX = fliplr(west_x);% distance from center of ice
DIST_XX(1,1:edge) = edge*1000;

Delta_hb = hb2;
for tt = END+20:TKM
Delta_hb(END:TKM,tt) =hb2(END:TKM,END+20);
end

Delta_Hb = Hb2;
for tt = END+20:TKM
Delta_Hb(END:TKM,tt) = Hb2(END:TKM,(END+20)); 
end
figure
hold on
plot(west_x,Delta_hb)
plot(west_x,Big_LGM_bed1)
plot(west_x,-140,'.')
%% uplift if small ice reatreat 
DIST_SS = fliplr(west_x);% distance from center of ice
DIST_SS(1,1:(TKM-SKM)) = (TKM-SKM)*1000;

Delta_hs = hs2;
for tt = END+20:TKM
Delta_hs(END:TKM,tt) =hs2(END:TKM,END+20);
end

Delta_Hs = Hs2;
for tt = END+20:TKM
Delta_Hs(END:TKM,tt) = Hs2(END:TKM,(END+20)); 
end
figure
hold on
plot(west_x,Delta_hs(:,:))
plot(west_x,Small_LGM_bed1)
plot(west_x,-140,'.')

%% find different in ice sheet ice and make that the value of w
% Ttt = zeros(1,725); Ttt(1,1:152) = 26;
% Ttt(1,153:725) = 26:26:(18850-(26*edge));
T = importdata('flex_timing.csv');
T = T.data;
Ta = T(:,1); %time step
Tb = T(:,2); %time ka
Tc = T(:,3);%distance km (big)
Td = T(:,4); %distance km (small)
vis = 4.67;viscos = 1;
tr2= ((4*pi*(vis*10^20))./(3200*9.81*((725*2)*1000)))/(365*24*60*60);
%(i changed the width of the ice shee from TKM*2 to 1146 and the viscosiy
%form 4.67 to 4) 7/29/2020
rbn = zeros(725:725); % c is the time step, the row of flex_real is distance
RBN = zeros(725);
for c = 1:145
rbn(c,1:TKM) = flex(146,1:TKM);
end
for i = 146:650 % edge to 304 uplfit from 15-2 ka
rbn(i,1:304) = (rbn(145,1:304))*2.718^(-(15106-Tb(i))/(tr2));
end
for i = 146:650% edge to TKM uplfit from 15-4.4 ka (makes it easy)
rbn(i,1:TKM) = (rbn(145,1:TKM))*2.718^(-(15106-Tb(i))/(tr2));
end

for ii = 0:TKM-3
for i = 650:TKM % front of fore deep to TKM sink 4.4-0 ka
rbn(i,(Tc(i)-150):end) = (rbn(649,(Tc(i)-150):end))-abs((rbn(649,(Tc(i)-150):end)) ...
    -((rbn(649,(Tc(i)-150):end))*2.718^(-(1976-Tb(i))/(tr2))));

% rbn(i,1:(Tc(650)-151)) = (rbn(649,1:(Tc(650)-151)))-abs((rbn(649,1:(Tc(650)-151))) ...
%     -((rbn(649,1:(Tc(650)-151)))*2.718^(-(1976-Tb(i))/(tr2))));

rbn(i,1+ii) = (rbn(649,1+ii)-abs((rbn(649,1+ii))) ...
    -((rbn(649,1+ii))*2.718^(-(1976-Tb(i))/(tr2))));


rbn(i,1:(Tc(650)-407)) = (rbn(649,1:(Tc(650)-407)))+abs((rbn(649,1:(Tc(650)-407))) ...
    -((rbn(649,1:(Tc(650)-407)))*2.718^(-(1976-Tb(i))/(tr2))));
end
end

for c = 1:TKM
RBN(c,:) = rbn(c,:)+Sis_bed(1,1:TKM);
end

F_bed = RBN;
Time = (Tb); TIME = Time';

% find different in ice sheet ice and make that the value of w SMALL
% rbns = zeros(725:725);
for c = 1:145
rbns(c,1:TKM) = flexs(146,1:TKM);
end
for i = 146:650 % edge to 304 uplfit from 15-2 ka
rbns(i,1:304) = (rbns(145,1:304))*2.718^(-(15106-Tb(i))/(tr2));
end
for i = 146:650% edge to TKM uplfit from 15-4.4 ka (makes it easy)
rbns(i,1:TKM) = (rbns(145,1:TKM))*2.718^(-(15106-Tb(i))/(tr2));
end

for ii = 0:TKM-3
for i = 650:TKM % front of fore deep to TKM sink 4.4-0 ka
rbns(i,(Td(i)-150):end) = (rbns(649,(Td(i)-150):end))-abs((rbns(649,(Td(i)-150):end)) ...
    -((rbns(649,(Td(i)-150):end))*2.718^(-(1976-Tb(i))/(tr2))));

% rbn(i,1:(Tc(650)-151)) = (rbn(649,1:(Tc(650)-151)))-abs((rbn(649,1:(Tc(650)-151))) ...
%     -((rbn(649,1:(Tc(650)-151)))*2.718^(-(1976-Tb(i))/(tr2))));

rbns(i,1+ii) = (rbns(649,1+ii)-abs((rbns(649,1+ii))) ...
    -((rbns(649,1+ii))*2.718^(-(1976-Tb(i))/(tr2))));


rbns(i,1:(Td(650)-407)) = (rbns(649,1:(Td(650)-407)))+abs((rbns(649,1:(Td(650)-407))) ...
    -((rbns(649,1:(Td(650)-407)))*2.718^(-(1976-Tb(i))/(tr2))));
end
end

for c = 1:TKM
RBNs(c,:) = rbns(c,:)+Sis_bed(1,1:TKM);
end

F_beds = RBNs;
Time = (Tb); TIME = Time';

%%

hb3 = zeros(length(west_x),length(west_x)); % ice surface (m)
Hb3 = zeros(length(west_x),length(west_x)); % ice thickness (m)
Hb3(1:TKM,1) = 0;
Bb3 = zeros(length(west_x),length(west_x)); % b term of quadractic eqaution
Cb3 = zeros(length(west_x),length(west_x)); % c term of quadratic equation
Cb3 (1:TKM,1) = 0;


for uuu=(1:724)  
hb3(uuu,uuu) = (F_bed(uuu,uuu)); 
for j =(uuu+1:TKM)
   Bb3(j,uuu) = -(F_bed(j-1,j)+F_bed(j-1,j-1));
   Cb3(j,uuu) = (hb3(j-1,uuu)*(F_bed(j-1,j)-(Hb3(j-1,uuu))))-((2*(west_x(1,j)-west_x(1,j-1))*(Tyyy(j,uuu)+Tyyy(j-1,uuu)/2))/(rho_ice*gk));
   hb3(j,uuu) =(-Bb3(j,uuu)+abs((Bb3(j,uuu)^2)-(4*Cb3(j,uuu)))^0.5)/2;
   Hb3 (j,uuu) = (hb3(j,uuu)-F_bed(j-1,j));
       
end
end

for LL = (1:724)
hb3(1:LL,LL+1) = F_bed(LL,1:LL);  
end
fixh = zeros(TKM:edge);
for rr = 1:edge
    fixh(:,rr) = hb;
end
%  hb3(:,1:edge) = fixh;
fixH = zeros(TKM:edge);
for rr = 1:edge
    fixH(:,rr) = Hb;
end
Hb3(:,1:edge) = fixH;

Delta_hb3 = hb3;
for tt = END:END+54
Delta_hb3(END:tt,tt) = F_bed(tt,END:tt);   
Delta_hb3(tt:TKM,tt) = hb3(tt:TKM,tt); 
end
for tt = END+54:(TKM-54)
    Delta_hb3(END:END+54,tt) = F_bed(tt,END:END+54);
    Delta_hb3(END+54:TKM,tt) = hb3(END+54:TKM,END+54);
end

 hhb3 = hb3;
 hhb3(END:TKM,TKM-54:TKM) = fliplr(hb3(END:TKM,END:END+54));

for tt = 671:725
Delta_hb3(END+(725-tt):TKM,tt) = hhb3(END+(725-tt):TKM,tt); 
Delta_hb3(END:END+(725-tt),tt) = F_bed(tt,END:END+(725-tt));
end

Delta_Hb3 = Hb3;
for tt = END:TKM
Delta_Hb3(END:TKM,tt) = Hb3(END:TKM,END); 
end
Delta_hb3(1:edge,edge)=F_bed(edge,1:edge);
Delta_hb3(1:END,TKM)=F_bed(TKM,1:END);

SEAb = zeros(1,725)-126.82;
polyx  = [west_x fliplr(west_x)];         % x locations of my polygon
bbPolyys(1:1450,1) = [Delta_hb3(1:TKM,edge).' fliplr(F_bed(edge,1:TKM))]; % y locations of ice polygon
bbPolyyb(1:1450,1) = [F_bed(edge,1:TKM) -1200*ones(size(F_bed(edge,1:TKM)))];        % y locations of my bed polygon (rock)
bbpolywater(1,1:2*edge) = [1:edge fliplr(1:edge)];
bbPolyyw(1,1:2*edge) = [SEAb(1,1:edge) min(fliplr(F_beds(edge,1:edge)),SEAb(1,edge))];
figure
hold on    
fill(polyx,bbPolyys(:,1),icecolor); 
hold on; fill(polyx,bbPolyyb(:,1),rockcolor);hold on;
p=fill(bbpolywater(1,:)*1000,bbPolyyw(1,:),b);
p.LineStyle='none';
text(5.5e5,-500,'Bedrock','fontsize',12);
text(5.5e5,1000,'Ice Sheet','fontsize',12);
text(5e4,-500,'Ocean','fontsize',8);
title('Maximum LGM Ice Sheet', 'FontSize', 12);
xlabel('Distance (km)');ylabel('Elevation (m)');
ax = gca;ax.XLim = [0 TKM*1000];ax.YLim = [-1200 3500];
xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/1000) 
%% small
hs3 = zeros(length(west_x),length(west_x)); % ice surface (m)
Hs3 = zeros(length(west_x),length(west_x)); % ice thickness (m)
Hs3(1:TKM,1) = 0;
Bs3 = zeros(length(west_x),length(west_x)); % b term of quadractic eqaution
Cs3 = zeros(length(west_x),length(west_x)); % c term of quadratic equation
Cs3 (1:TKM,1) = 0;
% for rr = 1:TKM
%  hs2(1:TKM-SKM,rr) = Sis_bed(1,1:TKM-SKM);
%  Hs2(1:TKM-SKM,rr) = 0;
%  Bs2(1:TKM-SKM,rr) = 0;
%  Cs2(1:TKM-SKM,rr) = 0;
%  end

for uuu=(TKM-SKM:TKM)  
%  hs3(uuu,uuu) = (F_beds(uuu,uuu)); 
for j =(uuu+1:TKM)
   Bs3(j,uuu-srt) = -(F_beds(j-1,j)+F_beds(j-1,j-1));
   Cs3(j,uuu-srt) = (hs3(j-1,uuu-srt)*(F_beds(j-1,j)-(Hs3(j-1,uuu-srt))))-((2*(west_x(1,j)-west_x(1,j-1))*(Tyyys(j,uuu-srt)+Tyyys(j-1,uuu-srt)/2))/(rho_ice*gk));
   hs3(j,uuu-srt) =(-Bs3(j,uuu-srt)+abs((Bs3(j,uuu-srt)^2)-(4*Cs3(j,uuu-srt)))^0.5)/2;
   Hs3 (j,uuu-srt) = (hs3(j,uuu-srt)-F_beds(j-1,j));
       
end
end
szhs3 = size(hs3(:,152:190)); 
nn = 3;
hs3s=reshape(repmat(hs3(:,152:190),nn,1),szhs3(1),nn*szhs3(2));
hs3a(:,152:268) = hs3s(:,1:117);
F_BEDDS = F_beds.';
bed3 = size(F_BEDDS(:,152:190)); 
nn = 3;
bed3s=reshape(repmat(F_BEDDS(:,152:190),nn,1),bed3(1),nn*bed3(2));
bed3a(1:725,152:268) = bed3s(1:725,1:117);
for LL = (1:38)
hs3a(smll+1:LL+smll,(3*LL+edge:3*LL+2+edge)) = bed3a(smll+1:LL+smll,(3*LL+edge:3*LL+2+edge));  
end
for LL=smll:TKM
    hs3a(LL,1:edge-1) = hs3a(LL,edge);
end
    hs3a(1:smll,1:TKM) = F_BEDDS(1:smll,1:TKM);
for LL = 1:TKM
hs3(LL:TKM,:) = hs3a(LL:TKM,1:TKM);
end
szHs3 = size(Hs3(:,152:190)); 
nn = 3;
Hs3s=reshape(repmat(Hs3(:,152:190),nn,1),szHs3(1),nn*szHs3(2));
Hs3(:,152:268) = Hs3s(:,1:117);
for uuu=(land:TKM)  
   hs3(uuu,uuu) = (F_beds(uuu,uuu));
for j =(uuu+1:TKM)
   Bs3(j,uuu) = -(F_beds(j-1,j)+F_beds(j-1,j-1));
   Cs3(j,uuu) = (hs3(j-1,uuu)*(F_beds(j-1,j)-(Hs3(j-1,uuu))))-((2*(west_x(1,j)-west_x(1,j-1))*(Tyyys(j,uuu)+Tyyys(j-1,uuu)/2))/(rho_ice*gk));
   hs3(j,uuu) =(-Bs3(j,uuu)+abs((Bs3(j,uuu)^2)-(4*Cs3(j,uuu)))^0.5)/2;
   Hs3 (j,uuu) = (hs3(j,uuu)-F_beds(j-1,j));
end
end
for LL = (land+1:TKM)
hs3(231:268,LL) = F_beds(LL,231:268);
hs3(269:LL,LL) = F_beds(LL,269:LL);
end

Delta_hs3 = hs3;
for tt = END:END+54
Delta_hs3(END:tt,tt) = F_beds(tt,END:tt);   
Delta_hs3(tt:TKM,tt) = hs3(tt:TKM,tt); 
end
for tt = END+54:(TKM-54)
    Delta_hs3(END:END+54,tt) = F_beds(tt,END:END+54);
    Delta_hs3(END+54:TKM,tt) = hs3(END+54:TKM,END+54);
end

 hhs3 = hs3;
 hhs3(END:TKM,TKM-54:TKM) = fliplr(hs3(END:TKM,END:END+54));

for tt = 671:725
Delta_hs3(END+(725-tt):TKM,tt) = hhs3(END+(725-tt):TKM,tt); 
Delta_hs3(END:END+(725-tt),tt) = F_beds(tt,END:END+(725-tt));
end

Delta_Hs3 = Hs3;
for tt = END:TKM
Delta_Hs3(END:TKM,tt) = Hs3(END:TKM,END); 
end

SEAs = zeros(1,725)-126.82;
polyx  = [west_x fliplr(west_x)];         % x locations of my polygon
ssPolyys(1:1450,1) = [Delta_hs3(1:TKM,152).' fliplr(F_beds(152,1:TKM))]; % y locations of ice polygon
ssPolyyb(1:1450,1) = [F_beds(230,1:TKM) -1200*ones(size(F_beds(230,1:TKM)))];        % y locations of my bed polygon (rock)
sspolywater(1,1:2*230) = [1:230 fliplr(1:230)];
ssPolyyw(1,1:2*230) = [SEAs(1,1:230) min(fliplr(F_beds(230,1:230)),SEAs(1,152))];
figure
hold on    
fill(polyx,ssPolyys(:,1),icecolor); 
hold on; fill(polyx,ssPolyyb(:,1),rockcolor);hold on;
p=fill(sspolywater(1,:)*1000,ssPolyyw(1,:),b);
p.LineStyle='none';
text(5.5e5,-500,'Bedrock','fontsize',12);
text(5.5e5,1000,'Ice Sheet','fontsize',12);
text(5e4,-500,'Ocean','fontsize',8);
title('Minimium LGM Ice Sheet', 'FontSize', 12);
xlabel('Distance (km)');ylabel('Elevation (m)');
ax = gca;ax.XLim = [0 TKM*1000];ax.YLim = [-1200 3500];
xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/1000) 
%% global sea level
%TT = 14500-T;
ML= importdata('fix_sea_lvl.csv');
for qc = 1:TKM
    %y = -8E-07x2 + 0.0046x - 3.8686
    %y = 2E-14x4 - 4E-10x3 + 3E-06x2 - 0.0064x + 1.324
 %sea(1,qc) = -(((2*10^(-14))*((TIME(qc))^4))-((3.95*10^(-10))*((TIME(qc))^3))+((3*10^(-6))*((TIME(qc))^2))-(.0064*(TIME(qc)))+1.324);
%sea(1,qq) = -(.58*2.718^(.0005*(14500-T(1,qq))));
end
sea =(ML.data(:,2)).';
sea =fliplr(sea);
for rt = 1:TKM
sea(rt,1:TKM) = sea(1,1:TKM);
end
figure
hold on
sc1=xline(11.9,'-.r','linewidth',2);sc2=xline(8.5,'-.g','linewidth',2);
plot((TIME(1,1:TKM)/1000),sea(1,1:TKM),'linewidth',3)
xlabel('Time (KY)'); ylabel('Global Sea Level (m)');title('Global Sea Level since LGM')
legend([sc1 sc2],'Time of deglaciation for Sisimiut (11.9 ka)', ...
    'Time of deglaciation for Kangerlussuaq (8.5 ka)','location','northwest','fontsize',12);
ax = gca; % current axes 
ax.XLim = [0 TIME(1,1)/1000];ax.YLim = [-140 0];set(gca,'xdir','reverse')
%% Plots
polyx  = [west_x fliplr(west_x)];         % x locations of my polygon
for ab = 1:TKM  
bPolyy(1:1450,ab) = [Delta_hb3(1:TKM,ab).' fliplr(F_bed(ab,1:TKM))]; % y locations of ice polygon
bPolyyb(1:1450,ab) = [F_bed(ab,1:TKM) -1200*ones(size(F_bed(1,1:TKM)))];        % y locations of my bed polygon (rock)
end
Big_names = {'LGM ','12 Ka  ','11.6 Ka  ','10.4 Ka ','9.1 Ka  ','8.1 Ka  ','7.3 ka  ','4-2 Ka ','Current '};
Big_steps = [edge 268 300 328 375 402 427 488 TKM];
% Big_names = {'12 Ka  ','11.6 Ka  ','10.4 Ka '};
% Big_steps = [ 268 300 328];
bpolywater=zeros(1,1450);
bPolyyw =zeros(1,1450);
for bb = 1:9
%   for bb = 1:3
bpolywater(bb,1:2*Big_steps(bb)) = [1:Big_steps(bb) fliplr(1:Big_steps(bb))];
bPolyyw(bb,1:2*Big_steps(bb)) = [sea(1:Big_steps(bb),Big_steps(bb)).' min(fliplr(F_bed(Big_steps(bb),1:Big_steps(bb))),sea(1:Big_steps(bb),Big_steps(bb)).')];
end

figure
for bb = 1:9
%   for bb = 1:3
hold on    
subplot(3,3,bb);
fill(polyx,bPolyy(:,Big_steps(bb)),icecolor); 
hold on; fill(polyx,bPolyyb(:,Big_steps(bb)),rockcolor);hold on;
p=fill(bpolywater(bb,:)*1000,bPolyyw(bb,:),b);
p.LineStyle='none';
text(5.5e5,-500,'Bedrock','fontsize',12);
text(5.5e5,1000,'Ice Sheet','fontsize',12);
text(5e4,-500,'Ocean','fontsize',8);
title(Big_names{bb}, 'FontSize', 12);
xlabel('Distance (km)');ylabel('Elevation (m)');
ax = gca;ax.XLim = [0 TKM*1000];ax.YLim = [-1200 3500];
xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/1000) 
end
qq = annotation('textbox', [.33,.9, 0.5,.1], 'String', "Maximum Ice Sheet Extent from LGM to Today ",'fontsize',20);
qq.LineStyle ='none';
%% Plots small
polyx  = [west_x fliplr(west_x)];         % x locations of my polygon
for ab = 1:TKM  
sPolyys(1:1450,ab) = [Delta_hs3(1:TKM,ab).' fliplr(F_beds(ab,1:TKM))]; % y locations of ice polygon
sPolyyb(1:1450,ab) = [F_beds(ab,1:TKM) -1200*ones(size(F_beds(1,1:TKM)))];        % y locations of my bed polygon (rock)
end
Small_names = {'LGM ','12 Ka  ','11.6 Ka  ','10.4 Ka ','9.1 Ka  ','8.1 Ka  ','7.3 ka  ','4-2 Ka ','Current '};
Small_steps = [152 268 300 328 375 402 427 488 TKM];
Small_stepsx = [171 268 300 328 375 402 427 488 TKM];
% Small_names = {'12 Ka  ','11.6 Ka  ','10.4 Ka '};
% Small_steps = [268 300 328];
% Small_stepsx = [268 300 328];


polywater=zeros(1,1450);
sPolyyw =zeros(1,1450);
for bb = 1:9
% for bb = 1:3
polywater(bb,1:2*Small_stepsx(bb)) = [1:Small_stepsx(bb) fliplr(1:Small_stepsx(bb))];
sPolyyw(bb,1:2*Small_stepsx(bb)) = [sea(1:Small_stepsx(bb),Small_stepsx(bb)).' min(fliplr(F_beds(Small_stepsx(bb),1:Small_stepsx(bb))),sea(1:Small_stepsx(bb),Small_stepsx(bb)).')];
end

figure
for bb = 1:9
% for bb = 1:3
hold on    
subplot(3,3,bb);

fill(polyx,sPolyys(:,Small_steps(bb)),icecolor); 
hold on
p=fill(polywater(bb,:)*1000,sPolyyw(bb,:),b);
p.LineStyle='none';

hold on; fill(polyx,sPolyyb(:,Small_steps(bb)),rockcolor);hold on;


text(5.5e5,-500,'Bedrock','fontsize',12);
text(5.5e5,1000,'Ice Sheet','fontsize',12);
text(5e4,-500,'Ocean','fontsize',8);
title(Small_names{bb}, 'FontSize', 12);
xlabel('Distance (km)');ylabel('Elevation (m)');
ax = gca;ax.XLim = [0 TKM*1000];ax.YLim = [-1200 3500];
xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/1000) 
end
qq = annotation('textbox', [.33,.9, 0.5,.1], 'String', "Minimum Ice Sheet Extent from LGM to Today ",'fontsize',20)
qq.LineStyle ='none';
%% Retreating ice (change shear stress in big ice and hit run) 
hbf = zeros(length(west_x),length(west_x)); % ice surface (m)
Hbf = zeros(length(west_x),length(west_x)); % ice thickness (m)
Hbf(1:TKM,1) = 0;
Bbf = zeros(length(west_x),length(west_x)); % b term of quadractic eqaution
Cbf = zeros(length(west_x),length(west_x)); % c term of quadratic equation
Cbf (1:TKM,1) = 0;
end_bed_B = F_bed(1,:);
for u=(1:TKM-1)  
hbf(u,u) = (end_bed_B(1,u)); 
for j =(u+1:TKM)
   Bbf(j,u) = -(end_bed_B(1,j)+end_bed_B(1,j-1));
   Cbf(j,u) = (hbf(j-1,u)*(end_bed_B(1,j)-(Hbf(j-1,u))))-((2*(west_x(1,j)-west_x(1,j-1))*(Tyyy(j,u)+Tyyy(j-1,u)/2))/(rho_ice*gk));
   hbf(j,u) =(-Bbf(j,u)+abs((Bbf(j,u)^2)-(4*Cbf(j,u)))^0.5)/2;
   Hbf (j,u) = (hbf(j,u)-end_bed_B(1,j));
       
end
end

for LL = (1:TKM-1)
hbf(1:LL,LL+1) = end_bed_B(1,1:LL);  
end
for cvt = 1:edge
hbf(:,cvt) = hbf(:,edge);
end

fixfB = zeros(TKM:edge);
for rr = 1:edge
    fixfB(:,rr) = Hb;
end
Hbf(1:edge,1:edge) = 0;
BED = zeros(725);
for ff = 1:TKM
BED(ff,:) = F_bed(1,1:TKM);
end
figure
hold on
c = [.23 .71 .69];
rectangle('Position',[0 -1200 7250000 1092],'FaceColor',c,'EdgeColor','n');hold on;
plot(west_x,F_bed(1,:))
plot(west_x,hbf(:,[edge 268 303 324 374 410 427 435 488]),'k','linewidth',1.5);
polyx  = [west_x fliplr(west_x)];         % x locations of my polygon
for ab = 1:TKM  
sPolyys(1:1450,ab) = [hbf(1:TKM,ab).' fliplr(BED(ab,1:TKM))]; % y locations of ice polygon
end

ccPolyyb(1:1450,1) = [F_bed(1,1:TKM) -1200*ones(size(F_bed(1,1:TKM)))];        % y locations of my bed polygon (rock)
ccpolywater(bb,1:2*END) = [1:END fliplr(1:END)];
ccPolyyw(bb,1:2*END) = [SEAb(1,1:END) min(fliplr(F_bed(1,1:END)),SEAb(1,1:END))];  
hold on; fill(polyx,ccPolyyb(:,1),rockcolor);hold on;
big_steps = [edge 268 300 328 375 402 427 435 488 ];

% p=fill(ccpolywater(1,:)*1000,ccPolyyw(1,:),b);
% p.LineStyle='none';
for bb = 1:9
fill(polyx,sPolyys(:,big_steps(bb)),icecolor);
end
hold on;
plot(west_x,Sis_bed(1,1:TKM),'k','linewidth', 1.5);
% text(5.5e5,-500,'Bedrock','fontsize',12);
% text(5e4,-500,'Ocean','fontsize',8);
% text(2e5,900,'LGM','fontsize',12,'rotation',40);text(2.6e5,650,'12 Ka','fontsize',12,'rotation',45);
% text(2.75e5,570,'11.6 Ka','fontsize',12,'rotation',45);text(3.3e5,800,'10.4 Ka','fontsize',12,'rotation',45);
% text(3.72e5,800,'9.1 Ka','fontsize',12,'rotation',45);text(3.9e5,582,'8.1 Ka','fontsize',12,'rotation',65);
% text(4.25e5,582,'7.3 Ka','fontsize',12,'rotation',55);text(4.45e5,582,'Current','fontsize',12,'rotation',55);
% text(4.99e5,582,'4-2 Ka','fontsize',12,'rotation',55);text(5.5e5,1000,'Ice Sheet','fontsize',12);
% title('Ice Sheet Retreat from Maximum LGM to Today', 'FontSize', 12);
% xlabel('Distance (km)');ylabel('Elevation (m)');
ax = gca;ax.XLim = [0 TKM*1000];ax.YLim = [-1200 3500];
xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/1000) 

%% Retreating small ice (change shear stress in small ice and hit run) 
hsf = zeros(length(west_x),length(west_x)); % ice surface (m)
Hsf = zeros(length(west_x),length(west_x)); % ice thickness (m)
Hsf(1:TKM,1) = 0;
Bsf = zeros(length(west_x),length(west_x)); % b term of quadractic eqaution
Csf = zeros(length(west_x),length(west_x)); % c term of quadratic equation
Csf (1:TKM,1) = 0;
end_bed_S = F_beds(1,:);
for u=(1:TKM-1)  
hsf(u,u) = (end_bed_S(1,u)); 
for j =(u+1:TKM)
   Bsf(j,u) = -(end_bed_S(1,j)+end_bed_S(1,j-1));
   Csf(j,u) = (hsf(j-1,u)*(end_bed_S(1,j)-(Hsf(j-1,u))))-((2*(west_x(1,j)-west_x(1,j-1))*(Tyyys(j,u)+Tyyys(j-1,u)/2))/(rho_ice*gk));
   hsf(j,u) =(-Bsf(j,u)+abs((Bsf(j,u)^2)-(4*Csf(j,u)))^0.5)/2;
   Hsf (j,u) = (hsf(j,u)-end_bed_S(1,j));
       
end
end

for LL = (1:TKM-1)
hsf(1:LL,LL+1) = end_bed_S(1,1:LL);  
end
for cvt = 1:edge
hsf(:,cvt) = hsf(:,edge);
end

fixfS = zeros(TKM:edge);
for rr = 1:edge
    fixfS(:,rr) = Hs;
end
Hsf(1:edge,1:edge) = 0;
BEDs = zeros(725);
for ff = 1:TKM
BEDs(ff,:) = F_beds(1,1:TKM);
end
figure
hold on
rectangle('Position',[0 -1200 7250000 1092],'FaceColor',c,'EdgeColor','n');hold on;
plot(west_x,F_beds(1,:))
plot(west_x,hsf(:,[230 268 303 324 374 410 427 435 488]),'k','linewidth',1.5);
polyx  = [west_x fliplr(west_x)];         % x locations of my polygon
for ab = 1:TKM  
csPolyys(1:1450,ab) = [hsf(1:TKM,ab).' fliplr(BEDs(ab,1:TKM))]; % y locations of ice polygon
end
csPolyyb(1:1450,1) = [F_beds(1,1:TKM) -1200*ones(size(F_beds(1,1:TKM)))];        % y locations of my bed polygon (rock)
% cspolywater(bb,1:2*END) = [1:END fliplr(1:END)];
% csPolyyw(bb,1:2*END) = [SEA(1,1:END) min(fliplr(Sis_bed(1,1:END)),SEA(1,END))];  
hold on; fill(polyx,csPolyyb(:,1),rockcolor);hold on;
small_steps = [230 268 300 328 375 402 427 435 488 ];
for bb = 1:9
fill(polyx,csPolyys(:,small_steps(bb)),icecolor);
end
hold on;
plot(west_x,Sis_bed(1,1:TKM),'k','linewidth', 1.5);
% p=fill(ccpolywater(1,:)*1000,ccPolyyw(1,:),b);
% p.LineStyle='none';
% text(5.5e5,-500,'Bedrock','fontsize',12);
% text(5e4,-500,'Ocean','fontsize',8);
% text(2.5e5,750,'LGM','fontsize',12,'rotation',40);text(2.62e5,650,'12 Ka','fontsize',12,'rotation',45);
% text(2.75e5,570,'11.6 Ka','fontsize',12,'rotation',45);text(3.3e5,800,'10.4 Ka','fontsize',12,'rotation',45);
% text(3.72e5,800,'9.1 Ka','fontsize',12,'rotation',45);text(3.9e5,582,'8.1 Ka','fontsize',12,'rotation',65);
% text(4.25e5,582,'7.3 Ka','fontsize',12,'rotation',55);text(4.45e5,582,'Current','fontsize',12,'rotation',55);
% text(4.99e5,582,'4-2 Ka','fontsize',12,'rotation',55);text(5.5e5,1000,'Ice Sheet','fontsize',12);
% title('Ice Sheet Retreat from Minimum LGM to Today', 'FontSize', 12);
% xlabel('Distance (km)');ylabel('Elevation (m)');
ax = gca;ax.XLim = [0 TKM*1000];ax.YLim = [-1200 3500];
xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/1000) 

%% sea level curves

loc =267; % study 268  KANG 400 SIS 275
F_bed1 = F_bed;

elv = F_bed1(edge:TKM,loc).';
uplift = zeros(1,725);
for xx = 1:TKM
uplift(viscos,xx) =F_bed1(TKM,loc)-F_bed1(xx,loc)+sea(1,xx);
end 

F_beds1 = F_beds;
uplifts = zeros(1,725);

for xx = 1:TKM
uplifts(viscos,xx) =F_beds1(TKM,loc)-F_beds1(xx,loc)+sea(1,xx);
end 

%data
Bennike_2011y = [129.9	94.9	94.9	98.4	98.4	98.4	98.4	39.8	18.4];
Bennike_2011x = [10193	10280	10421	9111.5	9784.5	9476.5	10844	8423	5075.5];
Long_2008by = [15.40	10.16	2.70	0.27	0.27	-1.88	-1.88	-1.98	-1.98]; 
Long_2008bx = [5743.5	4633	3488.5	2248	1181	1264	1393	2221.5	2153.5]; 

%plot
figure
hold on
rectangle('Position',[12000 -19 4000 500],'FaceColor',icecolor,'EdgeColor','n');hold on;
plot(TIME(1,1:TKM),(1*uplifts(1,:)),'r','linewidth',2);
plot(TIME(1,1:TKM),(uplift(1,:)),'b','linewidth',2);hold on;

err_ben = [78	150	166	349.5	384.5	72.5	213	101	207.5];%
errorbar(Bennike_2011x,Bennike_2011y,err_ben,'k','horizontal', 'linestyle','n', 'LineWidth', .7)
err_long = [137.5	184	89.5	95	100	81	93	101.5	151.5];
errorbar(Long_2008bx,Long_2008by,err_long,'r','horizontal', 'linestyle','n','LineWidth', .7)

set(gca,'TickDir','out','XMinorTick','on','YMinorTick','on');

scatter(Long_2008bx([1,2,4,6]),Long_2008by([1,2,4,6]),25,'<','filled','r','MarkerEdgeColor','k');
scatter(Long_2008bx([3,5,7,8,9]),Long_2008by([3,5,7,8,9]),25,'>','filled','r','MarkerEdgeColor','k');
scatter(Bennike_2011x([1:6,8,9]),Bennike_2011y([1:6,8,9]),25,'<','filled','k','MarkerEdgeColor','k');
scatter(Bennike_2011x([7]),Bennike_2011y([7]),25,'^','filled','k','MarkerEdgeColor','k');

xlabel('Age (cal yr BP)','Fontsize',20);ylabel('Elevation (m asl)','Fontsize',20)
% legend([s2 s1],'Bennike et al 2011 ^1^4C dates', ...
%     'Long et al 2008b ^1^4C dates','location','northeast','fontsize',12)

ax = gca; % current axes 
ax.XLim = [0 TIME(1,110)];ax.YLim = [-20 400];set(gca,'xdir','reverse')
pbaspect([1 1 1])

%%
% kang
loc =391; % study 268  KANG 400 SIS 275
% elv = F_bed(edge:TKM,loc).';
uplift2 = zeros(1,725);%
% for xx = 1:val
% uplift2(viscos,xx) =F_bed1(TKM,loc)-F_bed1(val,loc)+sea(1,xx);
% end
for xx = 1:TKM
uplift2(viscos,xx) =F_bed1(TKM,loc)-F_bed1(xx,loc)+sea(1,xx);
end 
figure;
hold on
zz1 =subplot(1,2,1); hold on
rectangle('Position',[8500 -20 7500 600],'FaceColor',icecolor);hold on;
plot(TIME(1,1:TKM),(uplift2(1,:)),'linewidth',2);
text(15400,200,'Location covered','fontsize',12);
text(15400,185,'by Ice Sheet','fontsize',12);
text(5000,200,'Location is','fontsize',12);
text(5000,185,'ice free','fontsize',12);
Storm_2012y =[12	12	25	26	13	46	2.89	1.51	1.58]; Storm_2012x = [6573.5	8056	6762	6541	6792.5	7949.5	198	268	563.5];
Ten_Brink_1975y = [34	15	10	14	4	24	10	1	0	31	95	20	37	55	47	57	37	20	19	105	115];
Ten_Brink_1975x = [6150	5070	5845	7260	5615	6045	7025	4335	119	6505	8250	6510 ...
    7220	7730	8460	6880	6333	6141	6060	8611	8670];
k1=scatter(Storm_2012x,Storm_2012y,'filled','g');k2=scatter(Ten_Brink_1975x,Ten_Brink_1975y,'filled','m');
k3=yline(74,'-.','linewidth',3); 
title('RSL Curve Near Kangerlussuaq Max ice extent');xlabel('Time(yrs)');ylabel('Elevation (m)')
legend([k1 k2 k3],'Storm et al 2012 ^1^4C dates', ...
     'Ten Brink et al 1975 ^1^4C dates',...
     'Local Marine Limit','location','northeast','fontsize',12)

ax = gca;  % current axes 
ax.XLim = [0 TIME(1,110)];ax.YLim = [-20 250];set(gca,'xdir','reverse')% current axes 
 %small SL curves sis

% kang
loc =409; % study 267  KANG 400 SIS 275
elvs = F_beds1((edge):TKM,loc).';
uplifts2 = zeros(1,725);
% for xx = 1:vals
% uplifts2(viscos,xx) =F_beds1(TKM,loc)-F_beds1(vals,loc)+sea(1,xx);
% end
for xx = 1:TKM
uplifts2(viscos,xx) =F_beds1(TKM,loc)-F_beds1(xx,loc)+sea(1,xx);
end 
zz4 =subplot(1,2,2); hold on
hold on
rectangle('Position',[8500 -20 7500 600],'FaceColor',icecolor);hold on; %400
plot(TIME(1,1:TKM),(1*uplifts2(1,:)),'linewidth',2);
text(15400,160,'Location covered','fontsize',12);
text(15400,145,'by Ice Sheet','fontsize',12);
text(5000,160,'Location is','fontsize',12);
text(5000,145,'ice free','fontsize',12);
k1=scatter(Storm_2012x,Storm_2012y,'filled','g');k2=scatter(Ten_Brink_1975x,Ten_Brink_1975y,'filled','m');
k3=yline(54,'-.','linewidth',3); 
title('RSL Curve Near Kangerlussuaq Minimum ice extent');xlabel('Time(yrs)');ylabel('Elevation (m)')
legend([k1 k2 k3],'Storm et al 2012 ^1^4C dates', ...
     'Ten Brink et al 1975 ^1^4C dates',...
     'Local Marine Limit','location','northeast','fontsize',12)

ax = gca; % current axes 
ax.XLim = [0 TIME(1,110)];ax.YLim = [-20 250];set(gca,'xdir','reverse')

% figure
% plot(west_x/1e3,Sis_bed(1, 1:TKM))


%  set(zz1,'position',[.1 .5 .36 .31])
%  set(zz2,'position',[.1 .1 .36 .31])
%  set(zz3,'position',[.5 .5 .36 .31])
%  set(zz4,'position',[.5 .1 .36 .31])

%%
videofilename = 'big_ice_flextest1.avi';
vid = VideoWriter(videofilename,'Uncompressed AVI');
vid.FrameRate = 1; % frames per second
open(vid)
polyx  = [west_x fliplr(west_x)];         % x locations of my polygon
for ab = 1:TKM  
fbPolyy(1:1450,ab) = [Delta_hb3(1:TKM,ab).' fliplr(F_bed(ab,1:TKM))]; % y locations of ice polygon
fbPolyyb(1:1450,ab) = [F_bed(ab,1:TKM) -1200*ones(size(F_bed(1,1:TKM)))];        % y locations of my bed polygon (rock)
end
fBig_names = {'15 ka ','12 Ka  ','11.6 Ka  ','10.4 Ka ','9.1 Ka  ','8.1 Ka  ','7.3 ka  ','4 Ka','2 Ka ','Current '};
fBig_steps = [edge 268 300 328 375 402 427 488 662 TKM];
fTime_steps = [149 264 280 326 376 414 445 572 649 725];
fbpolywater=zeros(1,1450);
fbPolyyw =zeros(1,1450);
for bb = 1:10
%   for bb = 1:3
fbpolywater(bb,1:2*fBig_steps(bb)) = [1:fBig_steps(bb) fliplr(1:fBig_steps(bb))];
fbPolyyw(bb,1:2*fBig_steps(bb)) = [sea(1:fBig_steps(bb),fBig_steps(bb)).' min(fliplr(F_bed(fBig_steps(bb),1:fBig_steps(bb))),sea(1:fBig_steps(bb),fBig_steps(bb)).')];
end
hh=figure; clf; hold on
for iii=1:10
   clf; hold on; 
        p=fill(fbpolywater((iii),:)*1000,fbPolyyw((iii),:),b);
        hhh(iii) =fill(polyx,fbPolyy(:,fBig_steps(iii)),icecolor); 
        hold on; fill(polyx,fbPolyyb(:,fBig_steps(iii)),rockcolor);hold on;
        p.LineStyle='none';
        title([sprintf('%d',(round(TIME(fTime_steps(iii)),-2))),' years before present for Maximum LGM Extent'])
        text(5.5e5,-500,'Bedrock','fontsize',12);
        text(5.5e5,1000,'Ice Sheet','fontsize',12);
        text(5e4,-500,'Ocean','fontsize',8);
        xlabel('Distance (km)');ylabel('Elevation (m)');
        axis([1 (7.25*10^5) -1000 4000]);
        xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/1000) 
          frame = getframe(hh);
    writeVideo(vid,frame);
end
%End video write
close(vid)

%%
videofilename = 'small_ice_flextest.avi';
vid = VideoWriter(videofilename,'Uncompressed AVI');
vid.FrameRate = 1; % frames per second
open(vid)

polyx  = [west_x fliplr(west_x)];         % x locations of my polygon
for ab = 1:TKM  
fsPolyys(1:1450,ab) = [Delta_hs3(1:TKM,ab).' fliplr(F_beds(ab,1:TKM))]; % y locations of ice polygon
fsPolyyb(1:1450,ab) = [F_beds(ab,1:TKM) -1200*ones(size(F_beds(1,1:TKM)))];        % y locations of my bed polygon (rock)
end
Small_names = {'LGM ','12 Ka  ','11.6 Ka  ','10.4 Ka ','9.1 Ka  ','8.1 Ka  ','7.3 ka  ','4 Ka ','2 Ka','Current '};
Small_steps = [152 268 300 328 375 402 427 488 662 TKM];
Small_stepsx = [230 268 300 328 375 402 427 488 662 TKM];
fTime_steps = [149 264 280 326 376 414 445 572 649 725];
fspolywater=zeros(1,1450);
fsPolyyw =zeros(1,1450);
for bb = 1:10
% for bb = 1:3
fspolywater(bb,1:2*Small_stepsx(bb)) = [1:Small_stepsx(bb) fliplr(1:Small_stepsx(bb))];
fsPolyyw(bb,1:2*Small_stepsx(bb)) = [sea(1:Small_stepsx(bb),Small_stepsx(bb)).' min(fliplr(F_beds(Small_stepsx(bb),1:Small_stepsx(bb))),sea(1:Small_stepsx(bb),Small_stepsx(bb)).')];
end
hh=figure; clf; hold on
for iii=1:10
   clf; hold on; 
        p=fill(fspolywater((iii),:)*1000,fsPolyyw((iii),:),b);
        hhh(iii) =fill(polyx,fsPolyys(:,Small_steps(iii)),icecolor); 
        hold on; fill(polyx,fsPolyyb(:,Small_stepsx(iii)),rockcolor);hold on;
        p.LineStyle='none';
        title([sprintf('%d',(round(TIME(fTime_steps(iii)),-2))),' years before present for Minimum LGM Extent'])
        text(5.5e5,-500,'Bedrock','fontsize',12);
        text(5.5e5,1000,'Ice Sheet','fontsize',12);
        text(5e4,-500,'Ocean','fontsize',8);
        xlabel('Distance (km)');ylabel('Elevation (m)');
        axis([1 (7.25*10^5) -1000 4000]);
        xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/1000) 
          frame = getframe(hh);
    writeVideo(vid,frame);
end
%End video write
close(vid)

%%
% videofilename = 'small_ice_flextest.avi';
% vid = VideoWriter(videofilename,'Uncompressed AVI');
% vid.FrameRate = 2; % frames per second
% open(vid)
% polyx  = [west_x fliplr(west_x)];         % x locations of my polygon
% for ab = 1:TKM  
% fsPolyy(1:1450,ab) = [Delta_hs3(1:TKM,ab).' fliplr(F_beds(ab,1:TKM))]; % y locations of ice polygon
% fsPolyyb(1:1450,ab) = [F_beds(ab,1:TKM) -1200*ones(size(F_beds(1,1:TKM)))];        % y locations of my bed polygon (rock)
% end
% %fBig_names = {'LGM ','12 Ka  ','11.6 Ka  ','10.4 Ka ','9.1 Ka  ','8.1 Ka  ','7.3 ka  ','4-2 Ka ','Current '};
% fSmall_steps = [26    51    76   101   126   151   176   201   226   251   276   301   326   351   376   401   426   451   476 501   526   551   576   601   626   651   676   701 725 725 725 725 725 725 725];
% fSMALL_steps = [151    151    151   151   151   151   176   201   226   251   276   301   326   351   376   401   426   451   476 501   526   551   576   601   626   651   676   701 725 725 725 725 725 725 725];
% %fSMALL_stepsy = [171    171    171   171   171   171  235 245  255   265   276   301   326   351   376   401   426   451   476 501   526   551   576   601   626   651   676   701 725 725 725 725 725 725 725];
% % Small_steps = [152 268 300 328 375 402 427 488 435];
% % Small_stepsx = [171 268 300 328 375 402 427 488 435];
% fspolywater=zeros(1,1450);
% fsPolyyw =zeros(1,1450);
% for bb = 1:35
% fspolywater(bb,1:2*fSMALL_steps(bb)) = [1:fSMALL_steps(bb) fliplr(1:fSMALL_steps(bb))];
% fsPolyyw(bb,1:2*fSMALL_steps(bb)) = [sea(1:fSMALL_steps(bb),fSMALL_steps(bb)).' min(fliplr(F_beds(fSMALL_steps(bb),1:fSMALL_steps(bb))),sea(1:fSMALL_steps(bb),fSMALL_steps(bb)).')];
% end
% hh=figure; clf; hold on
% for iii=1:35
%    clf; hold on; 
%         hhh(iii) =fill(polyx,fsPolyy(:,fSmall_steps(iii)),icecolor); 
%         hold on; fill(polyx,fsPolyyb(:,fSmall_steps(iii)),rockcolor);hold on;
%         p=fill(fspolywater((iii),:)*1000,fsPolyyw((iii),:),b);
%         p.LineStyle='none';
%         title([sprintf('%d',(round(TIME(fSmall_steps(iii)),-2))),' years before present for Minimum LGM Extent'])
%         text(5.5e5,-500,'Bedrock','fontsize',12);
%         text(5.5e5,1000,'Ice Sheet','fontsize',12);
%         text(5e4,-500,'Ocean','fontsize',8);
%         xlabel('Distance (km)');ylabel('Elevation (m)');
%         axis([1 (7.25*10^5) -1000 4000]);
%         xt = get(gca, 'XTick');                                 % 'XTick' Values
% set(gca, 'XTick', xt, 'XTickLabel', xt/1000) 
%           frame = getframe(hh);
%     writeVideo(vid,frame);
% end
% %End video write
% close(vid)
%%
% close all