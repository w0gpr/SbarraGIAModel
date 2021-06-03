% approximate GrIS as a series of line loads
% run plate1bvp.m (adapted from Turcotte & Schubert) 
% and scale its output (vertical displacement) to each line load, 
% add together the to get elastic forebulge from whole ice sheet
%
% kristin poinar april 16 2020 for chris sbarra MS work
%
% is this going to be meaningfully better than treating the ice sheet as a
% rectangular block?  we'll see.
%
% clear variables
clear
% constants
dx = 1000;    % [m] choose some suitably small delta x

Lmax  = (725-152)*1000*2; % [m] full width of the ice sheet at maximum (edge to edge)
Lmin = (725-455)*1000*2; % ]m] full width of the ice sheet at minimum (min margin)
Lmarg = Lmax:-dx:Lmin;
% up = zeros(607,11461);
% surf = zeros(1,11461);
% Sis_surf_flip= importdata('Sis_surf_flip.xlsx')';
% surf(5441:6021) =Sis_surf_flip; 
%%
for i = 1:607
Lk = Lmax:-dx:Lmax-((i-1)*1000);
Lk = Lk';

icetau = 90e3;  % [Pa], yield strength of ice 100 kPa
rho_ice = 917; % [kg/m3] ice density
rho_mantle = 3200; % [kg/m3] mantle density
gk = 9.8;    % [m/s2] gravity
hk = 80e3;  % [m] Lithosphere elastic thickness: Te ~ h ~ 100 km for cratons
Ek = 70e9;   % [Pa] Young's Modulus, Turcotte & Schubert p. 152 Problem 3.19
nuk = 0.25;  % [] Poisson's Ratio, p. 152 Problem 3.19
Dk = Ek*hk^3/12/(1-nuk^2); % [Pa m^3] Flexural rigidity parameter for lithoplate
alpha = (4*Dk/rho_mantle/gk)^0.25; % [m] length scale of bulging
%
% horizontal space
xk = -5*Lmax : dx : 5*Lmax;

% parabolic ice sheet profile
Hk = real(sqrt(2*icetau / rho_ice / gk * (Lk(i)/2-xk)));
Hk(xk<0) = fliplr(Hk(xk>0));
up(i,:) = Hk(:);
end
Hk = up(1:2:end,:);
%%
% here you put in you for loop for creating the animation
% first initialize your output matrices to save speed
% then run the for loop

% load V(x)
out = zeros(304, 11461);
for i = 1:304
Vk = Hk(i,:)*dx;  % [Pa/m] spatially varying line load V(x)

% Flexure from a generic load
gensol = plate1bvp_v2(false);  % false will not make any plots
% gensol.x is x in units of alpha
% gensol.y is w in units of (V_0 alpha^3) / (8 D)
xscaled = gensol.x * alpha;  % now this is in meters

% Step through all the line loads V and scale the flexure solution to them.
% They can then be superimposed to get the full flexure from the ice sheet.
nx = length(xk);
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

%%
%out = out-(out(304,:));
valk = 1:304;
figure; 
plot(xk/1e3,[Hk(valk,:); out(valk,:)])
xlabel('kilometers')
ylabel('meters')
legend('ice sheet','lithosphere plate')

