%% Constants used in the model
% This script is used to store constants that are used in the model. They are assumed accurate based on Chris Sbarra's research. They are input in this manner to allow editing them to do sensitivity study if desired.

icetau = 110e3;                     % [Pa], yield strength of ice 100 kPa
rho_ice = 917;                      % [kg/m3] ice density
rho_mantle = 3200;                  % [kg/m3] mantle density
gk = 9.81;                           % [m/s2] gravity
hk =40e3;				% [m] Lithosphere elastic thickness: Te ~ h ~ 100 km for cratons
Ek = 70e9;                          % [Pa] Young's Modulus, Turcotte & Schubert p. 152 Problem 3.19
nuk = 0.25;                         % [] Poisson's Ratio, p. 152 Problem 3.19
Dk = Ek*hk^3/12/(1-nuk^2);          % [Pa m^3] Flexural rigidity parameter for lithoplate
alpha = (4*Dk/rho_mantle/gk)^0.25;  % [m] length scale of bulging
