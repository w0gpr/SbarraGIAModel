% Bending of the Lithosphere under a Triangular Load
% Turcotte & Schubert Section 12.1
% example code plate1bvp.m

function sol = plate1bvp_v2(makeplot)


infinity = 20;  % This means that the half-width of the plate, L, is 20x the flexural parameter, alpha
% For a line load like the Greenland Ice Sheet, L ?~ 200 km and alpha is
% given by Turcotte & Schubert Equation 3.127,
% alpha = (4D/d_rho/g)^0.25
% for density difference d_rho = rho_mantle - rho_ice (ice is what the
% plate is being loaded by and the mantle is what it is displacing)
% for flexural rigitidy D = (E*h^3 / 12 / (1-nu^2)) Equation 3.72
% for plate thickness h (the thickness of lithosphere or crust)
% and Young's Modulus E and Poisson's ratio nu
% Look these up for the lithosphere...
solinit = bvpinit(linspace(0,infinity,1001),[1 0 0 1]);
options = bvpset('stats','off');

sol = bvp4c(@plateode,@platebc,solinit,options);


if makeplot,
    eta = sol.xk;
    fk = sol.yk;
    clf reset
    figure(1); clf;
    subplot(2,1,1)
    plot(eta,-fk(1,:))
    axis([0 infinity -1.1 0.1]);
    title('Plate bending due to a line load')
    xlabel('x - units: $\alpha$')
    ylabel('y - units: $V_0 \alpha^3 / 8 / D$')
    % the line load has magnitude V_0, which is in Pa*m = kg/s2 = N/m
    % V_0 = rho_ice * g * length^2, is it ice thickness^2?
    shg % show graph window
    
    hk = 40e3;  % m; Lithosphere elastic thickness: Te ~ h ~ 100 km for cratons
    Ek = 70e9;   % Pa; Young's Modulus, p. 152 Problem 3.19
    nu = 0.25;  % []; Poisson's Ratio, p. 152 Problem 3.19
    Dk = Ek*hk^3/12/(1-nu^2); % Pa m^3; Flexural rigidity parameter
    rho_mantle = 3200; % kg/m3 rho is density
    rho_ice = 917; % kg/m3
    d_rho = rho_mantle - rho_ice;
    gk = 9.8;    % m/s2; gravity acceleration
    alpha = (4*Dk/d_rho/gk)^0.25;
%     alpha = (4*Dk/rho_mantle/gk)^0.25;
    
     H_k = 1574;  % m; average ice thickness guess

%     L = 1500e3; % m; ice sheet width guess; this is from Chris's slides
    V0k = rho_ice*gk*H_k*Lk(1,:);  % kg/s2; line load
    w0k = V0k*alpha^3/8/Dk;
    fprintf('ice causes elastic depression w0 %1.1f km and forebulge %1.1f m at %1.1f km from center of ice sheet\n',w0k/1e3,w0k*exp(-pi),pi*alpha/1e3)
    
    subplot(2,1,2)
    plot(eta*alpha/1e3,-fk(1,:)*w0k)
    %axis([0 5*alpha/1e3 -0.1*w0 0.1*w0]);
    title('Plate bending due to a line load')
    xlabel('x (km)')
    ylabel('y (meters)')
    
end


function dfdeta = plateode(eta,fk)
    dfdeta = [ fk(2)
               fk(3)
               fk(4)
                -4*fk(1)];

function res = platebc(f0,finf)
    res =  [f0(2)
            f0(4)-4
            finf(2)
            finf(3)];