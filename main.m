%=========================================================================%
% 1D Non-equilibrium Richards Equation model                              %
% Copyright (C) 2018 Nicolas R. Leroux
% Author: Nicolas R. Leroux                                               %
% Affiliation: University of Saskatchewan                                 %
% Email: niccolas.r.leroux@gmail.com                                      %
%                                                                         %
% DESCRIPTION: This is a 1-D non-Equilibriums Richards equation solver    %
%                                                                         %
% The formulation is based on a forward Euler implementation of the       %
% mixed theta-head formulation of the 1-D non equilibrium Richards        %
% Equation. A finite difference scheme is used to discretize the PDE.     %
%                                                                         %
% This code requires the functions 'hysteresis.m' and 'Kfun.m'.           %                                                     %
%                                                                         %
% Please cite the following publication if you are using this code:       %
%    Leroux, N.R. and J.W. Pomeroy (2018), Simulation of Capillary        %
%    Overshoot in Snow Combining Trapping of the Wetting Phase with a     %
%    Non-Equilibrium Richards Equation Model, Water Resources Research    %               
%                                                                         %                 
%                                                                         %
%=========================================================================%
%=========================================================================%
% This program is free software: you can redistribute it and/or modify    %
%  it under the terms of the GNU General Public License as published      %
%  by the Free Software Foundation, either version 3 of the License,      %
%  or (at your option) any later version.                                 %
%                                                                         %
%  This program is distributed in the hope that it will be useful,        %
%  but WITHOUT ANY WARRANTY; without even the implied warranty of         %
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          %
%  GNU General Public License for more details.                           % 
%                                                                         %
%  You should have received a copy of the GNU General Public License      %
%  along with this program.  If not, see <https://www.gnu.org/licenses/>. %
%=========================================================================%


clc
clear
close all

% Model parameters:
gamma = 2;  % Scaling factor for boundary wetting curve [-]
tau_0 = 0.1; % Parameter for relaxation coefficient tau [m s]
lambda = 1; % Parameter for relaxation coefficient tau [-]
epsilon = 1e-6; % Initial value of water content if initially dry [m3 m-3]

% Space grid:
z0=0.0;  
zN=0.2;  % Snow depth [m]
dz=0.001;  % Space resolution [m]
z=z0+dz/2:dz:zN-dz/2;  % Spatial grid [m]
nz=numel(z);  % Number of cells

% Snow properties
density = 300;   % initial dry snow density [kg/m3]
Dgrain = 1e-3;   % initial grain size [m]


% Hydraulic parameters:
thetaR_dry = 0.039; % Maximum residual water content (for boundary drying curve) [m3 m-3]
thetaS =  1 - density/917;   % Water content at full saturation/pore scale [m3 m-3]

alpha = 4.4e6*(density/Dgrain)^(-0.98);  % van Genuchten parameter for boundary drying curve [m-1]
alpha_wetting = gamma * alpha;   % van Genuchten parameter for boundary wetting curve [m-1]
n = 1+2.7e-3*(density/Dgrain)^(0.61); % van Genuchten parameter
m = 1-1./n;  % van Genuchten parameter
Ks = 9.81 * 1e6 * 3 * (Dgrain*0.5)^2 * exp(-0.013 * density);



% Inital conditions:
thetaR = zeros(1,nz);   % Residual water content [m3 m-3]
theta = zeros(1,nz)+epsilon; % Initial water content [m3 m-3]
order = ones(1,nz);  % Order of the water retention curve (max = 6)
wrc = ones(1,nz);   % Index for water retention curve (1=boundary wetting curve, 2 = boundary drying curve, 3 = scanning drying curve, 4 = scanning wetting curve)
theta_s(1:nz,1:6) = thetaS; % Parameter to compute scanning curve [m3 m-3]
theta_r(1:nz,1:6) = repmat(thetaR',1,6); % Parameter to compute scanning curve [m3 m-3]
theta_id(1:nz,1:6) = thetaS;  % Water content at reversal point from wetting to drainage processes [m3 m-3]
theta_di(1:nz,1:6) = repmat(thetaR',1,6); % Water content at reversal point from drainage to wetting processes [m3 m-3]
Pdi = zeros(1,nz)-inf;  % Capillary pressure wt reversal point from drainage to wetting processes [m]
Pid = zeros(1,nz);    % Capillary pressure wt reversal point from wetting to drainage processes [m]
F = zeros(nz,1);   % Matrix containing gravitational and capillary forces 
psi = zeros(1,nz); % Capillary pressure [m]
Swf = zeros(1,nz); % Flowing saturation [-]



% Boundary conditions
qI = 10e-3/3600; % Input flux [m/s]


% Time grid:
tf = 2500; % [s]
Dt = 0.005;
nt = round(tf/Dt); % Number of timesteps
tPlot = 60; % Time step to save and plot outputs [s]


h = 0;

for i = 1:nz
    [psi(i),thetaR(i),theta_id(i,:),theta_di(i,:),theta_s(i,:),theta_r(i,:),order(i),wrc(i),Pdi(i),Pid(i), Swf(i)] = ...
        hysteresis(theta(i),theta(i),theta_s(i,:),theta_r(i,:),theta_id(i,:),theta_di(i,:),order(i),...
        n,m,alpha,alpha_wetting,thetaS,thetaR(i),wrc(i),psi(i), thetaR_dry,Pdi(i),Pid(i));
end



for t = 1:nt
    
    
    tau = tau_0 * (Swf).^(-lambda);
    
    K = Kfun(Swf,Ks,m);
    
    K_south = [K(1) (K(1:nz-1)+K(2:nz))*0.5] ;
    K_north = [(K(1:nz-1)+K(2:nz))*0.5 K(nz)];
    
    
    M = sparse([1,2:nz-1,nz,2:nz,1:nz-1],[1,2:nz-1,nz,1:nz-1,2:nz],...
        [1 + tau(1)/dz^2 .* K_north(1),...
        1 + tau(2:nz-1)/dz^2 .* (K_south(2:nz-1) + K_north(2:nz-1)),...
        1 + tau(nz)/dz^2 .* K_south(nz),...
        - K_south(2:nz) .* tau(1:nz-1) / dz^2,...
        -K_north(1:nz-1) .* tau(2:nz) / dz^2]);
    
    
    F(nz) = - Dt / dz^2 * K_south(nz) .* (psi(nz) - psi(nz-1)) + Dt / dz * (qI-K_south(nz));
    F(2:nz-1) = Dt / dz^2 * (K_north(2:nz-1).*(psi(3:nz) - psi(2:nz-1)) - K_south(2:nz-1) .* (psi(2:nz-1) - psi(1:nz-2))) + ...
        Dt / dz * (K_north(2:nz-1)-K_south(2:nz-1));
    F(1) = Dt / dz^2 *  K_north(1) * (psi(2) - psi(1)) + Dt / dz * (K_north(1)-K_south(1));
    
    theta_new = (M\F)' + theta;
    
    for i = 1:nz
        [psi(i),thetaR(i),theta_id(i,:),theta_di(i,:),theta_s(i,:),theta_r(i,:),order(i),wrc(i),Pdi(i),Pid(i), Swf(i)] = ...
            hysteresis(theta(i),theta_new(i),theta_s(i,:),theta_r(i,:),theta_id(i,:),theta_di(i,:),order(i),...
            n,m,alpha,alpha_wetting,thetaS,thetaR(i),wrc(i),psi(i), thetaR_dry,Pdi(i),Pid(i));
    end
    
    
    if (mod(t,tPlot/Dt)==0)
        disp(['Time = ', num2str(t*Dt),' s'])
        
        h = h+1;
        
        figure(1), clf;
        hold on
        plot(theta,z,'linewidth',3)
        xlabel('\theta')
        ylabel('Height [m]')
        set(gca,'fontsize',25)
        hold off
        drawnow
        
        
    end
    
    theta = theta_new;
    
end

