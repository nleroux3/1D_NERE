%=========================================================================%
% 1D Non-equilibrium Richards Equation model                              %
% Copyright (C) 2018 Nicolas R. Leroux
% Author: Nicolas R. Leroux                                               %
% Affiliation: University of Saskatchewan                                 %
% Email: niccolas.r.leroux@gmail.com                                      %
%=========================================================================%
%=========================================================================%
% This file is part of 1D_NERE.                                           %
%                                                                         %
%     1D_NERE is free software: you can redistribute it and/or modify     %
%     it under the terms of the GNU General Public License as published by%
%     the Free Software Foundation, either version 3 of the License, or   %
%     (at your option) any later version.
% 
%     Foobar is distributed in the hope that it will be useful,           %
%     but WITHOUT ANY WARRANTY; without even the implied warranty of      %
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       % 
%     GNU General Public License for more details.                        %
%                                                                         %
%     You should have received a copy of the GNU General Public License   %
%     along with Foobar.  If not, see <https://www.gnu.org/licenses/>.    %
%=========================================================================%

function [psi,thetaR,theta_id,theta_di,theta_s,theta_r,order,wrc,Pdi,...
    Pid, Swf]=hysteresis(theta_1,theta_2,theta_s,theta_r,theta_id,...
    theta_di,order,n,m,alpha,alpha_wetting,thetaS,thetaR,wrc,psi, thetaR_dry,Pdi,Pid)

wrc_ini = wrc;


if wrc_ini == 1  % On main wetting curve
    
    order = 1;
    wrc = 1;
    thetaR = 0;
    theta_s(order) = thetaS;
    theta_r(order) = thetaR;
    theta_id(order) = thetaS;
    theta_di(order) = thetaR;
    Pdi = -inf;
    Pid = 0;
    
    Swf = theta_2/thetaS;
    
    if theta_2 >= theta_1   % Stays on main wetting curve
        
        psi = -((((theta_2 - thetaR)/(thetaS-thetaR))^(-1/m) - 1)^(1/n)) / alpha_wetting ;
        
    else % Moves to scanning drying curve
        
        order = order + 1;
        wrc = 3;
        Pdi = -inf;
        Pid = psi;
        theta_id(order) = theta_1;
        
        thetaR = theta_1 / (1+(thetaS/thetaR_dry-1)*theta_1/thetaS);
        
        theta_di(order) = thetaR;
        
        Swf = 0.5 * ((theta_2/thetaS - thetaR/thetaS) + sqrt((theta_2/thetaS - thetaR/thetaS)^2 + 4/(thetaS/thetaR_dry-1) *(theta_2/thetaS - thetaR/thetaS) ));
        
        Sd_Pdi = 0;
        Sd_Pid = (1+(-alpha * Pid)^n)^(-m);
        
        theta_s(order) = (theta_di(order)*(1 - Sd_Pid) - theta_id(order) *(1 - Sd_Pdi))  /(Sd_Pdi - Sd_Pid);
        theta_r(order) = (Sd_Pdi * theta_s(order) - theta_di(order)) / (Sd_Pdi - 1);
        
        psi = - (((theta_2 - theta_r(order))/(theta_s(order) - theta_r(order)))^(-1/m)-1)^(1/n) /alpha;
        
        if theta_2 > thetaR_dry
            P_w = - (((theta_2 - thetaR_dry)/(thetaS - thetaR_dry))^ (-1/m)-1)^(1/n) /alpha;
            
            if abs(psi) > abs(P_w)  %  Moves to main drying curve
                wrc = 2;
                order = 1;
                thetaR = thetaR_dry;
                psi = P_w;
            end
        end
    end
    return
    
end


if wrc_ini == 3 % On drying scanning curve
    
    if theta_2 <= theta_1 %  stays on drying scanning curve
        
        
        if theta_2 > theta_di(2)
            while theta_2 < theta_di(order)
                order = order - 2;
            end
        else
            theta_2 = theta_di(2) + 1e-5;
        end
        
        psi = - (((theta_2 - theta_r(order))/(theta_s(order) - theta_r(order)))^(-1/m)-1)^(1/n)/alpha;
        
        Swf = 0.5 * ((theta_2/thetaS - thetaR/thetaS) + sqrt((theta_2/thetaS - thetaR/thetaS)^2 + 4/(thetaS/thetaR_dry-1) *(theta_2/thetaS - thetaR/thetaS) ));
        
        if theta_2 > thetaR_dry
            P_w = - (((theta_2 - thetaR_dry)/(thetaS - thetaR_dry))^ (-1/m)-1)^(1/n) /alpha;
            
            if abs(psi) > abs(P_w)  %  Moves to main drying curve
                wrc = 2;
                order = 1;
                thetaR = thetaR_dry;
                psi = P_w;
            end
        end
        
    else  % WRC moves to wetting scanning (di: drainage to imbibition)
        
        wrc = 4;
        
        order = order + 1 ;
        Pdi = psi;
        theta_di(order) = theta_1;
        theta_id(order) = theta_id(order-1);
        
        Si_Pdi = (1 + (-alpha_wetting * Pdi)^n)^(-m);
        Si_Pid = (1 + (-alpha_wetting * Pid)^n)^(-m);
        
        theta_s(order) = (theta_di(order)*(1 - Si_Pid) - theta_id(order) * (1 - Si_Pdi))  / (Si_Pdi - Si_Pid);
        
        theta_r(order) = (Si_Pdi * theta_s(order) - theta_di(order)) / (Si_Pdi - 1);
        
        Swf = 0.5 * ((theta_2/thetaS - thetaR/thetaS) + sqrt((theta_2/thetaS - thetaR/thetaS)^2 + 4/(thetaS/thetaR_dry-1) *(theta_2/thetaS - thetaR/thetaS) ));
        
        while theta_2 > theta_id(order)
            order = order - 2;
        end
        
        if order > 5
            order = order - 2;
        end
        
        psi = - (((theta_2 - theta_r(order))/(theta_s(order) - theta_r(order)))^(-1/m)-1)^(1/n)/alpha_wetting;
        
        P_w = - ((theta_2 /thetaS)^(-1/m)-1)^(1/n)/alpha_wetting;
        
        if abs(psi) < abs(P_w) || order == 1 % Moves to main wetting curve
            wrc = 1;
            order = 1 ;
            thetaR = 0;
            psi = P_w;
        end
    end
    return
end



if wrc_ini == 4 %  On wetting scanning scanning
    if theta_2 >= theta_1 %  stays on wetting scanning curve
        
        while theta_2 > theta_id(order)
            order = order - 2;
        end
        
        psi = - (((theta_2 - theta_r(order))/(theta_s(order) - theta_r(order)))^(-1/m)-1)^(1/n) /alpha_wetting;
        
        Swf = 0.5 * ((theta_2/thetaS - thetaR/thetaS) + sqrt((theta_2/thetaS - thetaR/thetaS)^2 + 4/(thetaS/thetaR_dry-1) *(theta_2/thetaS - thetaR/thetaS) ));
        
        P_w = - ((theta_2 /thetaS)^(-1/m)-1)^(1/n)/alpha_wetting;
        
        if abs(psi) < abs(P_w) || order == 1 %  Moves to main wetting curve
            wrc = 1;
            order = 1 ;
            thetaR = 0;
            psi = P_w;
        end
        
        
        
    else  % WRC moves to drying scanning curve (id : imbibition to dyring)
        
        order  = order + 1;
        wrc = 3;
        
        Pid = psi;
        theta_id(order) = theta_1 ;
        theta_di(order) = theta_di(order-1);
        
        if abs(Pid-Pdi) < 1e-9
           Pid = Pdi + 1e-5;
           
        end
        
        Swf = 0.5 * ((theta_2/thetaS - thetaR/thetaS) + sqrt((theta_2/thetaS - thetaR/thetaS)^2 + 4/(thetaS/thetaR_dry-1) *(theta_2/thetaS - thetaR/thetaS) ));
        
        Sd_Pdi = (1 + (-alpha * Pdi)^n)^(-m);
        Sd_Pid = (1 + (-alpha * Pid)^n)^(-m);
        
        
        
        theta_s(order) = (theta_di(order)*(1 - Sd_Pid) - theta_id(order) * (1 - Sd_Pdi))  / (Sd_Pdi - Sd_Pid);
        
        theta_r(order) = (Sd_Pdi * theta_s(order) - theta_di(order)) / (Sd_Pdi - 1);
        
        if theta_2 > theta_di(2)
            while theta_2 < theta_di(order)
                order = order - 2;
            end
        else
            theta_2 = theta_di(2) + 1e-5;
        end
        
        if order > 5
            order = order - 2;
        end
        
        psi = - (((theta_2 - theta_r(order))/(theta_s(order)-  theta_r(order)))^(-1/m)-1)^(1/n) /alpha;
        
        if theta_2 > thetaR_dry
            P_w = - (((theta_2 - thetaR_dry)/(thetaS - thetaR_dry))^ (-1/m)-1)^(1/n) /alpha;
            
            if abs(psi) > abs(P_w)  %  Moves to main drying curve
                wrc = 2;
                order = 1;
                thetaR = thetaR_dry;
                psi = P_w;
            end
        end
        
    end
    return
end

if wrc_ini == 2 % On main drying curve
    order = 1;
    
    thetaR = thetaR_dry;
    theta_s(order) = thetaS;
    theta_r(order) =  thetaR;
    theta_di(order) = thetaR;
    theta_id(order) = thetaS;
    Pdi = -inf;
    Pid = 0;
    
    if theta_2 > theta_1 % WRC moves to wetting scanning (di: drainage to imbibition)
        
        
        order = order + 1;
        wrc = 4;
        
        Pdi = psi;
        Pid = 0;
        theta_di(order) = theta_1 ;
        theta_id(order) = thetaS;
        
        Swf = 0.5 * ((theta_2/thetaS - thetaR/thetaS) + sqrt((theta_2/thetaS - thetaR/thetaS)^2 + 4./(thetaS/thetaR_dry-1) *(theta_2/thetaS - thetaR/thetaS) ));
        
        Si_Pdi = (1 + (-alpha_wetting * Pdi)^n)^(-m);
        Si_Pid = (1 + (-alpha_wetting * Pid)^n)^(-m);
        
        theta_s(order) = (theta_di(order)*(1 - Si_Pid) - theta_id(order) * (1. - Si_Pdi))  / (Si_Pdi - Si_Pid);
        
        theta_r(order) = (Si_Pdi * theta_s(order) - theta_di(order)) / (Si_Pdi - 1);
        
        psi = - (((theta_2 - theta_r(order))/(theta_s(order) - theta_r(order)))^(-1/m)-1)^(1/n) /alpha_wetting;
        
        P_w = - ((theta_2 /thetaS)^(-1/m)-1)^(1/n)/alpha_wetting;
        
        if abs(psi) < abs(P_w) || order == 1 % Moves to main wetting curve
            wrc = 1;
            order = 1 ;
            psi = - ((theta_2 /thetaS )^(-1/m)-1)^(1/n) /alpha_wetting;
        end
        
        
    else  % Stays on main drying curve
        
        Swf = 0.5 * ((theta_2/thetaS - thetaR/thetaS) + sqrt((theta_2/thetaS - thetaR/thetaS)^2 + 4/(thetaS/thetaR_dry-1) *(theta_2/thetaS - thetaR/thetaS) ));
        psi = -((((theta_2 - thetaR)/(thetaS-thetaR))^(-1/m) - 1) ^(1/n)) / alpha  ;
        
    end
    return
end


end