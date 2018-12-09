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

function K=Kfun(Swf,Ks,m)


K = Ks*Swf.^0.5.*(1-(1-Swf.^(1/m)).^m).^2;

end