function q = pumpFlow(t,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pumpFlow.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 3/29/2021
%
% PURPOSE/DESCRIPTION:
% Calculate the time dependent flow rate of the pump
%
% FILE DEPENDENCY: NA
%
% UPDATES:
% 3/29/2021 - created.
%
% Copyright (C) 2022  Jeremy W. Simmons II
% 
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program. If not, see <https://www.gnu.org/licenses/>.
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if par.irreg % irregular sea-state
        persistent Shat_sqrt w phi dw q_  %#ok<TLEV>
        if isempty(w)
            nw = floor((par.wend-par.w1)/par.dw)+1;
            dw = (par.wend-par.w1)/(nw-1);
            w = linspace(par.w1,par.wend,nw);
            phi = 2*pi*(rand(1,nw)-0.5);
            Shat_sqrt = sqrt(5*(pi/par.Tp)^4./w.^5.*exp(-20*(pi/par.Tp./w).^4));
            q_ = @(t) abs(par.Xq.*sum(w.*Shat_sqrt.*sqrt(dw).*sin(w*t + phi)));
        end
        q = q_(t);
        
    else % regular sea-state
        q = abs(par.Xq*sin(2*pi/par.Tp*t)); % rectified sinusoidal flow

    end

end
