function [q_in, q_out, dydt_pLine, pLsoln] = ...
         pipelinePi_uFric(t,y,p_in,p_out,par,lineID,postProcess)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pipelinePi_V01x00.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 5/20/2021
%
% PURPOSE/DESCRIPTION:
% 
%
% FILE DEPENDENCY: NA
%
% UPDATES:
% 5/20/2021 - created based on pipelinePi.m
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

    % model outputs
     % Echo state
    q = y(1);
    q_in = q;
    q_out = q;
    
     % load accumulated sums
    nyf = par.nyf;
    yf = y(2:nyf+1);
     
     % Steady friction term
    Js = flowR(q)*q;
    
     % Unsteady friction term
%     T = par.rho*par.d_line^2/4/par.mu; % unsteady friction time scale 
    
    Ju = 4*par.rho*(par.L_line/par.n_seg(lineID))/par.T*sum(par.m.*yf);
     
     % ODEs
      % accel. of flow
    dqdt = (p_in - p_out - Js - Ju)/par.I(lineID);

      % accumulated values unsteady friction
    dyfdt = zeros(1,nyf);
    for i = 1:nyf
        dyfdt(i) = dqdt/par.A_line - par.n(i)/par.T*yf(i);
    end
    
    dydt_pLine = [dqdt dyfdt];
    
     % Calculate friction loss
    pLsoln.PPs = abs(Js*q);
    pLsoln.PPu = abs(Ju*q);
    pLsoln.PPfric = pLsoln.PPs + pLsoln.PPu;  
            
    function R = flowR(q) 
    % Purpose: calculate flow resistance for either laminar or 
    % turbulent flow 
           
        % parameters of the pipe friction model 
        % with linear transisiton in transitional region of laminar and
        % turbulent regimes
        Re1 = 2300; % transition form laminar to transitional flow
        f1 = 64/Re1; % friction factor at transisiton
        Re2 = 4500; % transition form transitional to laminar flow
        f2  = 0.316*Re2^-0.25; % friction factor at transisiton

        % calculate the Re of flow (set floor to avoid zero result)
        Re = max(4*par.rho*abs(q)/(par.mu*pi*par.d_line),1e-3);

        % Calculate darcy friction factor based on Reynolds number
        if Re < Re1  % Laminar flow regime
            f = 64/max(Re,0.01);
        elseif Re < Re2 % transitional flow 
            f = f1 + (f2-f1)/(Re2 - Re1)*(Re - Re1);
        else % turbulent flow regime
            f = 0.316*Re^-.25; % Blasius correlation
        end
        % calculate flow resistance coefficient
        R = f*Re*par.mu*(par.L_line/par.n_seg(lineID)) ...
                          /(2*par.d_line^2*par.A_line);
    end
    

end