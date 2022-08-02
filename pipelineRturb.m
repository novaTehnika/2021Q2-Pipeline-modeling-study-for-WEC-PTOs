function [q_in, q_out, dydt_pLine, pLsoln] = ...
         pipelineRturb(t,y,p_in,p_out,par,lineID,postProcess)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pipelineRturb.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 3/25/2021
%
% PURPOSE/DESCRIPTION:
% 
%
% FILE DEPENDENCY: NA
%
% UPDATES:
% 3/25/2021 - created.
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
    
    % solve implicit function for flow (i.e. f(q,Re) = 0 ) 
     %set solver options
    options = optimoptions('fsolve','Display','off'); 

     % inital guess assuming laminar flow for inital guess
    q_o = (p_in-p_out)/(8*pi*par.L_line/par.A_line^2);
    
     % call solver and obtain result
    q = fsolve(@(x) x-(p_in-p_out)/flowR(x),q_o,options); 

    % model outputs
    q_in = q;
    q_out = q;
    dydt_pLine = [];

    % Calculate friction loss
    pLsoln.PPfric = flowR(q)*q^2;

    % Nested functions of pipeline(*)
    
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