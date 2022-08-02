function [q_in, q_out, dydt_pLine, pLsoln] = ...
         pipelineNPi(t,y,p_in,p_out,par,lineID,postProcess)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pipelineNPi.m function m-file
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

    % load arrays for easier indexing
    q = y(1:2:2*par.n_seg-1);
    p = [y(2:2:2*par.n_seg-1); p_out];

    % output boundary conditions to system
    q_in = q(1);
    q_out = q(end);

    % intialize dy/dt vector
    dydt_pLine = zeros(2*par.n_seg(lineID)-1,1); 

    % dq/dt at first pi-lump
    dydt_pLine(1) = (p_in - p(1) - flowR(q(1))*q(1))...
                /par.I(lineID);

    % dp/dt and dq/dt at remaining pi-lumps
    for j = 1:par.n_seg(lineID)-1 
        % dp/dt
        dydt_pLine(2*j) = (q(j) - q(j+1)) / lineCap(p(j)); 
        % dq/dt
        dydt_pLine(2*j+1) = (p(j) - p(j+1) - flowR(q(j+1))*q(j+1))...
                       /par.I(lineID); 
    end

    % Calculate friction loss
    if postProcess
        
        % intialize sketch variable
        PPfric = zeros(par.n_seg(lineID),1);
        
        % calc. friction loss in each segment
        for i = 1:par.n_seg(lineID)
            PPfric(i) = flowR(q(i))*q(i)^2;
        end

        % calc. total friction loss
        pLsoln.PPfric = sum(PPfric);
        
    else
        pLsoln = [];
        
    end

    %% Nested functions of pipeline(*)
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
    
    function C = lineCap(p)
        % calculate capacitance in pipeline lump
        
        V_line = par.A_line*par.L_line/par.n_seg(lineID); % volume of each lump
        
        % calculate effective bulk modulus 
        
         % via Cho method (as arranged in Yudell, 2017)
%         beta_eff = par.beta* ...
%         (((p/par.p_o)^(1/par.gamma)*exp((par.p_o-p)/par.beta)+par.R) / ...
%         (par.R/par.gamma*par.beta/p+(p/par.p_o)^(1/par.gamma)*...
%         exp((par.p_o-p)/par.beta))); 
         
         % isothermal bulk modulus
        beta_eff = par.beta/(1 + par.beta*(par.R*par.p_o/p^2));
    
        % calculate capacitance
        C = V_line/beta_eff;
    end

end