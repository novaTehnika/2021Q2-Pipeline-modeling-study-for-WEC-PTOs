function [q_in, q_out, dydt_pLine, pLsoln] = ...
         pipelineNPi_uFric(t,y,p_in,p_out,par,lineID,postProcess)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pipelineNPi.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 5/28/2021
%
% PURPOSE/DESCRIPTION:
% 
%
% FILE DEPENDENCY: NA
%
% UPDATES:
% 5/28/2021 - created based on pipelineNPi.m and the
% usteady friction model implimented in pipelinePi_uFric.m.
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

    % calculate state indices
    ny = (2 + par.nyf)*par.n_seg(lineID)-1; % number of states
    
    iq = 1:2:2*par.n_seg(lineID)-1;
    iyf = ones(par.n_seg(lineID),1)*(0:par.nyf-1) ...
          + (2*par.n_seg(lineID):par.nyf:ny)'*ones(1,par.nyf);
    ip = 2:2:2*par.n_seg(lineID)-1;
    
    % load state arrays for easier indexing
    q = y(iq);
    yf = y(iyf);
    p = [y(ip); p_out];
    
    % output boundary conditions to system
    q_in = q(1);
    q_out = q(end);

    % intialize
    dydt_pLine = zeros(ny,1); 
    Js = zeros(par.n_seg(lineID),1);
    Ju = zeros(par.n_seg(lineID),1);
    
    % T = par.rho*par.d_line^2/4/par.mu; % unsteady friction time scale

    % dp/dt and dq/dt at remaining pi-lumps
    for j = 1:par.n_seg(lineID)
        % dq/dt
         % Steady friction term
        Js(j) = flowR(q(j))*q(j);
         % Unsteady friction term
        Ju(j) = 4*par.rho*(par.L_line/par.n_seg(lineID))...
                /par.T*sum(par.m'.*yf(j,:));
        
        if j ~= 1
            dydt_pLine(iq(j)) = (p(j-1) - p(j) - Js(j) - Ju(j))...
                                /par.I(lineID); 
        else
            dydt_pLine(iq(j)) = (p_in - p(1) - Js(j) - Ju(j))...
                                /par.I(lineID);            
        end
        
        for k = 1:par.nyf
            dydt_pLine(iyf(j,k)) = dydt_pLine(iq(j))/par.A_line ...
                            - par.n(k)/par.T*yf(j,k);
        end
            
        % dp/dt
        if j ~= par.n_seg(lineID)
            dydt_pLine(ip(j)) = (q(j) - q(j+1)) / lineCap(p(j)); 
        end
    end

    % Calculate friction loss
    if postProcess
        % calc. total friction loss
        pLsoln.PPs = sum(abs(Js.*q));
        pLsoln.PPu = sum(abs(Ju.*q));
        pLsoln.PPfric = pLsoln.PPs + pLsoln.PPu; 
        
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