function pLsoln = mocSimple(t,pLsoln_old,p_in,p_out,ID,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mocSimple.m function m-file
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

    % calculate the number of nodes based on the number of segments
    nx = 2*par.n_seg(ID)+1; 
    % calculate the distance between nodes based on the number of nodes
    dx = par.L_line/(nx-1);
    
    I = par.a(ID)*par.rho/par.A_line; % line inertance, 
                                      % commonly used in MOC soln.
    R = par.rho*dx/(2*par.d_line*par.A_line^2); % constant variables in 
                                                %resistance term, commonly 
                                                % used in MOC soln.
    
    % intialize the grid points
    p = zeros(nx,1);
    q = zeros(nx,1);
    
    % load previous state into odd valued grid points (at t = t_last)
    p(1:2:end) = pLsoln_old.p;
    q(1:2:end) = pLsoln_old.q; 
    
    % calculate states for interior nodes at even valued grid points (j==2)
    % (at t = t - dt) 
    % and odd valued nodes(j==2) at time t = t_new  
    for j = [2 3]
        for i = j:2:nx-1
            [p(i), q(i)] = updateInterior(p(i-1),q(i-1),p(i+1),q(i+1));
        end
    end
    
    % calculate states at at t = t + 2dt for boundary nodes
     % inlet boundary (pressure BC)
        Cm_1 = p(2) - I*q(2);
        Bm_1 = I + R*fric(q(2))*abs(q(2));
        p(1) = p_in;
        q(1) = (p(1) - Cm_1)/Bm_1;
     % outlet boundary (pressure BC)
        Cp_nx = p(nx-1) + I*q(nx-1);
        Bp_nx = I + R*fric(q(nx-1))*abs(q(nx-1));
        p(nx) = p_out;
        q(nx) = (Cp_nx - p(nx))/Bp_nx;
        
        % load output structure
        pLsoln.q = q(1:2:nx);
        pLsoln.p = p(1:2:nx);
%         pLstate.PP = sum(R*fric(q(1:nx-1)) ...
%                          .*abs(q(1:nx-1)).*q(1:nx-1).^2);
        PPfric = zeros(nx-1,1);
        for i = 1:nx-1
            PPfric(i) = R*fric(q(i)).*abs(q(i)).*q(i).^2;
        end
        pLsoln.PPfric = sum(PPfric);
        
    function [pP, qP] = updateInterior(pA,qA,pB,qB)
        Cp = pA + I*qA;
        Bp = I + R*fric(qA)*abs(qA);
        Cm = pB - I*qB;
        Bm = I + R*fric(qB)*abs(qB);
        qP = (Cp - Cm)/(Bp + Bm);
        pP = Cp - Bp*qP;
    end

    function f = fric(q)
    % Purpose: calculate the Darcy friction factor for either laminar or 
    % turbulent flow   
    
        Re = par.rho*abs(q)/par.A_line*par.d_line/par.mu;
        
        % Parameters of the pipe friction model 
         % with linear transisiton in transitional region of laminar and
         % turbulent regimes
        Re1 = 2300; % transition form laminar to transitional flow
        f1 = 64/Re1; % friction factor at transisiton
        Re2 = 4500; % transition form transitional to laminar flow
        f2  = 0.316*Re2^-0.25; % friction factor at transisiton

        % Calculate darcy friction factor based on Reynolds number
        if Re < Re1  % Laminar flow regime
            f = 64/max(Re,1e-3);
        elseif Re < Re2 % transitional flow 
            f = f1 + (f2-f1)/(Re2 - Re1)*(Re - Re1);
        else % turbulent flow regime
            f = 0.316*Re^-.25; % Blasius correlation
        end
    end
    
end
