function pLsoln = mocSimple_uFric(t,pLsoln_old,p_in,p_out,ID,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mocSimple_uFric.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 5/17/2021
%
% PURPOSE/DESCRIPTION:
% Adapted from mocSimple.m. Differences are: Both updated
% timesteps are returned in pLsoln.p and pLsoln.q. Therefore pLsoln_old
% contains two previous time step solution so that the t-2dt flow rate and
% y array can be used in calculating the new y array.
%
% FILE DEPENDENCY: NA
%
% UPDATES:
% 5/17/2021 - created based on mocSimple.m
% 5/25/2021 - Added update to accumulated sums at boundary.
% 6/15/2021 - Corrected unsteady friction model for integration along
%           characteristics per derivation given in 
%           "Transient friction equation sheet.docx" 
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
% if t > 100 && ID == 1
%     aaa = 99; 
% end
    % calculate the number of nodes based on the number of segments
    nx = 2*par.n_seg(ID)+1; 
    % calculate the distance between nodes based on the number of nodes
    dx = par.L_line/(nx-1);
    
    I = par.a(ID)*par.rho/par.A_line; % line inertance, 
                                      % commonly used in MOC soln.
    R = par.rho*dx/(2*par.d_line*par.A_line^2); % constant variables in 
                                                %resistance term, commonly 
                                                % used in MOC soln.
    tau = par.dt_MOC(ID)/par.T;
    Bf = 2*par.a(ID)*par.rho*tau ...
            *sum(par.m.*((1 - exp(-par.n*tau)) ...
            ./(par.A_line*par.n*tau))); 
        
    % intialize the grid points
    
    PPs = zeros(nx,1); % ********************************
    PPu = zeros(nx,1);
    
    % load previous states into odd valued grid points (at t = t_last)
    p = pLsoln_old.p; 
    q = pLsoln_old.q;
    y = pLsoln_old.y;  % ********************************
    
    
    % calculate states for interior nodes at even valued grid points (j==2)
    % (at t = t - dt) 
    % and odd valued nodes(j==2) at time t = t_new  
    for j = [2 3]
        for i = j:2:nx-1
            [p(i), q(i), y(i,:), PPs(i), PPu(i)] = ...
                updateInterior(p(i-1),q(i-1),p(i+1),q(i+1),q(i),...
                y(i-1,:),y(i,:),y(i+1,:)); % ********************************
        end
    end
    
    % calculate states at at t = t + 2dt for boundary nodes
     % inlet boundary (pressure BC)
        Cm_1 = p(2) - I*q(2) + Cf_(q(1), y(2,:), y(1,:));
        Bm_1 = I + R*fric(q(2))*abs(q(2)) + Bf;
        p(1) = p_in;
        q_old = q(1);
        q(1) = (p(1) - Cm_1)/Bm_1;
        for k = 1:5
            dy = (1 - exp(-par.n(k)*tau))/(par.A_line*par.n(k)*tau)...
                 *(q(1) - q_old);
            y(1,k) = exp(-par.n(k)*tau)*y(1,k) + dy;
        end
     % outlet boundary (pressure BC)
        Cp_nx = p(nx-1) + I*q(nx-1) - Cf_(q(nx), y(nx-1,:), y(nx,:));
        Bp_nx = I + R*fric(q(nx-1))*abs(q(nx-1)) + Bf;
        p(nx) = p_out;
        q_old = q(nx);
        q(nx) = (Cp_nx - p(nx))/Bp_nx;
        for k = 1:5
            dy = (1 - exp(-par.n(k)*tau))/(par.A_line*par.n(k)*tau)...
                 *(q(nx) - q_old);
            y(nx,k) = exp(-par.n(k)*tau)*y(nx,k) + dy;
        end
        
     % load output structure
        pLsoln.q = q;
        pLsoln.p = p;
        pLsoln.y = y;
        
        pLsoln.PPs = sum(PPs(1:nx-1));
        pLsoln.PPu = sum(PPu(1:nx-1));
        pLsoln.PPfric = pLsoln.PPs + pLsoln.PPu;
        
    function [pP, qP, yP, PPs, PPu] = updateInterior(pA,qA,pB,qB,qPminus,...
                                    yA,yPminus,yB)
        Cpf = Cf_(qPminus, yA, yPminus);
        Cmf = Cf_(qPminus, yB, yPminus);
        Cp = pA + I*qA - Cpf;
        Bp = I + R*fric(qA)*abs(qA) + Bf;
        Cm = pB - I*qB + Cmf;
        Bm = I + R*fric(qB)*abs(qB) + Bf;
        qP = (Cp - Cm)/(Bp + Bm);
        pP = Cp - Bp*qP;
        
        yP = zeros(size(yPminus));
        for k = 1:5
            dyP = (1 - exp(-par.n(k)*tau))/(par.A_line*par.n(k)*tau)...
                 *(qP - qPminus);
            yP(k) = exp(-par.n(k)*tau)*yPminus(k) + dyP;
        end
        PPs = R*fric(qP).*abs(qP^3);
        PPu = (Cpf + Bf*qP)*qP;
    end

    function Cf = Cf_(qPminus, yI, yPminus)
        Cf = 2*par.a(ID)*par.rho*tau ...
            *sum(par.m.*(yI' + exp(-par.n*tau).*yPminus' ...
            - (1 - exp(-par.n*tau)) ...
            ./(par.A_line*par.n*tau)*qPminus));
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
