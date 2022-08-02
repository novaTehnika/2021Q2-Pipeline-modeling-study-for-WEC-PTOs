function pLsoln = mocDGCM(t,pLsoln_old,p_in,p_out,ID,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mocDGCM.m function m-file
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
% note: indexing is not intuitive.
%       -range of p is [1 nx] with indicies identical to node indices
%       -range of qd and qu are [1 nx-1] with first index for qd identical 
%       to first index of mesh nodes and the first index of qu identical to
%       the second index of the mesh nodes
%       -range of V is [1 nx-2] with the first index identical to the 
%       second index of the mesh nodes
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
    dt = dx/par.a(ID);    
    
    I = par.a(ID)*par.rho/par.A_line; % line inertance, 
                                      % commonly used in MOC soln.
    R = par.rho*dx/(2*par.d_line*par.A_line^2); % constant variables in 
                                                %resistance term, commonly 
                                                % used in MOC soln.
    % intialize the grid points
     % times t and t - dt
    p = zeros(nx,1);
    qd = zeros(nx-1,1);
    qu = zeros(nx-1,1);
    V = zeros(nx-2,1);
    
    % load previous states 
        % temporarily store pressure and flow rates at time t - 2dt into 
        % appropriate indices 
    p(1:2:nx) = pLsoln_old.p;
    qd(1:2:nx-1) = pLsoln_old.qd(1:2:nx-1);
    qu(2:2:nx-1) = pLsoln_old.qu(2:2:nx-1);
    
        % flow rates and gas volume at t - 2dt (odd valued indices) 
        % and t - 3dt  t - 2dt (even valued indices)
    qd_old = pLsoln_old.qd; % downstream flow rate
    qu_old = pLsoln_old.qu; % upstream flow rate
    V_old = pLsoln_old.V;
    
    % calculate states for interior nodes at even valued grid points (j==2)
    % (at t = t - dt) 
    % and odd valued nodes(j==2) at time t = t_new  
    for j = [2 3]
        for i = j:2:nx-1
            [p(i), qd(i), qu(i-1), V(i-1)] = updateInterior(p(i-1),qd(i-1),p(i+1),qu(i),...
                                          qd_old(i),qu_old(i-1),V_old(i-1));
        end
    end
    
    % calculate states at at t = t + 2dt for boundary nodes
     % inlet boundary (pressure BC)
        Cm_1 = p(2) - I*qu(1);
        Bm_1 = I + R*fric(qu(1))*abs(qu(1)); 
        p(1) = p_in;
        qd(1) = (p(1) - Cm_1)/Bm_1;
     % outlet boundary (pressure BC)
        Cp_nx = p(nx-1) + I*qd(nx-1);
        Bp_nx = I + R*fric(qd(nx-1))*abs(qd(nx-1));
        p(nx) = p_out;
        qu(nx-1) = (Cp_nx - p(nx))/Bp_nx;
        
        % load output structure
        pLsoln.qu = qu; % upstream flow rate
        pLsoln.qd = qd; % downstream flow rate
        pLsoln.V = V;
        %q = [qd(1);  qu(2:2:end)];
        pLsoln.q = [qd(1:2:nx-1);  qu(nx-1)];
        pLsoln.p = p(1:2:nx);
        
        % Calculate power loss due to steady friction
%         pLstate.PP = sum(R*fric(qd).*abs(qd).*qd.^2);
        PPfric = zeros(nx-1,1);
        for i = 1:nx-1
            PPfric(i) = R*fric(qd(i)).*abs(qd(i)).*qd(i).^2;
        end
        pLsoln.PPfric = sum(PPfric);

    function [pP, qdP, quP, V] = updateInterior(pA,qdA,pB,quB,q,qu,V_old)
        
        Cp = pA + I*qdA;
        Bp = I + R*fric(qdA)*abs(qdA);
        Cm = pB - I*quB;
        Bm = I + R*fric(quB)*abs(quB);
        Bv = (V_old/(2*dt) + (1-par.phi)*(q-qu))/par.phi;
        B2 = 0.5/(Bp+Bm);
        B1 = -B2*(Bp*Cm + Bm*Cp) + B2*Bm*Bp*Bv + par.p_vap/2;
        C4 = par.C1(ID)*Bp*Bm*B2/(par.phi*dt);
        Bb = C4/B1^2;
        
        if abs(Bb) < 0.001 % Bb is small, use linearize formulation for pP
            if B1 < 0
                pP = -2*B1 - C4/(2*B1) + par.p_vap;
            else % B1 > 0
                pP = C4/(2*B1) + par.p_vap;
            end
        else
            if B1 < 0
                pP = -B1*(1 + sqrt(1+Bb)) + par.p_vap;
            else % B1 > 0
                pP = -B1*(1 - sqrt(1+Bb)) + par.p_vap;
            end
        end
        
        qdP = (pP - Cm)/Bm;
        quP = (Cp - pP)/Bp;
        
%         V = par.C1(ID)/(pP-par.p_vap);
        V = V_old + (par.phi*(qdP - quP) + (1-par.phi)*(q - qu))*2*dt;

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
%         iLam = find(Re <= Re1);
%         iTrans = find(Re > Re1 & Re <= Re2);
%         iTurb = find(Re > Re2);

%         f(Re <= Re1) = 64./max(Re(Re <= Re1),1);
%         f(Re > Re1 & Re <= Re2) = f1 + (f2-f1)/(Re2 - Re1)*(Re(Re > Re1 & Re <= Re2) - Re1);
%         f(Re > Re2) = 0.316*Re(Re > Re2).^-.25; % Blasius correlation
        
        if Re < Re1  % Laminar flow regime
            f = 64./max(Re,1e-3);
        elseif Re < Re2 % transitional flow 
            f = f1 + (f2-f1)/(Re2 - Re1)*(Re - Re1);
        else % turbulent flow regime
            f = 0.316*Re.^-.25; % Blasius correlation
        end

    end
    
end