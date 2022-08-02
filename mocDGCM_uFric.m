function pLsoln = mocDGCM_uFric(t,pLsoln_old,p_in,p_out,ID,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mocDGCM_uFric.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 5/18/2021
%
% PURPOSE/DESCRIPTION:
% Adapted from mocDGCM.m. Differences are: Both updated
% timesteps are returned in pLsoln.p, pLsoln.qu, and pLsoln.qd. Therefore 
% pLsoln_old contains two previous time step solution so that the t-2dt 
% flow rate and yu and yd array can be used in calculating the new y array.
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
% 5/17/2021 - created based on mocDGCM.m
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
    tau = par.dt_MOC(ID)/par.T;
    Bf = 0*2*par.a(ID)*par.rho*tau ...
            *sum(par.m.*((1 - exp(-par.n*tau)) ...
            ./(par.A_line*par.n*tau))); 
        
    % intialize the grid points
     % times t and t - dt
    PPs = zeros(nx,1); 
    PPu = zeros(nx,1);

    
    % load previous states 
        % temporarily store pressure and flow rates at time t - 2dt into 
        % appropriate indices 
    p = pLsoln_old.p;
    qd = pLsoln_old.qd;
    qu = pLsoln_old.qu;
    yd = pLsoln_old.yd;  % ********************************
    yu = pLsoln_old.yu;  % ********************************
    V = pLsoln_old.V;
    
    % calculate states for interior nodes at even valued grid points (j==2)
    % (at t = t - dt) 
    % and odd valued nodes(j==2) at time t = t_new  
    for j = [2 3]
        for i = j:2:nx-1
            [p(i), qd(i), qu(i-1), yd(i,:), yu(i-1,:), V(i-1), PPs(i), PPu(i)] ...
                = updateInterior(p(i-1),qd(i-1),p(i+1),qu(i),...
                                qd(i),qu(i-1),...
                                yd(i-1,:),yd(i,:),...
                                yu(i-1,:),yu(i,:),...
                                V(i-1));
        end
    end
    
    % calculate states at at t = t + 2dt for boundary nodes
     % inlet boundary (pressure BC)
        Cmfd_1 = Cf_(qd(1), yu(1,:), yd(1,:));
        Cm_1 = p(2) - I*qu(1) + Cmfd_1;
        Bm_1 = I + R*fric(qu(1))*abs(qu(1)); 
        p(1) = p_in;
        qd_old = qd(1);
        qd(1) = (p(1) - Cm_1)/Bm_1;
        for k = 1:5
            dyd = (1 - exp(-par.n(k)*tau))/(par.A_line*par.n(k)*tau)...
                 *(qd(1) - qd_old);
            yd(1,k) = exp(-par.n(k)*tau)*yd(1,k) + dyd;
        end
        PPs(1) = R*fric(qd(1)).*abs(qd(1)^3);
        PPu(1) = (Cmfd_1 + Bf*qd(1))*qd(1);
     % outlet boundary (pressure BC)
        Cpfu_nx = Cf_(qu(nx-1), yd(nx-1,:), yu(nx-1,:)); %************** check index
        Cp_nx = p(nx-1) + I*qd(nx-1) - Cpfu_nx;
        Bp_nx = I + R*fric(qd(nx-1))*abs(qd(nx-1));
        p(nx) = p_out;
        qu_old = qu(nx-1);
        qu(nx-1) = (Cp_nx - p(nx))/Bp_nx;
        for k = 1:5
            dyu = (1 - exp(-par.n(k)*tau))/(par.A_line*par.n(k)*tau)...
                 *(qu(nx-1) - qu_old);
            yu(nx-1,k) = exp(-par.n(k)*tau)*yu(nx-1,k) + dyu;
        end
        
     % load output structure
        pLsoln.qd = qd; % downstream flow rate
        pLsoln.qu = qu; % upstream flow rate
        pLsoln.yd = yd; 
        pLsoln.yu = yu; 
        pLsoln.V = V;
        pLsoln.q = [qd(1:nx-1);  qu(nx-1)];
        pLsoln.p = p(1:nx);
        
        % Calculate power loss due to steady friction
        pLsoln.PPs = sum(PPs(1:nx-1));
        pLsoln.PPu = sum(PPu(1:nx-1));
        pLsoln.PPfric = pLsoln.PPs + pLsoln.PPu;

    function [pP, qdP, quP, ydP, yuP, V, PPs, PPu] ...
                    = updateInterior(pA,qdA,pB,quB,...
                                qdPminus,quPminus,...
                                ydA,ydPminus,...
                                yuPminus,yuB,...
                                V_old)

        Cpfu = Cf_(quPminus, ydA, yuPminus);
        Cmfd = Cf_(qdPminus, yuB, ydPminus);
        Cp = pA + I*qdA - Cpfu;
        Bp = I + R*fric(qdA)*abs(qdA) + Bf;
        Cm = pB - I*quB + Cmfd;
        Bm = I + R*fric(quB)*abs(quB) + Bf;
        Bv = (V_old/(2*dt) + (1-par.phi)*(qdPminus-quPminus))/par.phi;
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
        V = V_old + ...
            (par.phi*(qdP - quP) + (1-par.phi)*(qdPminus - quPminus))*2*dt;

        ydP = zeros(size(ydPminus));
        yuP = zeros(size(yuPminus));
        for k = 1:5
            dydP = (1 - exp(-par.n(k)*tau))/(par.A_line*par.n(k)*tau)...
                 *(qdP - qdPminus);
            ydP(k) = exp(-par.n(k)*tau)*yd(k) + dydP;
            
            dyuP = (1 - exp(-par.n(k)*tau))/(par.A_line*par.n(k)*tau)...
                 *(quP - quPminus);
            yuP(k) = exp(-par.n(k)*tau)*yu(k) + dyuP;
            
        end
        PPs = R*fric(qdP).*abs(qdP^3);
        PPu = (Cmfd + Bf*qdP)*qdP;
    end

    function Cf = Cf_(qPminus, yI, yPminus)
        Cf = 0*2*par.a(ID)*par.rho*tau ...
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