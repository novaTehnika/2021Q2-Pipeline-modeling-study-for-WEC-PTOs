function out = sim_PTOpipeline(par,getTDdata)
% sim_PTOpipeline.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 01/18/2021
%
% PURPOSE/DESCRIPTION:
% This script simulates a wave energy PTO with a long pipeline in the
% high-pressure branch and low-pressure branch. The flow through the
% WEC-driven pump is the input to the system and is determined from a power
% spectrum representing waves. The dynamics of the WEC are ignored and the
% the power spectrum is used directly to calculate the flow rate of the
% WEC-driven pump. A variety of pipeline models are implemented for
% comparison. Accumulators are modeled as linear capacitive elements. The
% low-pressure source node, at the inlet of the low-pressure pipeline, is
% modeled as ... The load on the system is modeled as a linear resistance;
% the resistance of the load is calculated from the mean flow of the
% WEC-driven pump (a predetermined value) and a specified nominal mean
% load pressure. (OR maybe just a specified load resistance)
%
% FILE DEPENDENCY:
% sys_PTOpipeline.m
% nonStateVarsPTOwPL
%   pumpFlow.m
% pipeline
%   pipelineRlam.m
%   pipelineRturb.m
%   pipelinePi.m
%   pipelineNPi.m
%   pipelineMOC.m
%     mocSimple.m
%     mocDGCM.m
%     mocSimple_uFric.m
%     mocDGCM_uFric.m
%   pLmodelParamSetup
% statsTimeVar
% statsTimeVar_cdf
% 
%
% UPDATES:
% 01/18/2021 - First version creation.
% 01/31/2021 - Implemented short line model with both laminar 
%           and turbulent friction using fsolve to acount for flow rate 
%           dependent Re. 
% 01/31/2021 - Rewrite staticVars() to select pipeline model.
%           Write seperate pipeline models for use by staticVars() instead
%           of using pipeline() as in previous versions.
% 02/05/2021 - Impliment time step checking for MOC model to
%           handle ODE solver steping back in time. Requires ODE solver to
%           have maximum time step.
% 02/31/2021 - Added irregular WEC motions with 
%           Pierson-Moskowitz spectrum and a normalizing factor X.
% 03/12/2021 - Changed mocSimple() model to deal with flow rate 
%           instead of mean flow velocity. Streamlined code in mocSimple by
%           including Re calc in fric() and implimenting inertance calc, 
%           I = a*rho/A.
% 03/12/2021 - Added discrete free-gas cavity model (DGCM) model as 
%           pipeline model 6.
% 03/17/2021 -  Added calculation of power loss in MOC pipeline 
%           models and in post-process.
% 03/17/2021 - Condensed friction term calculations using R = rho*dx/2dA^2
% 03/23/2021 - Added mass balance with calculation of the change in volume
%           of the entrained gas in the pipelines. Also added an effective
%           bulk modulus based on isothermal compresion to the capacitacne 
%           calculation for the lumped parameter model.
% 03/24/2021 - Some cleaning-up of code. Made mass and power balance work
%           for every model.
% 03/25/2021 - Removed pipeline functions to seperate m-files.
%           Seperated pipeline model switching and pipeline model/managers.
%           Added list of file depedancies. Some cleaning up of code and
%           changes to variable names. Added postprocessing input to sys().
% 03/29/2021 - Converted script to function for use in studies.
% 03/29/2021 - Added ODE solver error parameters to "par" struc.
%          	(an input to the function).
% 04/2/2021 - Fixed pipeline time domain solutions as output.
% 05/18/2021 - Updated model versions implimenting unsteady 
%           friction. Added unsteady friction model indices to switch
%           statements.
% 05/21/2021 - Updated model number assignments. Added initial conditions 
%           for accumulated sums in the unsteady friciton model. Added 
%           outputs for steady and unsteady friction power loss.
% 05/25/2021 - added conditional state for assignment of itemized pipeline 
%           losses (PPs and PPu).
% 06/01/2021 - Added logic to set ode solver MaxStep to
%           fraction of smalled swithing period.
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
clear pumpFlow pipeline pipelineMOC

%% %%%%%%%%%%%%   SOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solve for states of the dynamic system
 % Set-up solver
    tspan = [0 par.tend];   % time interval
    switch par.pLmodel  % inital conditions
        case {1 2 6 7 10 11}
            y0 = [  par.deltaP_mean; ...
                    par.deltaP_mean; ...
                    par.p_tank];
        case {3 4 5 8 9}
            y0 = [  par.deltaP_mean; ...
                    par.deltaP_mean; ...
                    par.p_tank; ...
                    par.p_tank*mod(2:2*par.n_seg(1),2)'; ...
                    zeros(1,par.nyf*par.n_seg(1))'; ...
                    par.deltaP_mean*mod(2:2*par.n_seg(2),2)'; ...
                    zeros(1,par.nyf*par.n_seg(2))'];
    end

 % Solver options
    options = odeset('RelTol',par.odeSolverRelTol,...
                     'AbsTol',par.odeSolverAbsTol); 
    switch par.pLmodel
        case {6 7 10 11}
            options.MaxStep = min([par.dt_MOC(:)' ...
                                0.1*[par.Tsw_pump par.Tsw_load]]);
        otherwise
            options.MaxStep = min(0.1*[par.Tsw_pump par.Tsw_load]);
    end

 % Run solver
    [t,y] = ode15s(@(t,y) sys_PTOpipeline(t,y,par,0), tspan, y0, options);

%% %%%%%%%%%%%%   POST-PROCESS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % determine number of time steps
    nt = length(t); % number of time steps

    % Extract system states (non-pipeline states/variables)
    p_hin = y(:,1); % [Pa] pressure at outlet of HP pipeline
    p_hout = y(:,2); % [Pa] pressure at inlet of HP pipeline
    p_lout = y(:,3); % [Pa] pressure at outlet of LP pipeline

    % Post process non-state, time and state-dependent variables
     % Initialize arrays
    q_p = zeros(nt,1); % [m^3/s] pump flow rate
    q_L = zeros(nt,1); % [m^3/s] load flow rate
    p_lin = zeros(nt,1); % [Pa] tank pressure
    dp_hindt = zeros(nt,1); % time deriv. of press. at outlet of HP pipeline
    dp_houtdt = zeros(nt,1); % time deriv. of press. at inlet of HP pipeline
    dp_loutdt = zeros(nt,1); % time deriv. of press. at outlet of LP pipeline

     % iterate through time
    for it = 1:nt
        nonState = nonStateVarsPTOwPL(t(it),y(it,:),par); % recover nonstate variables
        q_p(it) = nonState.q_p;
        q_L(it) = nonState.q_L;
        p_lin(it) = nonState.p_lin;

        dydt = sys_PTOpipeline(t(it),y(it,:)',par,1);
        dp_hindt(it) = dydt(1);
        dp_houtdt(it) = dydt(2);
        dp_loutdt(it) = dydt(3);
    end

    %% Extract pipeline states/variables
      % Initialize arrays
    q_lin = zeros(nt,1);
    q_lout = zeros(nt,1);
    q_hin = zeros(nt,1);
    q_hout = zeros(nt,1); 

    qLP = zeros(nt,par.n_seg(1)+1);
    pLP = zeros(nt,par.n_seg(1)+1);
    PPLP_sfriction = zeros(nt,1);
    PPLP_ufriction = zeros(nt,1);
    PPLP_friction = zeros(nt,1);

    qHP = zeros(nt,par.n_seg(2)+1);
    pHP = zeros(nt,par.n_seg(2)+1);
    PPHP_sfriction = zeros(nt,1);
    PPHP_ufriction = zeros(nt,1);
    PPHP_friction = zeros(nt,1);


    switch par.pLmodel
        case {1 2}
            for it = 1:nt
                % Get states of LP pipeline (in post process mode)
                [q_lin(it), q_lout(it),~,pLsoln] = ...
                    pipeline(t(it),[],p_lin(it),p_lout(it),par,1,1);
                PPLP_friction(it) = pLsoln.PPfric;
                % Get states of HP pipeline (in post process mode)
                [q_hin(it), q_hout(it),~,pLsoln] = ...
                    pipeline(t(it),[],p_hin(it),p_hout(it),par,2,1);
                PPHP_friction(it) = pLsoln.PPfric;
            end
        case {3 4 5 8 9}
            % LP pipeline states
            iqLP = 1:2:2*par.n_seg(1)-1;
            ipLP = 2:2:2*par.n_seg(1)-1;
            yLP = y(:,par.iy_LP);
            qLP = yLP(:,iqLP);
            pLP = yLP(:,ipLP);
            q_lin = qLP(:,1);
            q_lout = qLP(:,end);
            for it = 1:nt
                [~,~,~,pLsoln] = ...
                        pipeline([],yLP(it,:)',p_lin(it),p_lout(it),par,1,1);
                if sum(par.pLmodel == [8 9 10 11])
                    PPLP_sfriction(it) = pLsoln.PPs;
                    PPLP_ufriction(it) = pLsoln.PPu;
                end
                PPLP_friction(it) = pLsoln.PPfric;
            end
            % HP pipeline states
            iqHP = 1:2:2*par.n_seg(2)-1;
            ipHP = 2:2:2*par.n_seg(2)-1; 
            yHP = y(:,par.iy_HP);
            qHP = yHP(:,iqHP);
            pHP = yHP(:,ipHP);
            q_hin = qHP(:,1);
            q_hout = qHP(:,end);
            for it = 1:nt
                [~,~,~,pLsoln] = ...
                        pipeline([],yHP(it,:)',p_hin(it),p_hout(it),par,2,1);
                if sum(par.pLmodel == [8 9 10 11])
                    PPHP_sfriction(it) = pLsoln.PPs;
                    PPHP_ufriction(it) = pLsoln.PPu;
                end
                PPHP_friction(it) = pLsoln.PPfric;
            end
        case {6 7 10 11}
           for it = 1:nt
                % Get states of LP pipeline (in post process mode)
                [q_lin(it), q_lout(it),~,pLsoln] = ...
                    pipeline(t(it),[],[],[],par,1,1);
                qLP(it,:) = pLsoln.q;
                pLP(it,:) = pLsoln.p;
                if sum(par.pLmodel == [8 9 10 11])
                    PPLP_sfriction(it) = pLsoln.PPs;
                    PPLP_ufriction(it) = pLsoln.PPu;
                end
                PPLP_friction(it) = pLsoln.PPfric;
                % Get states of HP pipeline (in post process mode)
                [q_hin(it), q_hout(it),~,pLsoln] = ...
                    pipeline(t(it),[],[],[],par,2,1);
                qHP(it,:) = pLsoln.q;
                pHP(it,:) = pLsoln.p;
                if sum(par.pLmodel == [8 9 10 11])
                    PPHP_sfriction(it) = pLsoln.PPs;
                    PPHP_ufriction(it) = pLsoln.PPu;
                end
                PPHP_friction(it) = pLsoln.PPfric;
           end
    end

    %% Energy Balance
    % Calculate power exerted on pipelines (PPin - PPout)
    PPLP_total = p_lin.*q_lin - p_lout.*q_lout;
    PPHP_total = p_hin.*q_hin - p_hout.*q_hout;

    % Calculate energy loss
     % calc. kinetic energy in pipelines at end of simulation
    I_LP = par.rho*par.L_line/(par.n_seg(1)*par.A_line);
    I_HP = par.rho*par.L_line/(par.n_seg(2)*par.A_line);
    KELP_start = sum(0.5*I_LP*qLP(1,:).^2);
    KEHP_start = sum(0.5*I_HP*qHP(1,:).^2);
    KELP_end = sum(0.5*I_LP*qLP(end,:).^2);
    KEHP_end = sum(0.5*I_HP*qHP(end,:).^2);
     % calc. mean energy loss
    PPmean_LP_fric = trapz(t,PPLP_friction)/par.tend;
    PPmean_HP_fric = trapz(t,PPHP_friction)/par.tend;
      % calc. mean steady and unsteady friction seperately
    if sum(par.pLmodel == [8 9 10 11])
        PPmean_LP_sfric = trapz(t,PPLP_sfriction)/par.tend;
        PPmean_LP_ufric = trapz(t,PPLP_ufriction)/par.tend;
        PPmean_HP_sfric = trapz(t,PPHP_sfriction)/par.tend;  
        PPmean_HP_ufric = trapz(t,PPHP_ufriction)/par.tend;
    end

    PPmean_LP_total = (trapz(t,PPLP_total) - (KELP_end - KELP_start))/par.tend;
    PPmean_HP_total = (trapz(t,PPHP_total) - (KEHP_end - KEHP_start))/par.tend;

    PPbalLP = PPmean_LP_total - PPmean_LP_fric;
    PPbalHP = PPmean_HP_total - PPmean_HP_fric;
    
    %% Mass Balance
    % Calculate the voluem stored in the pipeline
     % low-pressure pipeline
    if par.n_seg(1) > 1
        if par.pLmodel == 5 || par.pLmodel == 9 % lumped parameter model;
                            % all pressure nodes are compressible
            iip = 1:par.n_seg(1)-1;
        else                % MOC model; external nodes are not compressible
            iip = 2:par.n_seg(1)-2;
        end
        VgLP_start = par.R*par.p_o*par.A_line*par.L_line/par.n_seg(1)*sum(1./pLP(1,iip));
        VgLP_end = par.R*par.p_o*par.A_line*par.L_line/par.n_seg(1)*sum(1./pLP(end,iip));
        deltaVgLP = VgLP_end - VgLP_start;
    else % no compressible volumes in the pipeline model
        VgLP_start = 0;
        VgLP_end = 0;
        deltaVgLP = VgLP_end - VgLP_start;
    end
     % high-pressure pipeline
    if par.n_seg(2) > 1
        if par.pLmodel == 5 || par.pLmodel == 9 % lumped parameter model;
                            % all pressure nodes are compressible
            iip = 1:par.n_seg(2)-1;
        else                % MOC model; external nodes are not compressible
            iip = 2:par.n_seg(2)-2;
        end
        VgHP_start = par.R*par.p_o*par.A_line*par.L_line/par.n_seg(2)*sum(1./pHP(1,iip));
        VgHP_end = par.R*par.p_o*par.A_line*par.L_line/par.n_seg(2)*sum(1./pHP(end,iip));
        deltaVgHP = VgHP_end - VgHP_start;
    else % no compressible volumes in the pipeline model
        VgHP_start = 0;
        VgHP_end = 0;
        deltaVgHP = VgHP_end - VgHP_start;
    end

    % Calculate the net volume flow into the pipeline
    VnetLP = trapz(t,q_lin - q_lout)/par.tend;
    VnetHP = trapz(t,q_hin - q_hout)/par.tend;
    % Total mass balance
    VbalLP = VnetLP - deltaVgLP;
    VbalHP = VnetHP - deltaVgHP;

    %% Pressure and rate of change in pressure metrics
    p_loutDist = statsTimeVar(t,p_lout);
    p_hinDist = statsTimeVar(t,p_hin);
    p_houtDist = statsTimeVar(t,p_hout);
    dp_loutdtDist = statsTimeVar_cdf(t,abs(dp_loutdt));
    dp_hindtDist = statsTimeVar_cdf(t,abs(dp_hindt));
    dp_houtdtDist = statsTimeVar_cdf(t,abs(dp_houtdt));
    
    %% Torque
    deltap_wec = p_hout-p_lout;

    deltap_wecDist = statsTimeVar(t,deltap_wec);

%% %%%%%%%%%%%%   FUNCTION OUTPUTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out.PPmean_LP_fric = PPmean_LP_fric;
out.PPmean_HP_fric = PPmean_HP_fric;
if sum(par.pLmodel == [8 9 10 11])
    out.PPmean_LP_sfric = PPmean_LP_sfric;
    out.PPmean_LP_ufric = PPmean_LP_ufric;
    out.PPmean_HP_sfric = PPmean_HP_sfric;
    out.PPmean_HP_ufric = PPmean_HP_ufric;
else
    out.PPmean_LP_sfric = PPmean_LP_fric;
    out.PPmean_LP_ufric = 0;
    out.PPmean_HP_sfric = PPmean_HP_fric;
    out.PPmean_HP_ufric = 0;
end
out.PPbalLP = PPbalLP;
out.PPbalHP = PPbalHP;
out.VbalLP = VbalLP;
out.VbalHP = VbalHP;
out.p_loutDist = p_loutDist;
out.p_hinDist = p_hinDist;
out.p_houtDist = p_houtDist;
out.dp_loutdtDist = dp_loutdtDist;
out.dp_hindtDist = dp_hindtDist;
out.dp_houtdtDist = dp_houtdtDist;
out.deltap_wecDist = deltap_wecDist;
out.par = par;

out.t = t;
out.qLP = qLP;
out.pLP = pLP;
out.qHP = qHP;
out.pHP = pHP;

if getTDdata
    out.t = t;
    out.p_hin = p_hin;
    out.p_hout = p_hout;
    out.p_lout = p_lout;
    out.q_p = q_p;
    out.q_L = q_L;
    out.p_lin = p_lin;
    out.dp_hindt = dp_hindt;
    out.dp_houtdt = dp_houtdt;
    out.dp_loutdt = dp_loutdt;
    out.q_lin = q_lin;
    out.q_lout = q_lout;
    out.q_hin = q_hin;
    out.q_hout = q_hout;
    out.qLP = qLP;
    out.pLP = pLP;
    out.PPLP_friction = PPLP_friction;
    out.qHP = qHP;
    out.pHP = pHP;
    out.PPHP_friction = PPHP_friction;
    out.deltap_wec = deltap_wec;
    
end

end
