function par = parameters_PTOpipelineStudy(par,designCase,designSubCase,pLmodel)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters_PTOpipelineStudy.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 4/10/2021
%
% PURPOSE/DESCRIPTION:
% 
%
% FILE DEPENDENCY: NA
%
% UPDATES:
% 4/10/2021 - created.
% 4/12/2021 - Added functions to calculate wave speed for the case of
% entrained air for the simple MOC model and no entrained air for the DGCM
% model.
% 4/16/2021 - Changed determination of number of segments and MOC 
% time step to be based on the number of segments (n_seg) instead of time 
% step (dt_moc).
% 5/18/2021 - Added coefficients for the unsteady friction from (Schohl
%           1993). Added unsteady friction model indices to switch
%           statements.
% 05/21/2021 - Updated model number assignments.
% 06/01/2021 - Added switching parameters for pump and load.
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
    
%% Case independent parameters
        % fluid and entrianed gas properties
    par.rho = 1023; % [kg/m3] density of air
    par.mu = 9.4e-4; % [Pa-s]  Kinematic (absolute) viscosity
    par.beta = 2.2e9; % [Pa]  Bulk Modulus of air free fluid
    par.p_vap = 0.037e5; % [Pa] vapour pressure of seawater
    par.R = 0.0001; % [-] fraction  Baseline fraction of air by volume entrained in the hydraulic fluid at atm
    par.p_o = 101.3e3; % [Pa]  Atmospheric pressure (reference)
    par.gamma = 1.4; % [-]ratio of specific heats for air
    
    
     % WEC-driven pump
    par.irreg = 1;
    par.w1 = 0.1;
    par.wend = 10;
    par.dw = 0.005; % for dw = 0.01: repeat at t=628 s
                    % for dw = 0.005: repeat at t=1257 s
                    
     % Swithcing parameters **** will probably move to individual design
                                % cases below
    par.Tsw_pump = 10;
    par.duty_pump = 0.0;
    
    par.Tsw_load = 10;
    par.duty_load = 0.0;
    
     % pipeline
    par.pLmodel = pLmodel;    % selection of pipeline model
                        % 1 - Lumped resistance using laminar friction
                        % 2 - Lumped resistance including turbulent 
                        %     friction
                        % 3 - Single Pi-lump resistance and inertance with
                        %     low inertance to spoof model 2. Same results 
                        %     but faster.
                        % 4 - Single Pi-lump resistance and inertance
                        % 5 - N Pi-lumps, multi-segments, resistance,
                        %     inertance, and capacitance
                        % 6 - Method of characteristics based solution 
                        %     assuming constant wave speed
                        % 7 - Method of characteristics based Discrete 
                        %     free-gas cavity model (DGCM)
                        % 8 - Single Pi-lump resistance and inertance, with
                        %     unsteady friction
                        % 9 - N Pi-lumps, multi-segments, resistance,
                        %     inertance, and capacitance, with unsteady
                        %     friction
                        % 10- Method of characteristics based solution 
                        %     assuming constant wave speed, with unsteady
                        %     friction
                        % 11- Method of characteristics based Descrete 
                        %     free-gas cavity model (DGCM, with unsteady
                        %     friction

    par.n_lines = 2;
    par.phi = 0.75; % weighting coefficient for DGCM model
    
    % Unsteady friction model coefficients from (Schohl 1993)
    switch par.pLmodel
        case {8 9 10 11}
            par.m = [   1.051;...
                        2.358;...
                        9.021;...
                        29.47;...
                        79.55   ];

            par.n = [   26.65;...
                        100;...
                        669.6;...
                        6497;...
                        57990  ];

            par.nyf = length(par.m); % number of accumulated sums used for unsteady 
                                 % friciton modeling
            
        otherwise
            par.nyf = 0;
    end
    
%% Case dependent parameters
    switch designCase
        
%% CASE 1 Tp = 6s, 100kW, 1000 meter pipe, 
% various capacitance levels
        case 1 
           % WEC-driven pump, flow input parameters
            par.Tp = 6; % [s] peak period of waves
            par.Xq = .103; % [?] some kind of parameter for changing the magnitude of the flow rate power spectrum

           % load
            par.deltaP_mean = 6e6; % [Pa] nominal mean pressure at load
            q_mean = 1.235*par.Xq/par.Tp;
            par.R_load = par.deltaP_mean/q_mean; % [Pa.s/m^3] resistance of load
            
           % low-pressure inlet
%             par.p_tank = 2e6; % [Pa] pressure at the inlet to the LP pipeline

           % pipeline
            par.L_line = 1000; % [m] pipeline length
            par.d_line = 0.15; % [m] internal pipeline diameter
            
            if sum(par.pLmodel == [5 9]); n_seg_ = 6; end

           % low-pressure inlet & accumulators
            switch designSubCase
                case 1
                    par.p_tank = 1.8e6; % [Pa] pressure at the inlet to the LP pipeline
                    
                    par.C_hin = 5e-08; %0.5/6e6; % [m^3/Pa] capacitance of accumulator inlet to HP pipeline
                    par.C_hout = 5e-08; %0.5/6e6; % [m^3/Pa] capacitance of accumulator outlet to HP pipeline
                    par.C_lout = 10e-08; %0.5/6e6; % [m^3/Pa] capacitance of accumulator outlet to LP pipeline
                case 2
                    par.p_tank = 1.35e6; 
                    
                    par.C_hin = 10e-08;
                    par.C_hout = 10e-08;
                    par.C_lout = 20e-08;
                case 3
                    par.p_tank = 1.1e6; 
                    
                    par.C_hin = 20e-08;
                    par.C_hout = 20e-08;
                    par.C_lout = 40e-08;
                case 4
                    par.p_tank = .95e6; 
                    
                    par.C_hin = 40e-08;
                    par.C_hout = 40e-08;
                    par.C_lout = 80e-08;
                case 5
                    par.p_tank = 1.35e6; 
        
                    par.C_hin = 1e-08;
                    par.C_hout = 19e-08;
                    par.C_lout = 20e-08;
                case 6
                    par.p_tank = 1.35e6; 
                    
                    par.C_hin = 19e-08;
                    par.C_hout = 1e-08;
                    par.C_lout = 20e-08;
            end
%% CASE 2 Tp = 6s, 100kW, 1000 meter pipe, 
% Shorter pipeline length varying accumulator capacitance

        case 2
           % WEC-driven pump, flow input parameters
            par.Tp = 6; % [s] peak period of waves
            par.irreg = 1;
            par.Xq = 0.103; % [?] some kind of parameter for changing the magnitude of the flow rate power spectrum
            
           % load
            par.deltaP_mean = 6e6; % [Pa] nominal mean pressure at load
            q_mean = 1.235*par.Xq/par.Tp;
            par.R_load = par.deltaP_mean/q_mean; % [Pa.s/m^3] resistance of load

           % pipeline
            par.L_line = 100; % [m] pipeline length
            par.d_line = 0.10; % [m] internal pipeline diameter
            if sum(par.pLmodel == [5 9]); n_seg_ = 3; end

            switch designSubCase
                case 1
                    par.p_tank = 1.35e6; 
                    
                    par.C_hin = 10e-08;
                    par.C_hout = 10e-08;
                    par.C_lout = 20e-08;
                case 2
                    par.p_tank = 1.35e6; 
                    
                    par.C_hin = 1e-08;
                    par.C_hout = 19e-08;
                    par.C_lout = 20e-08;
                case 3
                    par.p_tank = 1.35e6; 
                    
                    par.C_hin = 19e-08;
                    par.C_hout = 1e-08;
                    par.C_lout = 20e-08;
            end
            
%% CASE 3 Entrained air fraction
        case 3
           % WEC-driven pump, flow input parameters
            par.Tp = 6; % [s] peak period of waves
            par.irreg = 1;
            par.Xq = 0.103; % [?] some kind of parameter for changing the magnitude of the flow rate power spectrum
            
           % load
            par.deltaP_mean = 6e6; % [Pa] nominal mean pressure at load
            q_mean = 1.235*par.Xq/par.Tp;
            par.R_load = par.deltaP_mean/q_mean; % [Pa.s/m^3] resistance of load

           % pipeline
            par.L_line = 1000; % [m] pipeline length
            par.d_line = 0.15; % [m] internal pipeline diameter

            if sum(par.pLmodel == [5 9]); n_seg_ = 6; end
            
           % low-pressure inlet           
            par.p_tank = 1.35e6;  % [Pa] pressure at the inlet to the LP pipeline
            
           % accumulators
            par.C_hin = 10e-08;
            par.C_hout = 10e-08;
            par.C_lout = 20e-08;
                    
            switch designSubCase
                case 1
                    par.R = 0.0001; % [-] fraction  Baseline fraction of air by volume entrained in the hydraulic fluid at atm
                case 2
                    par.R = 0.001; % [-] fraction  Baseline fraction of air by volume entrained in the hydraulic fluid at atm

            end
 
            
%% CASE 4 Tp = 6s, 100kW, Resonant pipeline lengths
        case 4 
           % WEC-driven pump, flow input parameters
%             par.Tp = 10; % [s] peak period of waves
            par.Xq = .103; % [?] some kind of parameter for changing the magnitude of the flow rate power spectrum

            switch designSubCase
                case 1
                    par.Tp = 6; % [s] peak period of waves
                    par.L_line = 2200; % [m] pipeline length
                    fl = 1466/(2*par.L_line)
                    fp = 2/par.Tp
                    fp/fl
                    
                    if sum(par.pLmodel == [5 9]); n_seg_ = 13; end
                    
                case 2 %%%!!!
                    par.Tp = 6; % [s] peak period of waves
                    par.L_line = 1650; % [m] pipeline length
                    fl = 1466/(2*par.L_line)
                    fp = 2/par.Tp
                    fp/fl
                    
                    if sum(par.pLmodel == [5 9]); n_seg_ = 10; end 
                    
                case 3
                    par.Tp = 6; % [s] peak period of waves
                    par.L_line = 1100; % [m] pipeline length
                    fl = 1466/(2*par.L_line)
                    fp = 2/par.Tp
                    fp/fl
                    
                    if sum(par.pLmodel == [5 9]); n_seg_ = 7; end
                    
                case 4
                    par.Tp = 6; % [s] peak period of waves
                    par.L_line = 3300; % [m] pipeline length
                    fl = 1466/(2*par.L_line)
                    fp = 2/par.Tp
                    fp/fl
                    
                    if sum(par.pLmodel == [5 9]); n_seg_ = 19; end
                    
                case 5
                    par.Tp = 6; % [s] peak period of waves
                    par.L_line = 825; % [m] pipeline length
                    fl = 1466/(2*par.L_line)
                    fp = 2/par.Tp
                    fp/fl
                    
                    if sum(par.pLmodel == [5 9]); n_seg_ = 5; end
                    
                case 6
                    par.Tp = 6; % [s] peak period of waves
                    par.L_line = 1000; % [m] pipeline length
                    fl = 1466/(2*par.L_line)
                    fp = 2/par.Tp
                    fp/fl
                    
                    if sum(par.pLmodel == [5 9]); n_seg_ = 6; end
                    
                case 7
                    par.Tp = 6; % [s] peak period of waves
                    par.L_line = 2200*2; % [m] pipeline length
                    fl = 1466/(2*par.L_line)
                    fp = 2/par.Tp
                    fp/fl
                    
                    if sum(par.pLmodel == [5 9]); n_seg_ = 26; end
                    
            end
            
           % load
            par.deltaP_mean = 6e6; % [Pa] nominal mean pressure at load
            q_mean = 1.235*par.Xq/par.Tp;
            par.R_load = par.deltaP_mean/q_mean; % [Pa.s/m^3] resistance of load

           % low-pressure inlet
            par.p_tank = 1.35e6; % [Pa] pressure at the inlet to the LP pipeline

           % pipeline
%             par.L_line = 2000; % [m] pipeline length
            par.d_line = 0.15; % [m] internal pipeline diameter

           % accumulators                    
            par.C_hin = 10e-08;
            par.C_hout = 10e-08;
            par.C_lout = 20e-08;
                    
                    
%% CASE 5 Switching at pump
        case 5

            par.duty_pump = 0.2;
            
           % WEC-driven pump, flow input parameters
            par.Tp = 6; % [s] peak period of waves
            par.Xq = .103/(1-par.duty_pump); % [?] some kind of parameter for changing the magnitude of the flow rate power spectrum
                    
           % load
            par.deltaP_mean = 6e6; % [Pa] nominal mean pressure at load
            q_mean = 1.235*(1-par.duty_pump)*par.Xq/par.Tp;
            par.R_load = par.deltaP_mean/q_mean; % [Pa.s/m^3] resistance of load     

           % pipeline
            par.L_line = 1000; % [m] pipeline length
            par.d_line = 0.15; % [m] internal pipeline diameter
            
            if sum(par.pLmodel == [5 9]); n_seg_ = 6; end

            switch designSubCase
                case 1
                    par.Tsw_pump = 0.01;

                    par.p_tank = 1.8e6; % [Pa] pressure at the inlet to the LP pipeline
                    
                    par.C_hin = 5e-08; %0.5/6e6; % [m^3/Pa] capacitance of accumulator inlet to HP pipeline
                    par.C_hout = 5e-08; %0.5/6e6; % [m^3/Pa] capacitance of accumulator outlet to HP pipeline
                    par.C_lout = 10e-08; %0.5/6e6; % [m^3/Pa] capacitance of accumulator outlet to LP pipeline

                case 2
                    par.Tsw_pump = 0.1;

                    par.p_tank = 1.35e6; % [Pa] pressure at the inlet to the LP pipeline
                    
                    par.C_hin = 10e-08; %0.5/6e6; % [m^3/Pa] capacitance of accumulator inlet to HP pipeline
                    par.C_hout = 10e-08; %0.5/6e6; % [m^3/Pa] capacitance of accumulator outlet to HP pipeline
                    par.C_lout = 20e-08; %0.5/6e6; % [m^3/Pa] capacitance of accumulator outlet to LP pipeline

                case 3
                    par.Tsw_pump = 0.1;

                    par.p_tank = 1.35e6; 
                    
                    par.C_hin = 1e-08;
                    par.C_hout = 19e-08;
                    par.C_lout = 20e-08;
                    
            end
            
%% CASE 6 Switching at load
        case 6
            
            par.duty_load = 0.2;
            
           % WEC-driven pump, flow input parameters
            par.Tp = 6; % [s] peak period of waves
            par.Xq = .103; % [?] some kind of parameter for changing the magnitude of the flow rate power spectrum
                   
           % load
            par.deltaP_mean = 6e6; % [Pa] nominal mean pressure at load
            q_mean = 1.235*par.Xq/par.Tp;
            par.R_load = (1-par.duty_load)*par.deltaP_mean/q_mean; % [Pa.s/m^3] resistance of load
            
           % pipeline
            par.L_line = 1000; % [m] pipeline length
            par.d_line = 0.15; % [m] internal pipeline diameter
            
            if sum(par.pLmodel == [5 9]); n_seg_ = 6; end

            switch designSubCase
                case 1
                    
                    par.Tsw_load = 0.1;

                    par.p_tank = 1.8e6; % [Pa] pressure at the inlet to the LP pipeline
                    
                    par.C_hin = 5e-08; %0.5/6e6; % [m^3/Pa] capacitance of accumulator inlet to HP pipeline
                    par.C_hout = 5e-08; %0.5/6e6; % [m^3/Pa] capacitance of accumulator outlet to HP pipeline
                    par.C_lout = 10e-08; %0.5/6e6; % [m^3/Pa] capacitance of accumulator outlet to LP pipeline
                
                case 2
                    
                    par.Tsw_load = 0.1;

                    par.p_tank = 1.35e6; % [Pa] pressure at the inlet to the LP pipeline
                    
                    par.C_hin = 10e-08; %0.5/6e6; % [m^3/Pa] capacitance of accumulator inlet to HP pipeline
                    par.C_hout = 10e-08; %0.5/6e6; % [m^3/Pa] capacitance of accumulator outlet to HP pipeline
                    par.C_lout = 20e-08; %0.5/6e6; % [m^3/Pa] capacitance of accumulator outlet to LP pipeline
                
            end
            
%% CASE 99 Regular waves
        case 99 
           % WEC-driven pump, flow input parameters
%             par.Tp = 10; % [s] peak period of waves
            par.irreg = 0;
            par.Xq = .05; % [?] some kind of parameter for changing the magnitude of the flow rate power spectrum
            
            switch designSubCase
                case 1
                    par.Tp = 6; % [s] peak period of waves
                    par.L_line = 2190; % [m] pipeline length
                    fl = 1460/(2*par.L_line)
                    fp = 2/par.Tp
                    fp/fl
                case 2
                    par.Tp = 6; % [s] peak period of waves
                    par.L_line = 1095; % [m] pipeline length
                    fl = 1460/(2*par.L_line)
                    fp = 2/par.Tp
                    fp/fl
                case 3
                    par.Tp = 6; % [s] peak period of waves
                    par.L_line = 3285; % [m] pipeline length
                    fl = 1460/(2*par.L_line)
                    fp = 2/par.Tp
                    fp/fl
                case 4
                    par.Tp = 6; % [s] peak period of waves
                    par.L_line = 730; % [m] pipeline length
                    fl = 1460/(2*par.L_line)
                    fp = 2/par.Tp
                    fp/fl
                case 5
                    par.Tp = 6; % [s] peak period of waves
                    par.L_line = 2*730; % [m] pipeline length
                    fl = 1460/(2*par.L_line)
                    fp = 2/par.Tp
                    fp/fl
                case 6
                    par.Tp = 6; % [s] peak period of waves
                    par.L_line = 2190+730; % [m] pipeline length
                    fl = 1460/(2*par.L_line)
                    fp = 2/par.Tp
                    fp/fl
            end
            
           % load
            par.deltaP_mean = 6e6; % [Pa] nominal mean pressure at load
            q_mean = 5*par.Xq/par.Tp;
            par.R_load = par.deltaP_mean/q_mean; % [Pa.s/m^3] resistance of load

           % low-pressure inlet
            par.p_tank = 2e6; % [Pa] pressure at the inlet to the LP pipeline

           % pipeline
%             par.L_line = 2000; % [m] pipeline length
            par.d_line = 0.15; % [m] internal pipeline diameter

           % accumulators
            par.p_tank = 1.6e6; 

            par.C_hin = 5e-08;
            par.C_hout = 10e-08;
            par.C_lout = 20e-08;
                    
                    
                    % wavelength/pipeline length = a*t_p/2/L
                    % 1460*(12)/2/1000

            


    end
    
%% Case dependent parameters but case independent code

    
     % pipeline
    par.A_line = pi/4*par.d_line^2; % [m^2] pipeline cross-sectional flow area
                 
    switch par.pLmodel
        case {5 9}
            n_seg = n_seg_;
        case {6 10}
%             % MOC solver time step (total timestep)
%             par.dt_MOC(1) = 0.01; % [s] LP line
%             par.dt_MOC(2) = 0.01; % [s] HP line

            % wave speed in each line
%             par.a(1) = 1200; % [m/s] LP wave speed 
%             par.a(2) = 1200; % [m/s] HP wave speed
     
            p_refLP = par.p_tank;
            par.a(1) = waveSpeedEntAir(p_refLP,par); % [m/s] LP wave speed 
            
            p_refHP = par.deltaP_mean;
            par.a(2) = waveSpeedEntAir(p_refHP,par); % [m/s] HP wave speed    
            
            switch designCase
                case 4 % CASE 4 Tp = 6s, 100kW, Resonant pipeline lengths
                    n_seg = 100;
                case 2 % CASE 2 Tp = 6s, 100kW, 1000 meter pipe, 
                        % Shorter pipeline length varying accumulator capacitance
                    switch designSubCase
                        case 2
                            n_seg = 50;
                        otherwise
                            n_seg = 10;
                    end
                otherwise
                    n_seg = 50;
            end  
          
        case {7 11}    
            % wave speed in each line
            par.a(1) = waveSpeed(par); % [m/s] LP wave speed 
            par.a(2) = waveSpeed(par); % [m/s] HP wave speed   
            
            % MOC solver time step (total timestep)
%             par.dt_MOC(1) = 0.01; % [s] LP line
%             par.dt_MOC(2) = 0.01; % [s] HP line            

            switch designCase
                case 2
                    n_seg = 100;
                case 4
                    switch designSubCase
                        case 2
                            n_seg = 50;
                        otherwise
                            n_seg = 10;
                    end
                otherwise
                    n_seg = 50;
            end
            
        otherwise
            n_seg = [];
    end
    
    if sum(par.pLmodel == [8 9 10 11])
        par.T = par.rho*par.d_line^2/4/par.mu; % unsteady friction time scale
    end
    
    par = pLmodelParamSetup(par,n_seg);

    %% %%%%%%%%%%%%   FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function a = waveSpeed(par)
         % wave speed of fluid in thick walled tube
    	a = sqrt(par.beta/par.rho);
end

function a = waveSpeedEntAir(p,par)
    	 % isothermal bulk modulus
        beta_eff = par.beta/(1 + par.beta*(par.R*par.p_o/p^2));
         % wave speed of fluid in thick walled tube
    	a = sqrt(beta_eff/par.rho);
    end

end