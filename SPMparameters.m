%% Initialize Single Particle Model Parameterss (SI Units)
% Defining a structure named "SPMparams" that contains all parameter values

SPMparams.Rg = 8.314;                 % Gas constant
SPMparams.F = 96487;                  % Faraday's constant
SPMparams.alpha_cell = 0.5;           % Anode/Cathode transfer coefficient
SPMparams.ce0 = 1200;                 % Average electrolyte concentration [mol/m^3]
SPMparams.capacity = 2;               % Noimnal cell capacity (Ah)   

SPMparams.c_n_max = 27265;            % Maximum anode concentration
SPMparams.c_p_max = 49389;            % Maximum cathode concentration
SPMparams.Ds_n = 2.87e-14;            % Anode diffusion coefficient
SPMparams.Ds_p = 4.85e-14;            % Cathode diffusion coefficient
SPMparams.Rs_n = 5e-6;                % Anode particle radius
SPMparams.Rs_p = 5e-6;                % Cathode particle radius
SPMparams.A = 0.077;                  % Area
SPMparams.Ln = 67e-6;                 % Anode thickness
SPMparams.Lp = 52e-6;                 % Cathode thickness
SPMparams.epsilon_n = 0.565;          % Anode solid phase volume fraction
SPMparams.epsilon_p = 0.583;          % Cathode solid phase volume fraction
SPMparams.theta100_n = 0.928;         % Anode stoichiometry at 100% SOC
SPMparams.theta100_p = 0.3504;        % Cathode stoichiometry at 100% SOC
SPMparams.k_n = 3.48e-10;             % Anode reaction rate constant
SPMparams.k_p = 4.164e-10;            % Cathode reaction rate constant
SPMparams.R_c = 0.029;                % Lumped resistance
SPMparams.theta0_n = 0.002;           % Anode stoichiometry at 0% SOC
SPMparams.theta0_p = 0.9986;          % Cathode stoichiometry at 0% SOC

% Specific interfacial (electroactive) surface area for anode and cathode
SPMparams.a_s_n = 3*SPMparams.epsilon_n/SPMparams.Rs_n;   
SPMparams.a_s_p = 3*SPMparams.epsilon_p/SPMparams.Rs_p; 

SPMparams.T_batt = 296.15;  