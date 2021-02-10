%% Energy 295 Spring 2019 - Homework 4
%  Single Particle Model (SPM) implementation using Finite Difference
%  Method (FDM).  
% 
%  Plots for model voltage, and bulk SOC at both electrodes is included.
%
%  Supporting files:
%  SPMparameters.m  : Values of all model parameters
%  ode_spm.m        : ODE system formulation
%  U_p.m            : Open circuit potential for cathode
%  U_n.m            : Open circuit potential for anode
%  soc_bound.m      : Saturate SOC values at its physical limits of 0 or 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
close all; 
clc; 
warning off

%% Experimental Data
% Load current data, define sampling time and time vector
% Variable for current data is: I_data
% I_data sign convention: Positive for discharge
% Load experimental data below (Comment or uncomment the desired data set)

% % 1C, 2C, 5C: Experimental current and voltage data
%load NMC_Cell_H1_T23_1C_CTID.mat; dt = 1;
% load NMC_Cell_H1_T23_2C_CTID.mat; dt = 1;
% load NMC_Cell_H1_T23_5C_CTID.mat; dt = 1;

% % UDDS: Experimental current and voltage data
 load US06_test.mat; dt = 0.1; 
 %I=[ones(3000,1).*4; ones(3000,1).*2; ones(3000,1).*-1; ones(3000,1).*2;ones(3000,1).*3; ones(3000,1).*-1.5;ones(3000,1).*2;ones(3000,1).*1;ones(3000,1).*0.5;ones(3000,1).*1.5; ones(3000,1).*-3; ones(3000,1).*1; ones(3000,1).*4; ones(3000,1).*1.5;ones(3000,1).*-3;ones(3000,1).*1.5;ones(3000,1).*4;ones(3000,1).*-2;ones(3000,1).*1.5 ];
%I=[ones(2000,1).*-3; ones(2000,1).*-2; ones(2000,1).*-1.5; ones(2000,1).*-1; ones(2000,1).*3; ones(2000,1).*0.5; ones(2000,1).*-3; ones(2000,1).*1.5; ones(2000,1).*-0.5; ones(2000,1).*-0.25; ones(2000,1).*-0.005; ones(2000,1).*4; ones(2000,1).*1; ones(2000,1).*-0.005; ones(2000,1).*-3; ones(2000,1).*-2; ones(2000,1).*-1; ones(2000,1).*0.5; ones(2000,1).*-3; ones(2000,1).*-4; ones(2000,1).*-0.5; ones(2000,1).*-3; ones(2000,1).*1.5];
 soc_init=0.55;
%I_data=I_expt;
dt=0.1;
 %load ('training_expt_data.mat');
% HPPC23=xlsread('NMC_Cell_H4_T23_HPPC_Test.xlsx');
% HPPC23_H4=HPPC23(5000:end,:);
%  dt=1;
%  t_expt=HPPC23_H4(:,2);
%  I_data=-1.*HPPC23_H4(:,3);
%  V_data=HPPC23_H4(:,4);
% 
%  I_data=-1.*US0623_H1(:,3); dt=0.1;
%  V_data=US0623_H1(:,4);

t_data = 0:dt:(length(I_data)-1)*dt;    % Load Time Vector
T_batt = 296.15;                           % Battery temperature [K]

%% Model Parameters 

run SPMparameters                       % Load SPM parameter values

%% Solid phase: Finite Difference Method Discretization

SPMparams.r_grid = 30;                  % Radial discretization grid points 

% Grid interval size for anode and cathode
SPMparams.delta_x_n = SPMparams.Rs_n/(SPMparams.r_grid-1);    
SPMparams.delta_x_p = SPMparams.Rs_p/(SPMparams.r_grid-1);  

% Coefficientes of discretized ODEs
SPMparams.alpha_n = SPMparams.Ds_n/(SPMparams.delta_x_n^2);
SPMparams.alpha_p = SPMparams.Ds_p/(SPMparams.delta_x_p^2);
SPMparams.beta_n = 1/(SPMparams.F*SPMparams.A*SPMparams.Ln*...
    SPMparams.a_s_n*SPMparams.delta_x_n);
SPMparams.beta_p = 1/(SPMparams.F*SPMparams.A*SPMparams.Lp*...
    SPMparams.a_s_p*SPMparams.delta_x_p);

% ODE: State-space formulation
A_mat = zeros(SPMparams.r_grid-1);
B_mat = zeros(SPMparams.r_grid-1,1);
for k = 1:SPMparams.r_grid-1
    if k == 1
        A_mat(k,k) = -2;
        A_mat(k,k+1) = 2;
        B_mat(k,1) = 0;
    elseif k == SPMparams.r_grid-1
        A_mat(k,k-1) = 2;
        A_mat(k,k) = -2;
        B_mat(k,1) = 2*(1+1/k);
    else
        A_mat(k,k-1) = 1*(1-1/k);
        A_mat(k,k) = -2;
        A_mat(k,k+1) = 1*(1+1/k);
        B_mat(k,1) = 0;
    end
end
SPMparams.A_mat_n = A_mat.*SPMparams.alpha_n;
SPMparams.B_mat_n = -B_mat.*SPMparams.beta_n;   
SPMparams.A_mat_p = A_mat.*SPMparams.alpha_p;
SPMparams.B_mat_p = B_mat.*SPMparams.beta_p;  

%% Initialize Concentration States

% Initialize concentration states in both electrodes
for i=1:SPMparams.r_grid
    cs_n0(i,1)= ((SPMparams.theta100_n-SPMparams.theta0_n)*(soc_init)+SPMparams.theta0_n)*SPMparams.c_n_max;
    cs_p0(i,1)= ((SPMparams.theta100_p-SPMparams.theta0_p)*(soc_init)+SPMparams.theta0_p)*SPMparams.c_p_max;
end

% Initial concentration for ODE solver (states: 2*(r_grid-1))
cs_0 = [cs_n0(2:SPMparams.r_grid,1);cs_p0(2:SPMparams.r_grid,1)];

%% Solve concentration ODEs

tspan = t_data;
reltol=1.0e-04; abstol=1.0e-04;
options=odeset('RelTol',reltol,'AbsTol',abstol,'MaxStep',5);
[t_out,cs] = ode23s(@(t_out,cs) ode_spm(t_out,cs,SPMparams,I_data,t_data), tspan, cs_0, options);
cs = cs';

%% Cell Voltage Calculation 

% Separate concentration as per the electrode
% Add concentration at the first grid where r=0 (not required, hence neglected)
cs_n = cs(1:SPMparams.r_grid-1,:);
cs_p = cs(SPMparams.r_grid:2*(SPMparams.r_grid-1),:);

% Surface concentration -> surface stoichiometry
theta_surf_n = cs_n(SPMparams.r_grid-1,:)/SPMparams.c_n_max;
theta_surf_p = cs_p(SPMparams.r_grid-1,:)/SPMparams.c_p_max;

% Open Circuit Potential (OCP) calculation
ocp_p = U_p(theta_surf_p);
ocp_n = U_n(theta_surf_n);

for i = 1:length(t_out)

% Exchange current density calculation (as a function of surface
% concentration)
i0_n(i) = SPMparams.k_n*SPMparams.F*(SPMparams.ce0^0.5)*...
    ((SPMparams.c_n_max*theta_surf_n(i))^0.5)*...
    ((SPMparams.c_n_max-SPMparams.c_n_max*theta_surf_n(i))^0.5);
i0_p(i) = SPMparams.k_p*SPMparams.F*(SPMparams.ce0^0.5)*...
    ((SPMparams.c_p_max*theta_surf_p(i))^0.5)*...
    ((SPMparams.c_p_max-SPMparams.c_p_max*theta_surf_p(i))^0.5);    

% Overpotential calculation
eta_n(i) = asinh(I_data(i)/(2*SPMparams.A*SPMparams.a_s_n*...
    SPMparams.Ln*i0_n(i)))*(SPMparams.Rg*T_batt)/...
    (SPMparams.F*SPMparams.alpha_cell);
eta_p(i) = asinh(I_data(i)/(2*SPMparams.A*SPMparams.a_s_n*...
    SPMparams.Ln*i0_p(i)))*(SPMparams.Rg*T_batt)/...
    (SPMparams.F*SPMparams.alpha_cell);

% Cell Voltage
V_cell(i) = ocp_p(i)-ocp_n(i)+eta_p(i)-eta_n(i)-I_data(i).*SPMparams.R_c;

end

%% SOC Calculation

% Reference SOC (Coulomb Counting)
SOC_cc = zeros(length(I_data),1);
SOC_cc(1) = soc_init;
for k=2:length(I_data)
    SOC_cc(k) = SOC_cc(k-1) - I_data(k)*dt/(SPMparams.capacity*3600);
    SOC_cc(k) = (SOC_cc(k));
end

% Bulk SOC
interval_n = [SPMparams.delta_x_n:SPMparams.delta_x_n:SPMparams.delta_x_n*(SPMparams.r_grid-1)]';
interval_p = [SPMparams.delta_x_p:SPMparams.delta_x_p:SPMparams.delta_x_p*(SPMparams.r_grid-1)]';
for i = 1:length(t_out)
cs_bulk_n(i) = trapz(interval_n,interval_n.^2.*cs_n(1:end,i))*3/(SPMparams.Rs_n^3);
cs_bulk_p(i) = trapz(interval_p,interval_p.^2.*cs_p(1:end,i))*3/(SPMparams.Rs_p^3);
soc_bulk_n(i) = (cs_bulk_n(i)/SPMparams.c_n_max - SPMparams.theta0_n)/(SPMparams.theta100_n - SPMparams.theta0_n);
soc_bulk_p(i) = (SPMparams.theta0_p - cs_bulk_p(i)/SPMparams.c_p_max)/(SPMparams.theta0_p - SPMparams.theta100_p); 

soc_bulk_n(i) = (soc_bulk_n(i));
soc_bulk_p(i) = (soc_bulk_p(i));

end

figure()
plot(t_data,SOC_cc,'b','LineWidth',2);hold on; grid on;
plot(t_out,soc_bulk_p,'--r','LineWidth',2); hold on
plot(t_out,soc_bulk_n,':g','LineWidth',2)
xlabel('Time [s]','FontSize', 16);ylabel('SOC [-]','FontSize', 16);
set(gca,'FontSize', 16); 
xlim([0, t_data(end)]); ylim([min(SOC_cc)-0.05 max(SOC_cc)+0.05])
legend('Reference','Cathode SOC','Anode SOC'); 
title('SOC profile','FontSize', 16);

%% For Plotting and Voltage RMS error (Plot in the range of 4.2V to 2.5V)

% volt_index = find(V_cell<=2.5,1);
% if isempty(volt_index)
%     volt_index = t_out(end)/dt;
% end
% V_spm = V_cell(1:volt_index);
% time_array = 0:dt:(length(V_spm)-1)*dt; % new time vector
% % select expt. voltage data for new time vector
% V_expt = interp1(t_data,V_data,time_array,'pchip'); 
% % RMS error with arrays of same length
% RMS_V_percentage = rms(V_expt-V_spm)*100/mean(V_expt) 
% 
% figure()
% plot(t_data,V_data,'b','LineWidth',2);hold on; grid on;
% plot(time_array,V_spm,'--r','LineWidth',2)
% xlabel('Time [s]','FontSize', 16);ylabel('Voltage [V]','FontSize', 16);
% set(gca,'FontSize', 16); 
% xlim([0, t_data(end)]); %ylim([min(V_data)-0.01 max(V_data)+0.01])
% legend('Experimental','SPM'); 
% title('Voltage Response','FontSize', 16);
