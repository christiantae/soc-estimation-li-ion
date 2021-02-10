%% SPM EKF Script
%  Extended Kalman Filter (EKF) implementation for State of Charge (SOC)
%  estimation using a second-order equivalent circuit battery model.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
%close all; 
clc; 
warning off

%% Load experimental data

%load HW5_expt_Data.mat
%load NMC_Cell_H1_T23_UDDS.mat; dt = 0.1;
load('US06_test1.mat'); dt = 0.1; 
%load('HW5_batt_params.mat')
capacity = 2;
%load UDDS23_H1_training.mat

% Define current, voltage, and time data
% Current sign convention: Positive for discharge
%I_data = I_data;
%V_data = V_data_spm;
 %I_data=I_expt;
 %V_data=V_expt;
t_data = 0:dt:(length(I_data)-1)*dt;   

%% Initialize model parameters

%load HW5_batt_params.mat
run SPMparameters.m

%% Calculate Reference SOC using Coulomb Counting

SOC_ref = zeros(length(I_data),1);
SOC_ref(1) = 1;
for k=2:length(I_data)
    SOC_ref(k) = SOC_ref(k-1)-I_data(k-1)*dt/(capacity*3600);
end

%% EKF Initialization

R = 8.432e-4; % measurement noise covariance
Q = diag([1000*R .1*R .01*R 0.01*R]); % process noise covariance
P_pred(:,:,1) = [.1 0 0 0;0 0.001 0 0;0 0 0.001 0;0 0 0 0.001]; % initial covariance

SPMparams.r_grid = 5;                  % Radial discretization grid points
soc_init= 0.55;

% Initialize concentration states in both electrodes
for i=1:SPMparams.r_grid
    cs_p0(i,1)= ((SPMparams.theta100_p-SPMparams.theta0_p)*(soc_init)+SPMparams.theta0_p)*SPMparams.c_p_max;
end

cs_bulk_p(1)=cs_p0(1,1);
soc_bulk_p(1)=soc_init;

x_pred = cs_p0(2:end); % initial predicted state
x_corr = [0;0;0;0]; % initial corrected state

delta_soc = 0.01; % delta soc for linearization

%% Solid phase: Finite Difference Method Discretization 

% Grid interval size for anode and cathode    
SPMparams.delta_x_p = SPMparams.Rs_p/(SPMparams.r_grid-1);  
interval_p = [SPMparams.delta_x_p:SPMparams.delta_x_p:SPMparams.delta_x_p*(SPMparams.r_grid-1)]';

% Coefficientes of discretized ODEs
SPMparams.alpha_p = SPMparams.Ds_p/(SPMparams.delta_x_p^2);
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
A = A_mat.*SPMparams.alpha_p;
B = B_mat.*SPMparams.beta_p; 
C= []; 
D=[SPMparams.R_c]; 

%change from matrices from continuous time to discrete time
SPM_CT_ss=ss(A,B,C,D);
SPM_DT_ss=c2d(SPM_CT_ss,dt);
A=SPM_DT_ss.a;
B=SPM_DT_ss.b;

%% EKF algorithm for the drive cycle (input profile)
for i = 2:length(I_data)
%% Predictor step (A-priori)
  % Predict state
    x_pred(:,i) = A*x_pred(:,i-1) + B*I_data(i-1); % 

  % Linearization around the predicted state
    etap_jb(i)=etap_jacobian(x_pred(end,i),I_data(i));
    up_jb(i)=Up_jacobain(x_pred(end,i));

  % C matrix of the State Space Model (continuous time)
    H = [0, 0, 0, (etap_jb(i)+up_jb(i))];
    %   Ob = obsv(A,C); % Check observability matrix
    %   unob(i) = length(A)-rank(Ob); % number of unobservable states

  % Predict covariance    
    P_pred(:,:,i) = A*P_pred(:,:,i-1)*A' + Q; 

%% Correction step (A-posteriori)
    
  % kalman gain    
    K(:,i) = P_pred(:,:,i)*H'/(H*P_pred(:,:,i)*H' + R);  
    
  % Innovation calculation
  csp(i)= x_pred(end,i);
  csn(i)= (((csp(i)/SPMparams.c_p_max-SPMparams.theta100_p)/(SPMparams.theta0_p-SPMparams.theta100_p))*(SPMparams.theta0_n-SPMparams.theta100_n)+SPMparams.theta100_n)*SPMparams.c_n_max;
  ocp_p(i) = U_p((csp(i)/SPMparams.c_p_max));
  ocp_n(i) = U_n((csn(i)/SPMparams.c_n_max));
  eta_p(i)=get_etap(csp(i),I_data(i));
  eta_n(i)=get_etan(csn(i),I_data(i));
  
  V_pred(i)= ocp_p(i)-ocp_n(i)+eta_p(i)-eta_n(i)-I_data(i).*SPMparams.R_c;

  Inn(i) = V_data(i)-V_pred(i);
    
  % update corrected state
    x_corr(:,i) = x_pred(:,i) + K(:,i)*Inn(i); 
    
  % update covariance    
    P_corr(:,:,i) = (eye(length(A))- K(:,i)*H)*P_pred(:,:,i); 
    
%% Output prediction
    
  % Cell voltage
  %re-calculate with updated x vector
  csp(i)=x_corr(end,i);
  csn(i)= (((csp(i)/SPMparams.c_p_max-SPMparams.theta100_p)/(SPMparams.theta0_p-SPMparams.theta100_p))*(SPMparams.theta0_n-SPMparams.theta100_n)+SPMparams.theta100_n)*SPMparams.c_n_max;
  ocp_p(i) = U_p((csp(i)/SPMparams.c_p_max));
  ocp_n(i) = U_n((csn(i)/SPMparams.c_n_max));
  eta_p(i)=get_etap(csp(i),I_data(i));
  eta_n(i)=get_etan(csn(i),I_data(i));
  
  V_est(i,:) = ocp_p(i)-ocp_n(i)+eta_p(i)-eta_n(i)-I_data(i).*SPMparams.R_c;
  
  cs_bulk_p(i) = trapz(interval_p,interval_p.^2.*x_corr(1:end,i))*3/(SPMparams.Rs_p^3);
  soc_bulk_p(i) = (SPMparams.theta0_p - cs_bulk_p(i)/SPMparams.c_p_max)/(SPMparams.theta0_p - SPMparams.theta100_p); 

  %soc_bulk_p(i) = soc_bound(soc_bulk_p(i));

    
%% Update state and covariance matrix for next time step
    P_pred(:,:,i) = P_corr(:,:,i);
    x_pred(:,i) = x_corr(:,i);

end

%just renamed it here for plotting
 SOC_est = soc_bulk_p;

%% RMS Error calculation and plots
% RMS_V_percentage = rms(V_data-V_est)*100/mean(V_data);
% RMS_soc_percentage = rms(SOC_ref-SOC_est')*100/mean(SOC_ref);
% 
% figure()
% plot(t_data,V_data,'b','LineWidth',2);hold on; grid on;
% plot(t_data,V_est,'--r','LineWidth',2)
% xlabel('Time [s]','FontSize', 16);ylabel('Voltage [V]','FontSize', 16);
% set(gca,'FontSize', 16); 
% xlim([0, t_data(end)]); ylim([min(V_data)-0.01 max(V_data)+0.01])
% legend('Experimental','EKF Estimated'); 
% title('Voltage Response','FontSize', 16);
% 
% figure()
% plot(t_data,SOC_ref,'b','LineWidth',2);hold on; grid on;
% plot(t_data,SOC_est,'--r','LineWidth',2)
% xlabel('Time [s]','FontSize', 16);ylabel('SOC [-]','FontSize', 16);
% set(gca,'FontSize', 16); 
% xlim([0, t_data(end)]); ylim([min(SOC_ref)-0.001 max(SOC_ref)+0.001])
% legend('Reference','EKF Estimated'); 
% title('SOC Profile','FontSize', 16);
% 
% figure()
% plot(t_data,SOC_ref-SOC_est,'b','LineWidth',2);hold on; grid on;
% xlabel('Time [s]','FontSize', 16);ylabel('SOC [-]','FontSize', 16);
% set(gca,'FontSize', 16); 
% xlim([0, t_data(end)]); 
% title('SOC Estimation Error (EKF)','FontSize', 16);
% 
% figure()
% plot(t_data,K(1,:),'b','LineWidth',2);hold on; grid on;
% xlabel('Time [s]','FontSize', 16);ylabel('Kalman Gain [-]','FontSize', 16);
% set(gca,'FontSize', 16); 
% xlim([0, t_data(end)]); 
% title('Kalman Gain for SOC state','FontSize', 16);