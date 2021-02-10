function [dcs_dt] = ode_spm(t_out,cs,SPMparams,I_data,t_data)
% Solving solid phase PDE using ODE solver

input_crt = interp1(t_data,I_data,t_out);   % Interpolate current

% Separate concentration as per the electrode
cs_n = cs(1:SPMparams.r_grid-1);
cs_p = cs(SPMparams.r_grid:2*(SPMparams.r_grid-1));

%% Solve solid phase ODEs

dcsn_dt = SPMparams.A_mat_n*cs_n + SPMparams.B_mat_n*input_crt;
dcsp_dt = SPMparams.A_mat_p*cs_p + SPMparams.B_mat_p*input_crt;

dcs_dt = [dcsn_dt; dcsp_dt];

end

