function eta_n =get_etan(csn,I_at_k)
run SPMparameters
T_batt=SPMparams.T_batt;
theta_surf_n=(csn)/SPMparams.c_n_max;
% Exchange current density calculation (as a function of surface
% concentration)
i0_n = SPMparams.k_n*SPMparams.F*(SPMparams.ce0^0.5)*...
    ((SPMparams.c_n_max*theta_surf_n)^0.5)*...
    ((SPMparams.c_n_max-SPMparams.c_n_max*theta_surf_n)^0.5);    

% Overpotential calculation
eta_n = asinh(I_at_k/(2*SPMparams.A*SPMparams.a_s_n*...
    SPMparams.Ln*i0_n))*(SPMparams.Rg*T_batt)/...
    (SPMparams.F*SPMparams.alpha_cell);
end