function eta_p =get_etap(csp,I_at_k)
run SPMparameters
T_batt=SPMparams.T_batt;
theta_surf_p=(csp)/SPMparams.c_p_max;
% Exchange current density calculation (as a function of surface
% concentration)
i0_p = SPMparams.k_p*SPMparams.F*(SPMparams.ce0^0.5)*...
    ((SPMparams.c_p_max*theta_surf_p)^0.5)*...
    ((SPMparams.c_p_max-SPMparams.c_p_max*theta_surf_p)^0.5) ;   

% Overpotential calculation
eta_p = asinh(I_at_k/(2*SPMparams.A*SPMparams.a_s_p*...
    SPMparams.Ln*i0_p))*(SPMparams.Rg*T_batt)/...
    (SPMparams.F*SPMparams.alpha_cell);
end