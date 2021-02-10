function y = etap_jacobian(csp,I_at_k)
%deltacsp for linearization= 0.01
delta=0.01;

eta_plus=get_etap(csp+delta,I_at_k);
eta_minus=get_etap(csp-delta,I_at_k);
y=(eta_plus-eta_minus)/(2*delta);
end