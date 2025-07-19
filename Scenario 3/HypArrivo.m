function [DeltaV2,a_H, e_H, Delta_H2] = HypArrivo(rPH,v_inf2,mu_a)

a_H = -mu_a / v_inf2^2;
r_oc = rPH;

% Velocità iperbolica al pericentro
v_pH = sqrt(v_inf2^2 + 2 * mu_a / r_oc);

% Velocità in orbita circolare (ipotetica cattura)
v_oc = sqrt(mu_a / r_oc);

% Delta-V per inserimento
DeltaV2 = abs(v_pH - v_oc);

% Parametri iperbolici
e_H = 1 - r_oc / a_H;
Delta_H2= -a_H * sqrt(e_H^2 - 1);

end