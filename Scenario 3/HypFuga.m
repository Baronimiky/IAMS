function [DeltaV1,a_H1, e_H1, Delta_H1] = HypFuga(rrf, vvf, v_inf1, mu_T)

% Calcola i parametri dell'iperbole di fuga da un'orbita di parcheggio (caso coplanare)
% INPUT:
%   rrf         - vettore posizione [km] nell'orbita di parcheggio (ECI)
%   vvf         - vettore velocità [km/s] nell'orbita di parcheggio (ECI)
%   v_inf1_vec  - vettore velocità di eccesso iperbolico [km/s] (ECI)
%   mu_T        - parametro gravitazionale del corpo centrale (es. Terra)
%   plotFlag    - true/false per generare il grafico 3D
% OUTPUT:
%   DeltaV1     - ∆v richiesto per la fuga [km/s]
%   v_inf1      - modulo velocità eccesso iperbolico [km/s]
%   v_H1        - velocità al pericentro dell'iperbole [km/s]
%   v_circ      - velocità orbitale nella circolare [km/s]
%   a_H1        - semiasse maggiore iperbole [km]
%   e_H1        - eccentricità iperbole
%   Delta_H1    - parametro di impatto [km]



% Parametri orbita di parcheggio
[a_f, e_f, i_f, om_f, OM_f, th_f] = car2par(rrf, vvf, mu_T);
rp_f = a_f * (1 - e_f);     % pericentro orbitale
r_op = rp_f;
v_op = sqrt(mu_T / r_op); % velocità orbitale nel pericentro

% Parametri iperbolici
v_H1 = sqrt(v_inf1^2 + 2 * mu_T / r_op); % velocità al pericentro iperbole
a_H1 = -mu_T / v_inf1^2;
e_H1 = 1 - r_op/a_H1;
Delta_H1 = -a_H1 * sqrt(e_H1^2 - 1);
DeltaV1 = abs(v_H1 - v_op);

% Output a console
fprintf('v_inf = %.3f km/s\n', v_inf1);
fprintf('r_pericentro = %.2f km\n', r_op);
fprintf('v_circ = %.3f km/s\n', v_op);
fprintf('v_iperbole = %.3f km/s\n', v_H1);
fprintf('DeltaV1 = %.3f km/s\n', DeltaV1);
fprintf('a_H = %.2f km, e_H = %.4f, Delta = %.2f km\n', a_H1, e_H1, Delta_H1);



end

 
