function [c, ceq] = vinc_nlin(x, a_i, e_i, i_i, OM_i, om_i, a_f, e_f, i_f, OM_f, om_f, mu_s)

    th_i = x(1);
    th_f = x(2);
    omega_T_n = x(3);

    % Stato iniziale e finale in coordinate cartesiane
    [r1, ~] = par2car(a_i, e_i, i_i, OM_i, om_i, th_i, mu_s);
    [r2, ~] = par2car(a_f, e_f, i_f, OM_f, om_f, th_f, mu_s);

    % Piano orbitale trasferimento
    h_T = cross(r1, r2) / norm(cross(r1, r2));
    k = [0, 0, 1];
    i_T = acos(dot(h_T, k));

    % Linea dei nodi
    N = cross(k, h_T);
    N = N / norm(N);

    NT_i = dot(N, [1, 0, 0]);
    NT_j = dot(N, [0, 1, 0]);

    if NT_j >= 0
        Omega_T = acos(NT_i);
    else
        Omega_T = 2*pi - acos(NT_i);
    end

    % Matrici di rotazione
    R_omega_T_n = [cos(omega_T_n), sin(omega_T_n), 0;
                  -sin(omega_T_n), cos(omega_T_n), 0;
                   0, 0, 1];

    R_i_T = [1, 0, 0;
             0, cos(i_T), sin(i_T);
             0, -sin(i_T), cos(i_T)];

    R_Omega_T = [cos(Omega_T), sin(Omega_T), 0;
                -sin(Omega_T), cos(Omega_T), 0;
                 0, 0, 1];

    T_Elio_PF_T_n = R_omega_T_n * R_i_T * R_Omega_T;

    r1_T_n = T_Elio_PF_T_n * r1;
    r2_T_n = T_Elio_PF_T_n * r2;

    theta1_T_n = atan2(r1_T_n(2)/norm(r1), r1_T_n(1)/norm(r1));
    theta2_T_n = atan2(r2_T_n(2)/norm(r2), r2_T_n(1)/norm(r2));

    theta1_T_n = PositiveAngles(theta1_T_n);
    theta2_T_n = PositiveAngles(theta2_T_n);

    % Calcolo eccentricità orbitale trasferimento
    e_T_n = (norm(r2) - norm(r1)) / ...
            (norm(r1) * r1_T_n(1)/norm(r1) - norm(r2) * r2_T_n(1)/norm(r2));

    % Vincoli di disuguaglianza: 0 ≤ e_T_n ≤ 1
    c = [0.001-e_T_n;         % e_T_n ≥ 0  → -e_T_n ≤ 0
         e_T_n - 0.999];     % e_T_n ≤ 1  → e_T_n - 1 ≤ 0

    ceq = []; % Nessun vincolo di uguaglianza
end

