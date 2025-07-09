function DeltaV_benchmark = benchmarkDeltaV(a_i, e_i, i_i, OM_i, om_i, ...
                                          a_f, e_f, i_f, OM_f, om_f, mu)

    % Calcolo raggio orbitale medio (circolare equivalente) per iniziale e finale
    r_i = a_i * (1 - e_i^2) / (1 + e_i);  % perigeo come approssimazione
    r_f = a_f * (1 - e_f^2) / (1 + e_f);  % perigeo orbita finale

    % Velocità circolare ai due raggi
    V_i = sqrt(mu / r_i);
    V_f = sqrt(mu / r_f);

    % Delta V trasferimento di Hohmann (da r_i a r_f)
    dV1 = abs(V_i * (sqrt(2 * r_f / (r_i + r_f)) - 1));
    dV2 = abs(V_f * (1 - sqrt(2 * r_i / (r_i + r_f))));
    DeltaV_Hohmann = dV1 + dV2;

    % Cambio inclinazione orbitale (usiamo la media tra le due velocità)
    V_avg = (V_i + V_f) / 2;
    delta_i = abs(i_f - i_i);
    DeltaV_plane = 2 * V_avg * sin(delta_i / 2);

    % Somma totale
    DeltaV_benchmark = DeltaV_Hohmann + DeltaV_plane;

    fprintf('Benchmark Hohmann + cambio piano = %.4f km/s\n', DeltaV_benchmark);
end
