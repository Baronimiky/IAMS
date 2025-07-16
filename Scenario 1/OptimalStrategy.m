function [optimal_time_strategies,optimal_deltav_strategies,minimum_time,minimum_deltav,all_times,all_deltav,strategy_parameters]=OptimalStrategy(ai,ei,Ii,OMi,omi,thi,af,ef,If,OMf,omf,thf,mu)

% function to calculate all the possible 192 (96 + bielliptical cases)
% standard strategies changing the order of the maneuvers and of the
% maneuver points.
% NB: the bielliptical cases are also considered in this function, taking
% the value of rb as the one associated to the minimum dv in a range from 0 to 5*rf
% OUTPUTS:
% St_num: array containing the strategies (in the form of numerical codes)
%         relative to the minimum total time
% Sv_num: array containing the strategies (in the form of numerical codes)
%         relative to the minimum total delta of velocity
% dt_MIN: minimum total time
% dv_MIN: minimum delta of velocity
% dt_vec: array with all the total times
% dv_vec: array with all the total delta v
% param_vec: array containing all the useful parameters for all the strategies
%            [rb_val,dv_cp,dv_ca,dv1,dv2,dv3,dt_cp,dt_ca,dt_trans,dt_cs_t,dt_fin,om_cp,om_ca,th_cp,thf_old,thf_new,th_cs]

% dt_trans: time to reach the maneuver point for the change of shape starting from the last point of the previous maneuver

%LEGEND:
% how to interpret the numerical codes assosciated to the strategies in St_num and Sv_num:
% - first number (c1):
%         1: plane-axes-shape
%         2: plane-shape-axes
%         3: axes-plane-shape
%         4: shape-plane-axes
%         5: axes-shape-plane
%         6: shape-axes-plane
% - second number (c2):
%         1: bitangent Pericenter-Apocenter
%         2: bitangent Pericenter-Pericenter
%         3: bitangent Apocenter-Pericenter
%         4: bitangent Apocenter-Apocenter
%         5: bielliptical Pericenter-Pericenter
%         6: bielliptical Pericenter-Apocenter
%         7: bielliptical Apocenter-Apocenter
%         8: bielliptical Apocenter-Pericenter
% - third number (c3):
%         maneuver point for change of plane (1 or 2)
% - fourth number (c4):
%         maneuver point for change of periapsis (1 or 2)


ra_f = af * (1 + ef);  % Apocentro finale

% Preallocazione matrici risultati
c_matr = zeros(192,4);
all_times = zeros(192,1);
all_deltav = zeros(192,1);
strategy_parameters = zeros(192,17);

% Precalcolo configurazioni per c2
phi_values = [0, pi, 0, pi, 0, pi, 0, pi];
P_values = [1, 1, 2, 2, 1, 1, 2, 2];
th_cs_1_values = [2*pi, 2*pi, pi, pi, 2*pi, 2*pi, pi, pi];

cont = 1;

for c1 = 1:6
    for c2 = 1:8
        phi = phi_values(c2);
        P = P_values(c2);
        th_cs_1 = th_cs_1_values(c2);
        is_bielliptic = c2 > 4;
        
        for c3 = 1:2
            for c4 = 1:2
                % Inizializzazione variabili
                rb_val = 0; dv1 = 0; dv2 = 0; dv3 = 0; th_cs = 0;
                dv_cp = 0; dv_ca = 0; dt_cp = 0; dt_ca = 0; 
                dt_trans = 0; dt_cs_t = 0; dt_fin = 0;
                om_cp = 0; om_ca = 0; thi_ca = 0; thf_ca = 0; th_cp = 0;
                
                switch c1
                    case 1  % P-A-S
                        [dv_cp, om_cp, th_cp, dt_cp] = changeOrbitalPlane(ai,ei,Ii,OMi,omi,If,OMf,thi,mu);
                        [dv_ca, om_ca, thi_ca, thf_ca, dt_ca] = changePericenterArg(ai,ei,om_cp,omf+phi,th_cp,P,mu);
                        dt_trans = TOF(ai,ei,thf_ca,th_cs_1,mu);
                        
                        if is_bielliptic
                            [min_dv, min_dt, best_rb, best_dv1, best_dv2, best_dv3, best_th_cs] = optimizeBielliptic(ai,ei,om_ca,af,ef,omf,P,ra_f,mu);
                            dv_cs = min_dv;
                            dt_cs_t = min_dt;
                            rb_val = best_rb;
                            dv1 = best_dv1; dv2 = best_dv2; dv3 = best_dv3; th_cs = best_th_cs;
                        else
                            [dv1,dv2,dv3,th_cs,dt_cs_temp] = changeOrbitShape(ai,ei,om_ca,af,ef,omf,P,mu);
                            dv_cs = abs(dv1) + abs(dv2) + abs(dv3);
                            dt_cs_t = dt_cs_temp;
                        end
                        dt_fin = TOF(af,ef,th_cs,thf,mu);
                        dt_cs = dt_cs_t + dt_trans;
                        
                    case 2  % P-S-A
                        [dv_cp, om_cp, th_cp, dt_cp] = changeOrbitalPlane(ai,ei,Ii,OMi,omi,If,OMf,thi,mu);
                        dt_trans = TOF(ai,ei,th_cp,th_cs_1,mu);
                        
                        if is_bielliptic
                            [min_dv, min_dt, best_rb, best_dv1, best_dv2, best_dv3, best_th_cs] = optimizeBielliptic(ai,ei,om_cp,af,ef,om_cp+phi,P,ra_f,mu);
                            dv_cs = min_dv;
                            dt_cs_t = min_dt;
                            rb_val = best_rb;
                            dv1 = best_dv1; dv2 = best_dv2; dv3 = best_dv3; th_cs = best_th_cs;
                        else
                            [dv1,dv2,dv3,th_cs,dt_cs_temp] = changeOrbitShape(ai,ei,om_cp,af,ef,om_cp+phi,P,mu);
                            dv_cs = abs(dv1) + abs(dv2) + abs(dv3);
                            dt_cs_t = dt_cs_temp;
                        end
                        dt_cs = dt_cs_t + dt_trans;
                        [dv_ca, om_ca, thi_ca, thf_ca, dt_ca] = changePericenterArg(af,ef,om_cp+phi,omf,th_cs,P,mu);
                        dt_fin = TOF(af,ef,thf_ca,thf,mu);
                        
                    case 3  % A-P-S
                        om = InverseCP(Ii,OMi,omf-phi,If,OMf);
                        [dv_ca, om_ca, thi_ca, thf_ca, dt_ca] = changePericenterArg(ai,ei,omi,om,thi,c4,mu);
                        [dv_cp, om_cp, th_cp, dt_cp] = changeOrbitalPlane(ai,ei,Ii,OMi,om_ca,If,OMf,thf_ca,mu);
                        dt_trans = TOF(ai,ei,th_cp,th_cs_1,mu);
                        
                        if is_bielliptic
                            [min_dv, min_dt, best_rb, best_dv1, best_dv2, best_dv3, best_th_cs] = optimizeBielliptic(ai,ei,om_cp,af,ef,omf,P,ra_f,mu);
                            dv_cs = min_dv;
                            dt_cs_t = min_dt;
                            rb_val = best_rb;
                            dv1 = best_dv1; dv2 = best_dv2; dv3 = best_dv3; th_cs = best_th_cs;
                        else
                            [dv1,dv2,dv3,th_cs,dt_cs_temp] = changeOrbitShape(ai,ei,om_cp,af,ef,omf,P,mu);
                            dv_cs = abs(dv1) + abs(dv2) + abs(dv3);
                            dt_cs_t = dt_cs_temp;
                        end
                        dt_cs = dt_cs_t + dt_trans;
                        dt_fin = TOF(af,ef,th_cs,thf,mu);
                        
                    case 4  % S-P-A
                        dt_trans = TOF(ai,ei,thi,th_cs_1,mu);
                        
                        if is_bielliptic
                            [min_dv, min_dt, best_rb, best_dv1, best_dv2, best_dv3, best_th_cs] = optimizeBielliptic(ai,ei,omi,af,ef,omi+phi,P,ra_f,mu);
                            dv_cs = min_dv;
                            dt_cs_t = min_dt;
                            rb_val = best_rb;
                            dv1 = best_dv1; dv2 = best_dv2; dv3 = best_dv3; th_cs = best_th_cs;
                        else
                            [dv1,dv2,dv3,th_cs,dt_cs_temp] = changeOrbitShape(ai,ei,omi,af,ef,omi+phi,P,mu);
                            dv_cs = abs(dv1) + abs(dv2) + abs(dv3);
                            dt_cs_t = dt_cs_temp;
                        end
                        dt_cs = dt_cs_t + dt_trans;
                        [dv_cp, om_cp, th_cp, dt_cp] = changeOrbitalPlane(af,ef,Ii,OMi,omi+phi,If,OMf,th_cs,mu);
                        [dv_ca, om_ca, thi_ca, thf_ca, dt_ca] = changePericenterArg(af,ef,om_cp,omf,th_cp,c4,mu);
                        dt_fin = TOF(af,ef,thf_ca,thf,mu);
                        
                    case 5  % A-S-P
                        om = InverseCP(Ii,OMi,omf,If,OMf);
                        [dv_ca, om_ca, thi_ca, thf_ca, dt_ca] = changePericenterArg(ai,ei,omi,om+phi,thi,c4,mu);
                        dt_trans = TOF(ai,ei,thf_ca,th_cs_1,mu);
                        
                        if is_bielliptic
                            [min_dv, min_dt, best_rb, best_dv1, best_dv2, best_dv3, best_th_cs] = optimizeBielliptic(ai,ei,om_ca,af,ef,om,P,ra_f,mu);
                            dv_cs = min_dv;
                            dt_cs_t = min_dt;
                            rb_val = best_rb;
                            dv1 = best_dv1; dv2 = best_dv2; dv3 = best_dv3; th_cs = best_th_cs;
                        else
                            [dv1,dv2,dv3,th_cs,dt_cs_temp] = changeOrbitShape(ai,ei,om_ca,af,ef,om,P,mu);
                            dv_cs = abs(dv1) + abs(dv2) + abs(dv3);
                            dt_cs_t = dt_cs_temp;
                        end
                        dt_cs = dt_cs_t + dt_trans;
                        [dv_cp, om_cp, th_cp, dt_cp] = changeOrbitalPlane(af,ef,Ii,OMi,om,If,OMf,th_cs,mu);
                        dt_fin = TOF(af,ef,th_cp,thf,mu);
                        
                    case 6  % S-A-P
                        dt_trans = TOF(ai,ei,thi,th_cs_1,mu);
                        
                        if is_bielliptic
                            [min_dv, min_dt, best_rb, best_dv1, best_dv2, best_dv3, best_th_cs] = optimizeBielliptic(ai,ei,omi,af,ef,omi+phi,P,ra_f,mu);
                            dv_cs = min_dv;
                            dt_cs_t = min_dt;
                            rb_val = best_rb;
                            dv1 = best_dv1; dv2 = best_dv2; dv3 = best_dv3; th_cs = best_th_cs;
                        else
                            [dv1,dv2,dv3,th_cs,dt_cs_temp] = changeOrbitShape(ai,ei,omi,af,ef,omi+phi,P,mu);
                            dv_cs = abs(dv1) + abs(dv2) + abs(dv3);
                            dt_cs_t = dt_cs_temp;
                        end
                        dt_cs = dt_cs_t + dt_trans;
                        om = InverseCP(Ii,OMi,omf,If,OMf);
                        [dv_ca, om_ca, thi_ca, thf_ca, dt_ca] = changePericenterArg(af,ef,omi+phi,om,th_cs,c4,mu);
                        [dv_cp, om_cp, th_cp, dt_cp] = changeOrbitalPlane(af,ef,Ii,OMi,om,If,OMf,thf_ca,mu);
                        dt_fin = TOF(af,ef,th_cp,thf,mu);
                end
                
                % Calcolo costi totali
                total_dv = abs(dv_cs) + abs(dv_cp) + abs(dv_ca);
                total_time = dt_cp + dt_cs + dt_ca + dt_fin;
                
                % Memorizzazione risultati
                strategy_parameters(cont,:) = [rb_val,dv_cp,dv_ca,dv1,dv2,dv3,dt_cp,dt_ca,dt_trans,dt_cs_t,dt_fin,om_cp,om_ca,th_cp,thi_ca,thf_ca,th_cs];
                all_deltav(cont) = total_dv;
                all_times(cont) = total_time;
                c_matr(cont,:) = [c1,c2,c3,c4];
                cont = cont + 1;
            end
        end
    end
end

% Ricerca strategie ottimali
minimum_deltav = min(all_deltav);
minimum_time = min(all_times);

deltav_optimal_indices = find(all_deltav == minimum_deltav);
time_optimal_indices = find(all_times == minimum_time);

optimal_deltav_strategies = [c_matr(deltav_optimal_indices,:), all_times(deltav_optimal_indices)];
optimal_time_strategies = [c_matr(time_optimal_indices,:), all_deltav(time_optimal_indices)];

end

% Funzione ottimizzata per trasferimenti biellittici con ricerca a due stadi
function [min_dv, min_dt, best_rb, best_dv1, best_dv2, best_dv3, best_th_cs] = optimizeBielliptic(ai,ei,om_input,af,ef,om_output,P,ra_f,mu)
    % Fase 1: Ricerca grossolana
    coarse_step = 5*ra_f / 100; % Passo iniziale ricerca grossolana
    rb_vect_coarse = 0:coarse_step:5*ra_f; % Vettore ricerca iniziale
    
    [min_dv, min_dt, best_rb, best_dv1, best_dv2, best_dv3, best_th_cs] = innerBiellipticSearch(ai,ei,om_input,af,ef,om_output,P,rb_vect_coarse,mu);
    
    % Fase 2: Ricerca raffinata (solo se abbiamo trovato un minimo valido)
    if best_rb > 0 && min_dv < Inf
        % Definisci l'intervallo di raffinamento (±10% del passo grossolano)
        refine_range = coarse_step * 0.5;
        low_bound = max(0, best_rb - refine_range);
        high_bound = min(5*ra_f, best_rb + refine_range);
        
        % Crea griglia più fine (1000 punti nell'intervallo ridotto)
        rb_vect_fine = linspace(low_bound, high_bound, 1000);
        
        % Esegui ricerca raffinata
        [min_dv_fine, min_dt_fine, best_rb_fine, best_dv1_fine, best_dv2_fine, best_dv3_fine, best_th_cs_fine] = ...
            innerBiellipticSearch(ai,ei,om_input,af,ef,om_output,P,rb_vect_fine,mu);
        
        % Se abbiamo trovato un minimo migliore, aggiorna i risultati
        if min_dv_fine < min_dv
            min_dv = min_dv_fine;
            min_dt = min_dt_fine;
            best_rb = best_rb_fine;
            best_dv1 = best_dv1_fine;
            best_dv2 = best_dv2_fine;
            best_dv3 = best_dv3_fine;
            best_th_cs = best_th_cs_fine;
        end
    end
end

% Funzione helper per ricerca biellittica
function [min_dv, min_dt, best_rb, best_dv1, best_dv2, best_dv3, best_th_cs] = innerBiellipticSearch(ai,ei,om_input,af,ef,om_output,P,rb_vect,mu)
    min_dv = inf;
    min_dt = inf;
    best_rb = 0;
    best_dv1 = 0; best_dv2 = 0; best_dv3 = 0; best_th_cs = 0;
    
    for i = 1:length(rb_vect)
        rb = rb_vect(i);
        [dv1,dv2,dv3,th_cs,dt_cs_temp] = changeOrbitShape(ai,ei,om_input,af,ef,om_output,P,mu,rb);
        total_dv = abs(dv1) + abs(dv2) + abs(dv3);
        
        if total_dv < min_dv
            min_dv = total_dv;
            min_dt = dt_cs_temp;
            best_rb = rb;
            best_dv1 = dv1;
            best_dv2 = dv2;
            best_dv3 = dv3;
            best_th_cs = th_cs;
        end
    end
end