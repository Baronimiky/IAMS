clear all
close all
clc

%% CARICAMENTO DATI
data= load("DatiSC1-2025.txt");
Or_i = data(27,2:7);
Or_f = data(27,8:13);

% orbita iniziale
a_i = Or_i(1,1); % semiasse maggiore[km] 
e_i = Or_i(1,2); % eccentricità [-]
i_i = Or_i(1,3); % inclinazione del piano orbitale[rad] inc_i
OM_i = Or_i(1,4); % ascensione retta del nodo ascendente[rad] RAAN
om_i = Or_i(1,5); % anomalia del pericentro [rad] w
th_i = Or_i(1,6); % posizione del satellite nell'orbita [rad] theta_i

% orbita finale
rx_f = Or_f(1,1); % [km]
ry_f = Or_f(1,2); % [km]
rz_f = Or_f(1,3); % [km]
vx_f = Or_f(1,4); % [km/s]
vy_f = Or_f(1,5); % [km/s]
vz_f = Or_f(1,6); % [km/s]

rr_f = [rx_f; ry_f; rz_f]; % raggio punto di arrivo orbita
vv_f = [vx_f; vy_f;vz_f]; % velocità punto di arrivo orbita

mu_t = 398600; % costante gravitazionale terra [km^3/s^2]
mu_s = 132712440018; % costante gravitazionale sole [km^3/s^2]

% Conversione deli parametri orbitali in coordinate cartesiane
[rr_i,vv_i] = par2car(a_i,e_i,i_i,OM_i,om_i,th_i,mu_t);

% parametri utili per plot orbite
dth = 0.01; % discriminante theta
Rt = 6378; % raggio terra
% Conversione delle coordinate cartesiane in parametri orbitali
[a_f,e_f,i_f,OM_f,om_f,th_f] = car2par(rr_f,vv_f,mu_t);

%% Strategia ottimale e plot
[St_num,Sv_num,dt_MIN,dv_MIN,dt_vec,dv_vec,parameters_all]=OptimalStrategy(a_i,e_i,i_i,OM_i,om_i,th_i,a_f,e_f,i_f,OM_f,om_f,th_f,mu_t);

%% FUNZIONE ANIMAZIONE
function animateTransfer(a,e,i,OM,om,th_start,th_end,mu,satellite,v,doRotate, view_start, view_end, time_offset, h_text_time, h_text_vel)
    dth_anim = 0.02;
    
    % Calcola l'anomalia eccentrica per th_start e th_end
    E_start = 2*atan(sqrt((1-e)/(1+e)) * tan(th_start/2));
    E_end = 2*atan(sqrt((1-e)/(1+e)) * tan(th_end/2));
    
    % Garantisce E_end > E_start
    if E_end < E_start
        E_end = E_end + 2*pi;
    end
    
    % Calcola il tempo esatto usando l'equazione di Kepler
    n = sqrt(mu/a^3); % moto medio
    DeltaT_exatto = (E_end - e*sin(E_end) - (E_start - e*sin(E_start))) / n;
    
    % Genera i punti lungo l'orbita
    th_vals = th_start:dth_anim:th_end;
    if th_end < th_start
        th_vals = [th_start:dth_anim:2*pi, 0:dth_anim:th_end];
    end
    N = length(th_vals);
    
    for k = 1:N
        th = th_vals(k);
        [rr, vv] = par2car(a,e,i,OM,om,mod(th, 2*pi),mu);
        set(satellite, 'XData', rr(1), 'YData', rr(2), 'ZData', rr(3));
        
        % Calcola velocità reale in quel punto
        v_norm = norm(vv);
        
        % Calcola tempo trascorso proporzionale al tempo esatto
        time_passed = time_offset + DeltaT_exatto * (k-1)/N;
        
        % Aggiorna i textbox
        set(h_text_time, 'String', sprintf('Tempo: %.1f s', time_passed));
        set(h_text_vel,  'String', sprintf('Velocità: %.2f km/s', v_norm));
        
        % Solo se è abilitata la rotazione della vista
        if doRotate
            az = view_start(1) + (view_end(1) - view_start(1)) * k / N;
            el = view_start(2) + (view_end(2) - view_start(2)) * k / N;
            view(az, el);
        end
        
        drawnow;
        frame = getframe(gcf);
        writeVideo(v, frame);
    end
    
    % Alla fine del blocco, fissa la vista
    if doRotate
        view(view_end);
    end
end

%% ANIMAZIONE STRATEGIA DELTA-V MINIMA

% Trova strategia ottimale per velocità minima
deltav_optimal_indices = find(dv_vec == dv_MIN & dt_vec == min(Sv_num(:,5)));
optimal_index_v = deltav_optimal_indices(1);

% Setup figura
fig1 = figure('Position', [100, 100, 1000, 600]);
hold on; grid on; axis equal;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title(['Strategia Ottimale Delta-V Minimo #' num2str(optimal_index_v)]);

% Plot orbite
plotOrbit(a_i,e_i,i_i,OM_i,om_i,0,2*pi,dth,mu_t, 'b'); % orbita iniziale
plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,dth,mu_t, 'r'); % orbita finale

% Calcolo parametri per strategia V_MIN
type = 'aa';
[DeltaV1_BT, DeltaV2_BT, DeltaT4, om_bt] = bitangentTransfer(a_i, e_i, a_f, e_f, om_i, type, mu_t);

% Plot punti di manovra
[rr_0] = par2car(a_i,e_i,i_i,OM_i,om_i,th_i,mu_t);
plot3(rr_0(1),rr_0(2),rr_0(3), '.', 'MarkerSize', 40, 'LineWidth', 4);

thi_bt = pi; thf_bt = pi; % per trasferimento aa
[rr_3] = par2car(a_i,e_i,i_i,OM_i,om_i,thi_bt,mu_t);
plot3(rr_3(1),rr_3(2),rr_3(3),'k+','LineWidth',2,'MarkerSize', 10);

% Parametri orbita di trasferimento
r_apo_i = a_i * (1 + e_i);
r_apo_f = a_f * (1 + e_f);
a_bt = (r_apo_i + r_apo_f) / 2;
e_bt = abs(r_apo_f - r_apo_i) / (r_apo_f + r_apo_i);
i_bt = i_i; OM_bt = OM_i;

% Plot orbite intermedie
plotOrbit(a_bt, e_bt, i_bt, OM_bt, om_i, pi, 2*pi, dth, mu_t, 'm');
plotOrbit(a_f, e_f, i_bt, OM_bt, om_bt, 0, 2*pi, dth, mu_t, 'black');

% Calcolo altre manovre
[DeltaV_CP, om_2, theta_CP, DeltaT1] = changeOrbitalPlane(a_f, e_f, i_bt, OM_bt, om_bt, i_f, OM_f, thf_bt, mu_t);
[rr_1] = par2car(a_f,e_f,i_bt,OM_bt,om_bt,theta_CP,mu_t);
plot3(rr_1(1),rr_1(2),rr_1(3),'k+','LineWidth',2,'MarkerSize', 10);

plotOrbit(a_f,e_f,i_f,OM_f,om_2,0,2*pi,dth,mu_t, 'g');

[DeltaV_CA, om_3, thi_CPA, thf_CPA, DeltaT2] = changePericenterArg(a_f, e_f, om_2, om_f, theta_CP, 1, mu_t);
[rr_2] = par2car(a_f,e_f,i_f,OM_f,om_2,thi_CPA,mu_t);
plot3(rr_2(1),rr_2(2),rr_2(3),'k+','LineWidth',2,'MarkerSize', 20);

plot3(rr_f(1),rr_f(2),rr_f(3), '.', 'MarkerSize', 40, 'LineWidth', 4,'Color','m');
Terra3d(Rt);

% Tempi di volo
DeltaT3 = TOF(a_i,e_i,th_i,thi_bt,mu_t);
DeltaT5 = TOF(a_f, e_f, thf_bt, th_f, mu_t);

% Legenda
legend('Orbita Iniziale', 'Orbita Finale', 'Posizione Iniziale Satellite', ...
      'Punto Manovra Trasferimento Bitangente', 'Orbita Trasferimento Bitangente', ...
      'Orbita dopo Trasferimento Bitangente', 'Orbita dopo Cambio Piano', ...
      'Punto Manovra Cambio Piano', 'Punto Manovra Cambio Argomento Pericentro', 'Posizione Finale Satellite','Location','northeastoutside');


% Setup animazione
view_start = [0, 20]; view_end = [120, 34];
pos_time = [0.82 0.13 0.15 0.05]; pos_vel = [0.82 0.13 0.15 0.05];
h_text_time = annotation('textbox', pos_time, 'String', '', 'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
h_text_vel = annotation('textbox', pos_vel, 'String', '', 'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');

satellite = plot3(NaN, NaN, NaN, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 8, 'HandleVisibility', 'off');

% Video writer per V_MIN
v1 = VideoWriter('Satellite_Transfer_VMIN.mp4', 'MPEG-4');
v1.Quality = 100; v1.FrameRate = 120; open(v1);

% Esecuzione animazione V_MIN
time_accum = 0;
animateTransfer(a_i, e_i, i_i, OM_i, om_i, th_i, thi_bt, mu_t, satellite, v1, true, view_start, view_end, time_accum, h_text_time, h_text_vel);
time_accum = time_accum + DeltaT3;

animateTransfer(a_bt, e_bt, i_bt, OM_bt, om_i, thi_bt, thf_bt+pi, mu_t, satellite, v1, false, view_start, view_end, time_accum, h_text_time, h_text_vel);
time_accum = time_accum + DeltaT4;

animateTransfer(a_f, e_f, i_bt, OM_bt, om_bt, thf_bt, PositiveAngles(theta_CP), mu_t, satellite, v1, false, view_start, view_end, time_accum, h_text_time, h_text_vel);
time_accum = time_accum + DeltaT1;

animateTransfer(a_f, e_f, i_f, OM_f, om_2, PositiveAngles(theta_CP), thi_CPA, mu_t, satellite, v1, false, view_start, view_end, time_accum, h_text_time, h_text_vel);
time_accum = time_accum + DeltaT2;

animateTransfer(a_f, e_f, i_f, OM_f, om_f, thf_CPA, th_f+2*pi, mu_t, satellite, v1, false, view_start, view_end, time_accum, h_text_time, h_text_vel);

close(v1);
disp('Animazione V_MIN completata. Video salvato come "Satellite_Transfer_VMIN.mp4"');

%% ANIMAZIONE STRATEGIA TEMPO MINIMO
pos_time = [0.82 0.13 0.15 0.05]; pos_vel = [0.82 0.07 0.15 0.05];
view_end = [125, 37];

% Trova strategia ottimale per tempo minimo
deltat_optimal_indices = find(dt_vec == dt_MIN & dv_vec == min(St_num(:,5)));
optimal_index_t = deltat_optimal_indices(1);

% Setup nuova figura
fig2 = figure('Position', [100, 100, 1000, 600]);
hold on; grid on; axis equal;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title(['Animazione Strategia Tempo Minimo #' num2str(optimal_index_t)]);

% Plot orbite
plotOrbit(a_i,e_i,i_i,OM_i,om_i,0,2*pi,dth,mu_t,'b');
plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,dth,mu_t,'r');

% Plot posizione iniziale
[rr_0] = par2car(a_i,e_i,i_i,OM_i,om_i,th_i,mu_t);
plot3(rr_0(1),rr_0(2),rr_0(3), '.', 'MarkerSize', 40, 'LineWidth', 4);

% 1. CAMBIO ASSI
[omi] = Inverse_changeplane(i_i,OM_i,om_f,i_f,OM_f);
[DeltaV_CA, om_ca_new, thi_CPA, thf_CPA, DeltaT_CA] = changePericenterArg(a_i, e_i, om_i, omi, th_i, 1, mu_t);

[rr_ca] = par2car(a_i,e_i,i_i,OM_i,om_i,thi_CPA,mu_t);
plot3(rr_ca(1),rr_ca(2),rr_ca(3),'k+','LineWidth',2,'MarkerSize', 10);
plotOrbit(a_i,e_i,i_i,OM_i,omi,0,2*pi,dth,mu_t, 'c');

% 2. CAMBIO PIANO
[DeltaV_CP, om_cp_new, theta_CP, DeltaT_CP] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, omi, i_f, OM_f, thf_CPA, mu_t);
[rr_cp] = par2car(a_i,e_i,i_i,OM_i,omi,theta_CP,mu_t);
plot3(rr_cp(1),rr_cp(2),rr_cp(3),'k+','LineWidth',2,'MarkerSize', 10);
plotOrbit(a_i,e_i,i_f,OM_f,om_cp_new,0,2*pi,dth,mu_t, 'g');

% 3. TRASFERIMENTO BITANGENTE AP
type = 'ap';
[DeltaV1_BT, DeltaV2_BT, DeltaT_BT, om_bt] = bitangentTransfer(a_i, e_i, a_f, e_f, om_cp_new, type, mu_t);

thi_bt = pi; thf_bt = 0; % ap transfer
[rr_bt_start] = par2car(a_i,e_i,i_f,OM_f,om_cp_new,thi_bt,mu_t);
plot3(rr_bt_start(1),rr_bt_start(2),rr_bt_start(3),'k+','LineWidth',2,'MarkerSize', 10);

% Parametri orbita trasferimento ap
r_apo_i = a_i * (1 + e_i);
r_peri_f = a_f * (1 - e_f);
a_bt = (r_apo_i + r_peri_f) / 2;
e_bt = abs(r_apo_i - r_peri_f) / (r_apo_i + r_peri_f);
i_bt = i_f; OM_bt = OM_f;

plotOrbit(a_bt, e_bt, i_bt, OM_bt, om_bt, pi, 2*pi, dth, mu_t, 'm');

plot3(rr_f(1),rr_f(2),rr_f(3), '.', 'MarkerSize', 40, 'LineWidth', 4, 'Color','m');
Terra3d(Rt);

% Legenda
legend('Orbita Iniziale', 'Orbita Finale', 'Posizione Iniziale Satellite', ...
      'Punto Manovra Cambio Argomento Pericentro', 'Orbita Cambio Argomento Pericentro', ...
      'Punto Manovra Cambio Piano', 'Orbita dopo Cambio Piano', ...
      'Punto Manovra Trasferimento Bitangente', 'Orbita Trasferimento Bitangente', 'Posizione Finale Satellite','Location','northeastoutside');

% Setup animazione T_MIN
h_text_time2 = annotation('textbox', pos_time, 'String', '', 'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
h_text_vel2 = annotation('textbox', pos_vel, 'String', '', 'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');

satellite2 = plot3(NaN, NaN, NaN, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 8, 'HandleVisibility', 'off');

% Video writer per T_MIN
v2 = VideoWriter('Satellite_Transfer_TMIN.mp4', 'MPEG-4');
v2.Quality = 100; v2.FrameRate = 120; open(v2);

% Esecuzione animazione T_MIN
time_accum = 0;
animateTransfer(a_i, e_i, i_i, OM_i, om_i, th_i, thi_CPA, mu_t, satellite2, v2, true, view_start, view_end, time_accum, h_text_time2, h_text_vel2);
time_accum = time_accum + DeltaT_CA;

animateTransfer(a_i, e_i, i_i, OM_i, omi, thf_CPA, PositiveAngles(theta_CP), mu_t, satellite2, v2, false, view_start, view_end, time_accum, h_text_time2, h_text_vel2);
time_accum = time_accum + DeltaT_CP;

animateTransfer(a_i, e_i, i_f, OM_f, om_cp_new, PositiveAngles(theta_CP), thi_bt, mu_t, satellite2, v2, false, view_start, view_end, time_accum, h_text_time2, h_text_vel2);
time_accum = time_accum + TOF(a_i,e_i,PositiveAngles(theta_CP),thi_bt,mu_t);

animateTransfer(a_bt, e_bt, i_bt, OM_bt, om_bt, thi_bt, thf_bt, mu_t, satellite2, v2, false, view_start, view_end, time_accum, h_text_time2, h_text_vel2);
time_accum = time_accum + DeltaT_BT;

animateTransfer(a_f, e_f, i_f, OM_f, om_bt, thf_bt, th_f, mu_t, satellite2, v2, false, view_start, view_end, time_accum, h_text_time2, h_text_vel2);

close(v2);
disp('Animazione T_MIN completata. Video salvato come "Satellite_Transfer_TMIN.mp4"');