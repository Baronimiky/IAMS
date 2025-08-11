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

%% STRATEGIA STANDARD

% PLOT ORBITA INIZIALE, ORBITA FINALE
fig = figure('Position', [100, 100, 1000, 600]);
hold on;
plotOrbit(a_i,e_i,i_i,OM_i,om_i,0,2*pi,dth,mu_t, 'b'); % orbita iniziale in blu
plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,dth,mu_t, 'r'); % orbita finale in rosso

% Cambio del piano orbitale
[DeltaV_CP, om_2, theta_CP, DeltaT1] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, th_i, mu_t);

% plot punto di partenza satellite
[rr_0] = par2car(a_i,e_i,i_i,OM_i,om_i,th_i,mu_t); % punto partenza manovra
plot3(rr_0(1),rr_0(2),rr_0(3), '.', 'MarkerSize', 40, 'LineWidth', 4);

% plot punto di arrivo satellite
plot3(rr_f(1),rr_f(2),rr_f(3), '.', 'MarkerSize', 40, 'LineWidth', 4);

% plot punto di partenza manovra cambio di piano
[rr_1] = par2car(a_i,e_i,i_i,OM_i,om_i,theta_CP,mu_t); % plot punto partenza manovra cambio piano
plot3(rr_1(1),rr_1(2),rr_1(3),'k+','LineWidth',2,'MarkerSize', 10);

plotOrbit(a_i,e_i,i_f,OM_f,om_2,0,2*pi,dth,mu_t, 'g');  % plot orbita di cambio piano

% Cambio dell'anomalia del pericentro
[DeltaV_CA, om_3, thi_CPA, thf_CPA, DeltaT2] = changePericenterArg(a_i, e_i, om_2, om_f, theta_CP, 1, mu_t);

DeltaT3 = TOF(a_i,e_i,thf_CPA,0,mu_t); % Tempo di volo per ulteriore trasferimento --> arrivo al punto del trasferimento successivo

[rr_2] = par2car(a_i,e_i,i_f,OM_f,om_2,thi_CPA,mu_t); % plot punto partenza manovra cambio argomento pericentro
plot3(rr_2(1),rr_2(2),rr_2(3),'k+','LineWidth',2,'MarkerSize', 10);

plotOrbit(a_i,e_i,i_f,OM_f,om_f,0,2*pi,dth,mu_t, 'c'); % plot orbita di cambio argomento pericentro

% Trasferimento bitangente pa

type = 'pa';

[DeltaV1_BT, DeltaV2_BT, DeltaT4, om_f] = bitangentTransfer(a_i, e_i, a_f, e_f, om_f, type, mu_t);

switch type
    case 'pa'
        thi_bt = 0;
        thf_bt = pi;
    case 'ap'
        thi_bt = pi;
        thf_bt = 0;
    case 'pp'
        thi_bt = 0;
        thf_bt = 0;
    case 'aa'
        thi_bt = pi;
        thf_bt = pi;
    otherwise
        error('Tipo di trasferimento bitangente non riconosciuto');
end

% Ultimo tempo di volo, dal punto di arrivo al punto finale
DeltaT5 = TOF(a_f, e_f, thf_bt, th_f, mu_t);

[rr_3] = par2car(a_i,e_i,i_f,OM_f,om_f,thi_bt,mu_t); % punto partenza manovra trasferimento bitangente
plot3(rr_3(1),rr_3(2),rr_3(3),'k+','LineWidth',2,'MarkerSize', 10);

% CALCOLO DELTA_V, DELTA_T
Delta_VTot_bitangente = abs(DeltaV_CP) + abs(DeltaV_CA) + abs(DeltaV1_BT) + abs(DeltaV2_BT); % deltaV totale bitangente [m/s]
Delta_TTot_bitangente = DeltaT1 + DeltaT2 + DeltaT3 + DeltaT4 + DeltaT5; % deltaT totale bitangente [s]

% Calcolo perigeo e apogeo
r_peri = a_i * (1 - e_i);   % perigeo orbita iniziale
r_apo  = a_f * (1 + e_f);   % apogeo orbita finale

% Parametri orbita di trasferimento bitangente
a_bt = (r_peri + r_apo) / 2;
e_bt = (r_apo - r_peri) / (r_apo + r_peri);

% Stessa inclinazione, OM e om dell’orbita finale (assunzione valida)
i_bt  = i_f;
OM_bt = OM_f;
om_bt = om_f;

% Plot orbita intermedia in magenta
plotOrbit(a_bt, e_bt, i_bt, OM_bt, om_bt, 0, 2*pi, dth, mu_t, 'm');

Terra3d(Rt); % plot terra 3d

% PLOT LEGENDA
xlabel('Asse X [Km]')
ylabel('Asse Y [Km]')
zlabel('Asse Z [Km]')
legend('Orbita iniziale','Orbita finale','Punto di partenza satellite','Punto di arrivo satellite','Punto di partenza cambio di piano','Cambio piano','Punto partenza cambio argomento pericentro','Cambio argomento del pericentro','Punto partenza trasferimento bitangente','Orbita Intermedia Tr. Bitangente','Location', 'northeastoutside');

% CREA I TEXTBOX FUORI DAGLI ASSI, SOTTO LA LEGENDA
pos_time = [0.82 0.13 0.15 0.05]; % [left bottom width height], normalized fig units
pos_vel  = [0.82 0.07 0.15 0.05];

h_text_time = annotation('textbox', pos_time, 'String', '', ...
    'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');

h_text_vel = annotation('textbox', pos_vel, 'String', '', ...
    'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');

% === ANIMAZIONE DEL SATELLITE ===
satellite = plot3(NaN, NaN, NaN, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 8, 'HandleVisibility', 'off');

% Imposta il video writer
v = VideoWriter('Satellite_Transfer.mp4', 'MPEG-4');
v.Quality = 100;
v.FrameRate = 60;
open(v);

% Parametri vista
view_start = [0, 20];
view_end   = [94, 40];

% Funzione per animare il movimento del satellite su una tratta
function animateTransfer(a,e,i,OM,om,th_start,th_end,mu,satellite,v,doRotate, view_start, view_end, time_offset, h_text_time, h_text_vel)
    dth_anim = 0.01;
    
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


% ESECUZIONE ANIMAZIONE

time_accum = 0; % tempo totale

% 1. Orbita iniziale fino al cambio piano (th_i -> theta_CP)
animateTransfer(a_i, e_i, i_i, OM_i, om_i, th_i, theta_CP, mu_t, satellite, v, true, view_start, view_end, time_accum, h_text_time, h_text_vel);
time_accum = time_accum + DeltaT1;

% 2. Cambio piano: da theta_CP a thi_CPA
animateTransfer(a_i, e_i, i_f, OM_f, om_2, theta_CP, thi_CPA, mu_t, satellite, v, false, view_start, view_end, time_accum, h_text_time, h_text_vel);
time_accum = time_accum + DeltaT2;

% 3. Cambio pericentro: da thi_CPA a thf_CPA
animateTransfer(a_i, e_i, i_f, OM_f, om_f, thf_CPA, thi_bt, mu_t, satellite, v, false, view_start, view_end, time_accum, h_text_time, h_text_vel);
time_accum = time_accum + DeltaT3;

% 4. Bitangente transfer: da thi_bt -> thf_bt su orbita intermedia
animateTransfer(a_bt, e_bt, i_bt, OM_bt, om_bt, thi_bt, thf_bt, mu_t, satellite, v, false, view_start, view_end, time_accum, h_text_time, h_text_vel);
time_accum = time_accum + DeltaT4;

% 5. Ultimo tratto su orbita finale: da thf_bt -> th_f (punto finale)
animateTransfer(a_f, e_f, i_f, OM_f, om_f, thf_bt, th_f, mu_t, satellite, v, false, view_start, view_end, time_accum, h_text_time, h_text_vel);
time_accum = time_accum + DeltaT5;

set(h_text_time, 'String', sprintf('Tempo: %.1f s', 72080.0));

% Forza il rendering
drawnow limitrate nocallbacks;  % Più efficiente per video
frame = getframe(gcf);
pause(0.1);  % 100ms per assicurare il rendering
writeVideo(v, frame);
drawnow limitrate nocallbacks;  % Più efficiente per video
frame = getframe(gcf);
pause(0.1);  % 100ms per assicurare il rendering
writeVideo(v, frame);
drawnow limitrate nocallbacks;  % Più efficiente per video
frame = getframe(gcf);
pause(0.1);  % 100ms per assicurare il rendering
writeVideo(v, frame);
drawnow limitrate nocallbacks;  % Più efficiente per video
frame = getframe(gcf);
pause(0.1);  % 100ms per assicurare il rendering
writeVideo(v, frame);
drawnow limitrate nocallbacks;  % Più efficiente per video
frame = getframe(gcf);
pause(0.1);  % 100ms per assicurare il rendering
writeVideo(v, frame);

% Chiudi video
close(v);
disp('Animazione completata. Video salvato come "Satellite_Transfer.mp4"');