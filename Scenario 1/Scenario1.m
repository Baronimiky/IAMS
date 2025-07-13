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
om_i = Or_i(1,5); % anomalia del pericentro[rad] w
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

[rr_i,vv_i] = par2car(a_i,e_i,i_i,OM_i,om_i,th_i,mu_t);

% Conversione delle coordinate cartesiane in parametri orbitali
[a_f,e_f,i_f,OM_f,om_f,th_f] = car2par(rr_f,vv_f,mu_t);

dth = 0.01; % discriminante theta
th0 = 0; % theta iniziale
thf = 2*pi; %theta finale

Rt = 6378; % raggio terra

%% PLOT TERRA, ORBITA INIZIALE, ORBITA FINALE

Terra3d(Rt); % plot terra 3d
plotOrbit(a_i,e_i,i_i,OM_i,om_i,th0,thf,dth,mu_t, 'b'); % orbita iniziale in blu
plotOrbit(a_f,e_f,i_f,OM_f,om_f,th0,thf,dth,mu_t, 'r'); % orbita finale in rosso


%% STRATEGIA STANDARD

%% Cambio del piano orbitale
[DeltaV_CP, om_2, theta_CP] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, th_i, mu_t);
DeltaT1 = TOF(a_i, e_i, th_i, theta_CP, mu_t); % tempo manovra

% plot punto di partenza manovra
[rr_1] = par2car(a_i,e_i,i_i,OM_i,om_i,theta_CP,mu_t); % punto partenza manovra
plot3(rr_1(1),rr_1(2),rr_1(3),'k+','LineWidth',2);

%% Cambio dell'anomalia del pericentro
[DeltaV_CA, om_3, thi_CPA, thf_CPA, DeltaT2] = changePericenterArg(a_i, e_i, om_2, om_f, theta_CP, 1, mu_t);

DeltaT3 = TOF(a_i,e_i,thf_CPA,0,mu_t); % Tempo di volo per ulteriore trasferimento

[rr_2] = par2car(a_i,e_i,i_i,OM_f,om_3,thi_CPA,mu_t); % punto partenza manovra
plot3(rr_2(1),rr_2(2),rr_2(3),'k+','LineWidth',2);

%% Trasferimento bitangente

% Identifica tipo manovra 'pa', 'ap', ecc.
if e_i > e_f
    type = 'pa';
elseif e_i < e_f
    type = 'ap';
elseif e_i == 0 && e_f == 0
    type = 'aa';
else
    type = 'pp';
end

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
        error('Tipo di trasferimento bitangente non riconosciuto.');
end

% Ultimo tempo di volo
DeltaT5 = TOF(a_f, e_f, thi_bt, thf_bt, mu_t);


% non mi funziona plot punto partenza manovra
[rr_3] = par2car(a_i,e_i,i_i,OM_f,om_f,thi_bt,mu_t); % punto partenza manovra
plot3(rr_3(1),rr_3(2),rr_3(3),'k+','LineWidth',2);

%% Trasferimento biellittico bitangente

%[DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2] = biellipticTransfer(a_i, e_i, a_f, e_f, apocenter_radius, mu_t);

%DeltaT6 = TOF(a_f, e_f, pi, th_f, mu_t);

%% PLOT CAMBIO PIANO E PERICENTRO

plotOrbit(a_i,e_i,i_f,OM_f,om_2,0,2*pi,dth,mu_t, 'g');  % cambio piano
plotOrbit(a_i,e_i,i_f,OM_f,om_f,0,2*pi,dth,mu_t, 'c'); % cambio pericentro

%% PLOT LEGENDA E CALCOLO DELTA_V, DELTA_T

plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,dth,mu_t) % 
legend('Orbita iniziale','Orbita finale','Cambio piano','Cambio angolo','Trasferimento bitangente');

Delta_VTot = abs(DeltaV_CP) + abs(DeltaV_CA) + abs(DeltaV1_BT) + abs(DeltaV2_BT); % deltav totale [m/s]
Delta_TTot = DeltaT1 + DeltaT2 + DeltaT3 + DeltaT4 + DeltaT5; % delta Tempo totale [s]






%%

[St_num,Sv_num,dt_MIN,dv_MIN,dt_vec,dv_vec,param_vec]=StrategyTestSTD_Optimized(a_i,e_i,i_i,OM_i,om_i,th_i,a_f,e_f,i_f,OM_f,om_f,th_f);
dv_MIN



% tra quelle presenti ora seleziono quella con il deltat minore



%% plot

% Trova l'indice della strategia con delta-v minimo
[~, idx_min] = min(dv_vec);

% Estrai i parametri rilevanti
param = param_vec(idx_min, :);

ra_t      = param(1);  % raggio apocentro per biellittico (solo se applicabile)
dv1       = param(4);  % primo impulso
dv2       = param(5);  % secondo impulso
dv3       = param(6);  % terzo impulso (solo biellittico)
om_plane  = param(13); % argomento pericentro dopo prima manovra
om_peri   = param(14); % argomento pericentro dopo seconda manovra
th_plane  = param(15); % anomalia vera manovra nel piano
th_shape  = param(17); % anomalia vera manovra di forma (cambio forma orbitale)

% Controlla se è un trasferimento biellittico o bitangente
isBiellittico = (dv3 > 0);  % se c'è il terzo impulso, è biellittico

% === Plot orbita iniziale ===
plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, dth, mu_t, 'g');  % verde
hold on

% === Plot orbita finale ===
plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, dth, mu_t, 'r');  % rossa

% === Punto partenza ===
[rr_start] = par2car(a_i,e_i,i_i,OM_i,om_i,th_plane,mu_t);
plot3(rr_start(1), rr_start(2), rr_start(3), 'k+','LineWidth',2);

% === Orbite di trasferimento ===
if isBiellittico
    % Primo arco di trasferimento: da orbita iniziale a ra_t
    rp1 = a_i * (1 - e_i);
    a_t1 = (rp1 + ra_t)/2;
    e_t1 = (ra_t - rp1)/(ra_t + rp1);
    plotOrbit(a_t1, e_t1, i_f, OM_f, om_plane, 0, pi, dth, mu_t, 'm');  % viola

    % Secondo arco: da ra_t a orbita finale
    rp2 = a_f * (1 - e_f);
    a_t2 = (rp2 + ra_t)/2;
    e_t2 = (ra_t - rp2)/(ra_t + rp2);
    plotOrbit(a_t2, e_t2, i_f, OM_f, om_peri, pi, 2*pi, dth, mu_t, 'c');  % ciano

else
    % Bitangente: una sola orbita di trasferimento
    % Calcola rpT e raT dal tipo (usa P se serve)
    rpT = a_i * (1 - e_i);
    raT = a_f * (1 + e_f);
    a_T = (rpT + raT)/2;
    e_T = (raT - rpT)/(raT + rpT);

    plotOrbit(a_T, e_T, i_f, OM_f, om_plane, 0, 2*pi, dth, mu_t, 'b');  % blu
end

% Legenda
legend('Orbita Iniziale','Orbita Finale','Punto manovra','Trasferimento 1','Trasferimento 2');
axis equal
grid on
title(sprintf('Strategia a Delta-v minimo (%.3f km/s)', dv_MIN))














%% Trasferimento biellittico bitangente

[DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2] = biellipticTransfer(a_i, e_i, a_f, e_f, apocenter_radius, mu_t);

%% Cambio del piano orbitale
[DeltaV_CP, om_2, theta_CP] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, th_i, mu_t);

%% Cambio dell'anomalia del pericentro
[DeltaV_CA, thi_CPA, thf_CPA] = changePericenterArg(a_i, e_i, om_2, om_f, mu_t);

deltav_tot = DeltaV1 + DeltaV2 + DeltaV3 + DeltaV_CP + DeltaV_CA