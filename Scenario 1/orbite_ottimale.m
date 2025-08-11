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

% tra quelle presenti per il deltaV minore ora seleziono quella con il deltat minore

[~, idx_min_dv_dt] = min(Sv_num(:,5));

% plot

% Trova gli indici delle strategie con delta-v minimo e tra questi quello
% con delta-t minimo e punto di manovra più vicino a quello attuale

deltav_optimal_indices = find(dv_vec == dv_MIN & dt_vec == min(Sv_num(:,5)));
optimal_index = deltav_optimal_indices(1);

% Estraggo i parametri della strategia selezionata
params = parameters_all(optimal_index,:);
rb_val = params(1);
th_cp = params(14);
thi_ca = params(15);
thf_ca = params(16);
th_cs = params(17);
om_cp = params(12);
om_ca = params(13);









%% PLOT TRASFERIMENTO OTTIMALE V_MIN

% Configura il plot
figure;
hold on;
grid on;
axis equal;
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title(['Strategia Ottimale per Velocità Minima #' num2str(optimal_index)]);


plotOrbit(a_i,e_i,i_i,OM_i,om_i,0,2*pi,dth,mu_t, 'b'); % orbita iniziale in blu
plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,dth,mu_t, 'r'); % orbita finale in rosso


% Trasferimento bitangente aa
type = 'aa';
[DeltaV1_BT, DeltaV2_BT, DeltaT4, om_bt] = bitangentTransfer(a_i, e_i, a_f, e_f, om_i, type, mu_t);
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
% plot punto di partenza satellite
[rr_0] = par2car(a_i,e_i,i_i,OM_i,om_i,th_i,mu_t);
plot3(rr_0(1),rr_0(2),rr_0(3), '.', 'MarkerSize', 40, 'LineWidth', 4);
[rr_3] = par2car(a_i,e_i,i_i,OM_i,om_i,thi_bt,mu_t); % punto partenza manovra trasferimento bitangente
plot3(rr_3(1),rr_3(2),rr_3(3),'k+','LineWidth',2,'MarkerSize', 10);
% Calcolo perigeo e apogeo - CORRETTI
r_apo_i = a_i * (1 + e_i); % apogeo orbita iniziale
r_apo_f = a_f * (1 + e_f); % apogeo orbita finale
% Parametri orbita di trasferimento bitangente
a_bt = (r_apo_i + r_apo_f) / 2;
e_bt = abs(r_apo_f - r_apo_i) / (r_apo_f + r_apo_i);
% Stessa inclinazione, OM dell'orbita iniziale
i_bt = i_i;
OM_bt = OM_i;
% Plot orbita intermedia in magenta
plotOrbit(a_bt, e_bt, i_bt, OM_bt, om_i, pi, 2*pi, dth, mu_t, 'm');

% Plot orbita dopo trasferimento bitangente
plotOrbit(a_f, e_f, i_bt, OM_bt, om_bt, 0, 2*pi, dth, mu_t, 'black');


% Cambio del piano orbitale
[DeltaV_CP, om_2, theta_CP, DeltaT1] = changeOrbitalPlane(a_f, e_f, i_bt, OM_bt, om_bt, i_f, OM_f, thf_bt, mu_t);

plotOrbit(a_f,e_f,i_f,OM_f,om_2,0,2*pi,dth,mu_t, 'g');  % plot orbita di cambio piano


% plot punto di partenza manovra cambio di piano
[rr_1] = par2car(a_f,e_f,i_bt,OM_bt,om_bt,theta_CP,mu_t); % plot punto partenza manovra cambio piano
plot3(rr_1(1),rr_1(2),rr_1(3),'k+','LineWidth',2,'MarkerSize', 10);

% Cambio dell'anomalia del pericentro
[DeltaV_CA, om_3, thi_CPA, thf_CPA, DeltaT2] = changePericenterArg(a_f, e_f, om_2, om_f, theta_CP, 1, mu_t);

%DeltaT3 = TOF(a_i,e_i,thf_CPA,0,mu_t); % Tempo di volo per ulteriore trasferimento --> arrivo al punto del trasferimento successivo

[rr_2] = par2car(a_f,e_f,i_f,OM_f,om_2,thi_CPA,mu_t); % plot punto partenza manovra cambio argomento pericentro
plot3(rr_2(1),rr_2(2),rr_2(3),'k+','LineWidth',2,'MarkerSize', 20);

%plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,dth,mu_t, 'c'); % plot orbita di cambio argomento pericentro

% plot punto di arrivo satellite
plot3(rr_f(1),rr_f(2),rr_f(3), '.', 'MarkerSize', 40, 'LineWidth', 4,'Color','m');

Terra3d(Rt); % plot terra 3d

% Legenda
legend('Orbita Iniziale', 'Orbita Finale', 'Posizione Iniziale Satellite', ...
      'Punto Manovra Trasferimento Bitangente', 'Orbita Trasferimento Bitangente', ...
      'Orbita dopo Trasferimento Bitangente', 'Orbita dopo Cambio Piano', ...
      'Punto Manovra Cambio Piano', 'Punto Manovra Cambio Argomento Pericentro', 'Posizione Finale Satellite');

%% PLOT TRASFERIMENTO OTTIMALE T_MIN

% tra quelle presenti per il deltaT minore ora seleziono quella con il deltaV minore
[~, idx_min_dt] = min(St_num(:,5));
% plot
% Trova gli indici delle strategie con delta-t minimo e, tra questi, quello
% con delta-v minimo e punto di manovra più vicino a quello attuale
deltat_optimal_indices = find(dt_vec == dt_MIN & dv_vec == min(St_num(:,5)));
optimal_index = deltat_optimal_indices(1);
% Estraggo i parametri della strategia selezionata
params = parameters_all(optimal_index,:);
rb_val = params(1);
th_cp = params(14);
thi_ca = params(15);
thf_ca = params(16);
th_cs = params(17);
om_cp = params(12);
om_ca = params(13);

% inizio plot
figure;
hold on;
grid on;
axis equal;
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title(['Strategia Ottimale per Tempo Minimo #' num2str(optimal_index)]);

% Plot orbita iniziale (blu)
plotOrbit(a_i,e_i,i_i,OM_i,om_i,0,2*pi,dth,mu_t,'b');

% Plot orbita finale (rossa)
plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,dth,mu_t,'r');

% Plot posizione iniziale satellite
[rr_0] = par2car(a_i,e_i,i_i,OM_i,om_i,th_i,mu_t);
plot3(rr_0(1),rr_0(2),rr_0(3), '.', 'MarkerSize', 40, 'LineWidth', 4);

% 1. CAMBIO ASSI (Cambio dell'anomalia del pericentro)
[omi] = Inverse_changeplane(i_i,OM_i,om_f,i_f,OM_f);

[DeltaV_CA, om_ca_new, thi_CPA, thf_CPA, DeltaT_CA] = changePericenterArg(a_i, e_i, om_i, omi, th_i, 1, mu_t);

% Plot punto manovra cambio argomento pericentro
[rr_ca] = par2car(a_i,e_i,i_i,OM_i,om_i,thi_CPA,mu_t);
plot3(rr_ca(1),rr_ca(2),rr_ca(3),'k+','LineWidth',2,'MarkerSize', 10);

% Plot orbita dopo cambio argomento pericentro (ciano)
plotOrbit(a_i,e_i,i_i,OM_i,omi,0,2*pi,dth,mu_t, 'c');

% 2. CAMBIO PIANO ORBITALE
[DeltaV_CP, om_cp_new, theta_CP, DeltaT_CP] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, omi, i_f, OM_f, thf_CPA, mu_t);

% Plot punto manovra cambio piano
[rr_cp] = par2car(a_i,e_i,i_i,OM_i,omi,theta_CP,mu_t);
plot3(rr_cp(1),rr_cp(2),rr_cp(3),'k+','LineWidth',2,'MarkerSize', 10);

% Plot orbita dopo cambio piano (verde)
plotOrbit(a_i,e_i,i_f,OM_f,om_cp_new,0,2*pi,dth,mu_t, 'g');

% 3. TRASFERIMENTO BITANGENTE AP (apocentro-pericentro)
type = 'ap';
[DeltaV1_BT, DeltaV2_BT, DeltaT_BT, om_bt] = bitangentTransfer(a_i, e_i, a_f, e_f, om_cp_new, type, mu_t);

% Per trasferimento ap: si parte dall'apocentro (pi) e si arriva al pericentro (0)
thi_bt = pi;  % parte dall'apocentro orbita iniziale
thf_bt = 0;   % arriva al pericentro orbita finale

% Plot punto partenza trasferimento bitangente (apocentro orbita corrente)
[rr_bt_start] = par2car(a_i,e_i,i_f,OM_f,om_cp_new,thi_bt,mu_t);
plot3(rr_bt_start(1),rr_bt_start(2),rr_bt_start(3),'k+','LineWidth',2,'MarkerSize', 10);

% Calcolo parametri orbita di trasferimento bitangente ap
r_apo_i = a_i * (1 + e_i);  % apocentro orbita iniziale (punto di partenza)
r_peri_f = a_f * (1 - e_f); % pericentro orbita finale (punto di arrivo)

% Parametri orbita di trasferimento
a_bt = (r_apo_i + r_peri_f) / 2;
e_bt = abs(r_apo_i - r_peri_f) / (r_apo_i + r_peri_f);

% Stessa inclinazione e RAAN dell'orbita dopo cambio piano
i_bt = i_f;
OM_bt = OM_f;

% Plot orbita di trasferimento bitangente (magenta)
plotOrbit(a_bt, e_bt, i_bt, OM_bt, om_bt, 0, pi, dth, mu_t, 'm');

% Plot orbita finale raggiunta (nera)
%plotOrbit(a_f, e_f, i_f, OM_f, om_bt, 0, 2*pi, dth, mu_t, 'k');

% Plot punto di arrivo finale
plot3(rr_f(1),rr_f(2),rr_f(3), '.', 'MarkerSize', 40, 'LineWidth', 4, 'Color','m');

% Plot Terra
Terra3d(Rt);

% Legenda corretta
legend('Orbita Iniziale', 'Orbita Finale', 'Posizione Iniziale Satellite', ...
       'Punto Manovra Cambio Argomento Pericentro', 'Orbita dopo Cambio Argomento Pericentro', ...
       'Punto Manovra Cambio Piano', 'Orbita dopo Cambio Piano', ...
       'Punto Partenza Trasferimento Bitangente', 'Orbita Trasferimento Bitangente', ...
       'Orbita Finale Raggiunta', 'Posizione Finale Satellite', 'Terra', ...
       'Location', 'best');