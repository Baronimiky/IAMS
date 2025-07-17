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

% Conversione delle coordinate cartesiane in parametri orbitali --> vedi
% sopra

% PLOT ORBITA INIZIALE, ORBITA FINALE
figure('Position', [100, 100, 1000, 600])
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

Terra3d(Rt); % plot terra 3d

% PLOT LEGENDA
xlabel('Asse X [Km]')
ylabel('Asse Y [Km]')
zlabel('Asse Z [Km]')
legend('dati', 'Location', 'northwest')
legend('Orbita iniziale','Orbita finale','Punto di partenza satellite','Punto di arrivo satellite','Punto di partenza cambio di piano','Cambio piano','Punto partenza cambio argomento pericentro','Cambio argomento del pericentro','Punto partenza trasferimento bitangente');

%% Trasferimento biellittico bitangente pa

P = 1;
thi_bielliptic = 0;
thf_bielliptic = pi;
apocenter_radius = 1.09 * a_f * (1 + e_f);    % [km]

DeltaT4 = TOF(a_i, e_i, thf_bt, thi_bielliptic, mu_t);

[DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2, om_f] = biellipticTransfer(a_i, e_i, a_f, e_f, apocenter_radius, om_f, om_f, P, mu_t);

DeltaT5 = TOF(a_f, e_f, thf_bielliptic, th_f, mu_t);

% CALCOLO DELTA_V, DELTA_T
Delta_VTot_biellittico = abs(DeltaV_CP) + abs(DeltaV_CA) + abs(DeltaV1) + abs(DeltaV2) + abs(DeltaV3); % deltaV totale biellittico [m/s]
Delta_TTot_biellittico = DeltaT1 + DeltaT2 + DeltaT3 + DeltaT4 + DeltaT5 + Deltat1 + Deltat2; % deltaT totale biellittico [s]






%% Strategia ottimale e plot
[St_num,Sv_num,dt_MIN,dv_MIN,dt_vec,dv_vec,param_vec]=OptimalStrategy(a_i,e_i,i_i,OM_i,om_i,th_i,a_f,e_f,i_f,OM_f,om_f,th_f,mu_t);

% tra quelle presenti per il deltaV minore ora seleziono quella con il deltat minore

[~, idx_min_dv_dt] = min(Sv_num(:,5));

% plot

% Trova gli indici delle strategie con delta-v minimo e tra questi quello
% con delta-t minimo e punto di manovra più vicino a quello attuale

deltav_optimal_indices = find(dv_vec == dv_MIN & dt_vec == min(Sv_num(:,5)));
optimal_index = deltav_optimal_indices(1);

% Estraggo i parametri della strategia selezionata
params = param_vec(optimal_index,:);
rb_val = params(1);
th_cp = params(14);
thi_ca = params(15);
thf_ca = params(16);
th_cs = params(17);
om_cp = params(12);
om_ca = params(13);

% Configura il plot
figure;
hold on;
grid on;
axis equal;
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title(['Strategia Ottimale #' num2str(optimal_index)]);

% Plot orbita iniziale (blu)
plotOrbit(a_i,e_i,i_i,OM_i,om_i,0,2*pi,0.01,mu_t,'b');

% Plot orbita finale (rossa)
plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,0.01,mu_t,'r');

% Plot orbite di trasferimento intermedie (linee tratteggiate)
% Orbita dopo cambio piano
plotOrbit(a_i,e_i,i_f,OM_f,om_cp,0,2*pi,0.01,mu_t,'g:');

% Orbita dopo cambio pericentro
plotOrbit(a_i,e_i,i_f,OM_f,om_ca,0,2*pi,0.01,mu_t,'m:');

% Se è una strategia biellittica, plotto anche l'orbita intermedia
if rb_val > 0
    % Calcola parametri orbita biellittica intermedia
    a_intermediate = (a_i + rb_val) / 2;
    e_intermediate = abs(rb_val - a_i) / (rb_val + a_i);
    plotOrbit(a_intermediate,e_intermediate,i_f,OM_f,om_ca,0,2*pi,0.01,mu_t,'k:');
end

% Plot punti iniziale e finale
rr_start = par2car(a_i,e_i,i_i,OM_i,om_i,th_i,mu_t);
plot3(rr_start(1),rr_start(2),rr_start(3),'bo','MarkerSize',10,'MarkerFaceColor','b');
rr_end = par2car(a_f,e_f,i_f,OM_f,om_f,th_f,mu_t);
plot3(rr_end(1),rr_end(2),rr_end(3),'ro','MarkerSize',10,'MarkerFaceColor','r');

% Plot punti di manovra
rr_cp = par2car(a_i,e_i,i_i,OM_i,om_i,th_cp,mu_t);
plot3(rr_cp(1),rr_cp(2),rr_cp(3),'ks','MarkerSize',10,'LineWidth',2);
rr_ca = par2car(a_i,e_i,i_f,OM_f,om_cp,thf_ca,mu_t);
plot3(rr_ca(1),rr_ca(2),rr_ca(3),'kd','MarkerSize',10,'LineWidth',2);
rr_cs = par2car(a_i,e_i,i_f,OM_f,om_ca,th_cs,mu_t);
plot3(rr_cs(1),rr_cs(2),rr_cs(3),'k^','MarkerSize',10,'LineWidth',2);

% Legenda
if rb_val > 0
    legend('Orbita iniziale','Orbita finale',...
        'Orbita post-cambio piano','Orbita post-cambio pericentro','Orbita biellittica intermedia',...
        'Punto iniziale','Punto finale',...
        'Cambio piano','Cambio pericentro','Cambio forma',...
        'Location','best');
else
    legend('Orbita iniziale','Orbita finale',...
        'Orbita post-cambio piano','Orbita post-cambio pericentro',...
        'Punto iniziale','Punto finale',...
        'Cambio piano','Cambio pericentro','Cambio forma',...
        'Location','best');
end

Terra3d(Rt); % plot terra 3d

hold off;

% Stampo informazioni sulla strategia
fprintf('\n=== Strategia Ottimale #%d ===\n', optimal_index);
fprintf('Delta-V totale: %.4f km/s\n', dv_MIN);
fprintf('Tempo totale: %.4f s\n', dt_vec(optimal_index));
fprintf('Parametri manovre:\n');
fprintf(' - Cambio piano: th=%.2f rad, dv=%.4f km/s\n', PositiveAngles(th_cp), params(2));
fprintf(' - Cambio pericentro: th=%.2f rad, dv=%.4f km/s\n', thf_ca, params(3));
fprintf(' - Cambio forma: th=%.2f rad, dv=%.4f, %.4f, %.4f km/s\n', th_cs, params(4:6));
if rb_val > 0
    fprintf(' - Raggio biellittico: %.2f km\n', rb_val);
end
fprintf('Parametri: c1: %1d\tc2: %1d\t c3: %1d\t c4: %1d\n', Sv_num(idx_min_dv_dt, 1), Sv_num(idx_min_dv_dt, 2), Sv_num(idx_min_dv_dt, 3), Sv_num(idx_min_dv_dt, 4));


% delta t minimo

% tra quelle presenti per il deltaT minore ora seleziono quella con il deltaV minore
[~, idx_min_dt] = min(St_num(:,5));

% plot
% Trova gli indici delle strategie con delta-t minimo e, tra questi, quello
% con delta-v minimo e punto di manovra più vicino a quello attuale
deltat_optimal_indices = find(dt_vec == dt_MIN & dv_vec == min(St_num(:,5)));
optimal_index = deltat_optimal_indices(1);

% Estraggo i parametri della strategia selezionata
params = param_vec(optimal_index,:);
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
plotOrbit(a_i,e_i,i_i,OM_i,om_i,0,2*pi,0.01,mu_t,'b');

% Plot orbita finale (rossa)
plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,0.01,mu_t,'r');

% Plot orbite di trasferimento intermedie (linee tratteggiate)
% Orbita dopo cambio piano
plotOrbit(a_i,e_i,i_f,OM_f,om_cp,0,2*pi,0.01,mu_t,'g:');

% Orbita dopo cambio pericentro
plotOrbit(a_i,e_i,i_f,OM_f,om_ca,0,2*pi,0.01,mu_t,'m:');

% Se è una strategia biellittica, plotto anche l'orbita intermedia
if rb_val > 0
    % Calcola parametri orbita biellittica intermedia
    a_intermediate = (a_i + rb_val) / 2;
    e_intermediate = abs(rb_val - a_i) / (rb_val + a_i);
    plotOrbit(a_intermediate,e_intermediate,i_f,OM_f,om_ca,0,2*pi,0.01,mu_t,'k:');
end

% Plot punti iniziale e finale
rr_start = par2car(a_i,e_i,i_i,OM_i,om_i,th_i,mu_t);
plot3(rr_start(1),rr_start(2),rr_start(3),'bo','MarkerSize',10,'MarkerFaceColor','b');

rr_end = par2car(a_f,e_f,i_f,OM_f,om_f,th_f,mu_t);
plot3(rr_end(1),rr_end(2),rr_end(3),'ro','MarkerSize',10,'MarkerFaceColor','r');

% Plot punti di manovra
rr_cp = par2car(a_i,e_i,i_i,OM_i,om_i,th_cp,mu_t);
plot3(rr_cp(1),rr_cp(2),rr_cp(3),'ks','MarkerSize',10,'LineWidth',2);

rr_ca = par2car(a_i,e_i,i_f,OM_f,om_cp,thf_ca,mu_t);
plot3(rr_ca(1),rr_ca(2),rr_ca(3),'kd','MarkerSize',10,'LineWidth',2);

rr_cs = par2car(a_i,e_i,i_f,OM_f,om_ca,th_cs,mu_t);
plot3(rr_cs(1),rr_cs(2),rr_cs(3),'k^','MarkerSize',10,'LineWidth',2);

% Legenda
if rb_val > 0
    legend('Orbita iniziale','Orbita finale',...
           'Orbita post-cambio piano','Orbita post-cambio pericentro','Orbita biellittica intermedia',...
           'Punto iniziale','Punto finale',...
           'Cambio piano','Cambio pericentro','Cambio forma',...
           'Location','best');
else
    legend('Orbita iniziale','Orbita finale',...
           'Orbita post-cambio piano','Orbita post-cambio pericentro',...
           'Punto iniziale','Punto finale',...
           'Cambio piano','Cambio pericentro','Cambio forma',...
           'Location','best');
end

Terra3d(Rt); % plot terra 3d

hold off;

% Stampo informazioni sulla strategia
fprintf('\n=== Strategia Ottimale per Tempo Minimo #%d ===\n', optimal_index);
fprintf('Tempo totale: %.4f s\n', dt_MIN);
fprintf('Delta-V totale: %.4f km/s\n', dv_vec(optimal_index));
fprintf('Parametri manovre:\n');
fprintf(' - Cambio piano: th=%.2f rad, dv=%.4f km/s\n', th_cp, params(2));
fprintf(' - Cambio pericentro: th=%.2f rad, dv=%.4f km/s\n', thf_ca, params(3));
fprintf(' - Cambio forma: th=%.2f rad, dv=%.4f, %.4f, %.4f km/s\n', th_cs, params(4:6));
if rb_val > 0
    fprintf(' - Raggio biellittico: %.2f km\n', rb_val);
end
fprintf('Parametri: c1: %1d\tc2: %1d\t c3: %1d\t c4: %1d\n', St_num(idx_min_dt, 1), St_num(idx_min_dt, 2), St_num(idx_min_dt, 3), St_num(idx_min_dt, 4));

%% Trasferimento diretto con punti iniziali e finali fissati

delta_v_ott = 0;
delta_v_vect = [];
e_T_n_vec = [];
i = 0;


[r1,vv_i] = par2car(a_i,e_i,i_i,OM_i,om_i,th_i,mu_t);

[r2,vv_f] = par2car(a_f,e_f,i_f,OM_f,om_f,th_f,mu_t);

% Versore normale al piano orbitale di trasferimento
h_T = cross(r1, r2)/norm(cross(r1, r2));

% Calcolo dell'inclinazione dell'orbita di trasferimento
k = [0, 0, 1]; % Versore asse Z
i_T = acos(dot(h_T, k));

% Determinazione della linea dei nodi  
N = cross(k, h_T) / norm(cross(k, h_T));

% Calcolo della ascensione retta del nodo ascendente (RAAN)
i_hat = [1, 0, 0]; % Versore asse X
j_hat = [0, 1, 0]; % Versore asse Y

NT_i = dot(N, i_hat);
NT_j = dot(N, j_hat);

if NT_j >= 0
    Omega_T = acos(NT_i);
else
    Omega_T = 2*pi - acos(NT_i);
end

% Visualizzazione dei risultati
%disp(['Inclinazione i_T (rad): ', num2str(i_T)]);
%disp(['Ascensione retta del nodo ascendente Omega_T (rad): ', num2str(Omega_T)]);




% Creazione del ciclo for per ogni ω_T,n
for omega_T_n = 0:0.1:2*pi % Definire omega_values con i valori desiderati
    % Definizione delle matrici di rotazione
    R_omega_T_n = [cos(omega_T_n), sin(omega_T_n), 0;
                  -sin(omega_T_n), cos(omega_T_n), 0;
                   0, 0, 1]; % Rotazione intorno all'asse Z
    
    R_i_T = [1, 0, 0;
             0, cos(i_T), sin(i_T);
             0, -sin(i_T), cos(i_T)]; % Rotazione intorno all'asse X
    
    R_Omega_T = [cos(Omega_T), sin(Omega_T), 0;
                -sin(Omega_T), cos(Omega_T), 0;
                 0, 0, 1]; % Rotazione intorno all'asse Z

    % Matrice di trasformazione dal sistema inerziale al perifocale
    T_Elio_PF_T_n = R_omega_T_n * R_i_T * R_Omega_T;

    % Calcolo della posizione nei sistemi perifocali
    r1_T_n = T_Elio_PF_T_n * (r1);
    r2_T_n = T_Elio_PF_T_n * (r2);

    % Calcolo delle anomalie vere
    theta1_T_n = atan2(r1_T_n(2)/norm(r1), r1_T_n(1)/norm(r1));
    theta2_T_n = atan2(r2_T_n(2)/norm(r2), r2_T_n(1)/norm(r2));
    theta1_T_n = PositiveAngles(theta1_T_n);
    theta2_T_n = PositiveAngles(theta2_T_n);

    % Calcolo degli elementi orbitali
    e_T_n = (norm(r2) - norm(r1))/(norm(r1) * r1_T_n(1)/norm(r1) - norm(r2) * r2_T_n(1)/norm(r2));
    if e_T_n <=0 || e_T_n >= 1
        continue;
    end
        e_T_n_vec = [e_T_n_vec; e_T_n];
        %e_T_n = norm(e_T_n);
        a_T_n = norm(r1) * (1 + e_T_n * cos(theta1_T_n))/(1 - e_T_n^2);
    
        % Determinazione delle velocità nei punti di manovra
    
        %p_T_n = a_T_n * (1-e_T_n^2);
        %v1T_n = sqrt(mu_t./p_T_n) * [-sin(theta1_T_n); e_T_n+cos(theta1_T_n); 0];
        %v2T_n = sqrt(mu_t./p_T_n) * [-sin(theta2_T_n); e_T_n+cos(theta2_T_n); 0];
        [~,v1T_n] = par2car(a_T_n,e_T_n,i_T,Omega_T,omega_T_n,theta1_T_n,mu_t);
        [~,v2T_n] = par2car(a_T_n,e_T_n,i_T,Omega_T,omega_T_n,theta2_T_n,mu_t);
    
        % Calcolo del delta-v totale
        delta_v_n = norm(v1T_n - vv_i) + norm(vv_f - v2T_n);

        %calcolo tempo

        deltat = TOF(a_T_n,e_T_n,theta1_T_n,theta2_T_n,mu_t);
    
        if abs(delta_v_ott) > abs(delta_v_n) || delta_v_ott == 0
            delta_v_ott = delta_v_n;
            delta_v_ott1= norm(v1T_n - vv_i);
            delta_v_ott2 = norm(vv_f - v2T_n);
            th_i_ott = th_i;
            th_f_ott = th_f;
            omega_T_ott = omega_T_n;
            e_T_ott = e_T_n;
            a_T_ott = a_T_n;
            theta1_T_ott = theta1_T_n;
            theta2_T_ott = theta2_T_n;
            i_T_ott =i_T;
            Omega_T_ott = Omega_T;
            delta_t_ott = deltat;
        end
        i = i+1;
        delta_v_vect= [delta_v_vect; delta_v_n];
end

    % Visualizzazione dei risultati
            disp(['Per ω_T,n = ', num2str(omega_T_ott)]);
            disp(['  θ1_T,n (rad): ', num2str(theta1_T_ott)]);
            disp(['  θ2_T,n (rad): ', num2str(theta2_T_ott)]);
            disp(['  e_T,n: ', num2str(e_T_ott)]);
            disp(['  a_T,n: ', num2str(a_T_ott)]);
            disp(['  Δv_n: ', num2str(delta_v_ott)]);
            disp('----------------------------------');
  j = 1:1:i;
   plot(j,delta_v_vect);
   legend('\DeltaV');
   xlabel('Iterazioni')
   ylabel('\DeltaV [Km/s]')
   xlim([1,i])
grid on
%exportgraphics(gcf, 'delta_v_for.pdf', 'ContentType', 'vector', 'BackgroundColor','none')
    
%% Plot Orbite con trasferimento Grid Search
figure('Position', [100, 100, 1000, 600]);
hold on;
plotOrbit(a_i,e_i,i_i,OM_i,om_i,0,2*pi,dth,mu_t,'r');
plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,dth,mu_t,'b');
plotOrbit(a_T_ott,e_T_ott,i_T_ott,Omega_T_ott,omega_T_ott,theta1_T_ott,theta2_T_ott,dth,mu_t,'y');
[rr1,~] = par2car(a_T_ott,e_T_ott,i_T_ott,Omega_T_ott,omega_T_ott,theta1_T_ott,mu_t);
plot3(rr1(1),rr1(2),rr1(3), '.', 'MarkerSize', 20, 'LineWidth', 2);
[rr2] = par2car(a_T_ott,e_T_ott,i_T_ott,Omega_T_ott,omega_T_ott,theta2_T_ott,mu_t);
plot3(rr2(1),rr2(2),rr2(3), '.', 'MarkerSize', 20, 'LineWidth', 2);
Terra3d(Rt);
grid on
xlabel('Asse X [Km]')
ylabel('Asse Y [Km]')
zlabel('Asse Z [Km]')
view(-10,25)
legend('dati', 'Location', 'northeast')
legend('orbita geostazionaria','orbita parcheggio','trasferimento','punto iniziale','punto finale','Terra');

%% Trasferimento diretto scegliendo i punti iniziali e finali più favorevoli, e poi attendendo sull'orbita

delta_v_ott = 0;
delta_v_vect = [];
e_T_n_vec = [];
i = 0;

for th_i = 0:0.1:2*pi
    for th_f = 0:0.1:2*pi

        [r1,vv_i] = par2car(a_i,e_i,i_i,OM_i,om_i,th_i,mu_t);
    
        [r2,vv_f] = par2car(a_f,e_f,i_f,OM_f,om_f,th_f,mu_t);
    
        % Versore normale al piano orbitale di trasferimento
        h_T = cross(r1, r2)/norm(cross(r1, r2));
        
        % Calcolo dell'inclinazione dell'orbita di trasferimento
        k = [0, 0, 1]; % Versore asse Z
        i_T = acos(dot(h_T, k));
        
        % Determinazione della linea dei nodi  
        N = cross(k, h_T) / norm(cross(k, h_T));
        
        % Calcolo della ascensione retta del nodo ascendente (RAAN)
        i_hat = [1, 0, 0]; % Versore asse X
        j_hat = [0, 1, 0]; % Versore asse Y
        
        NT_i = dot(N, i_hat);
        NT_j = dot(N, j_hat);
        
        if NT_j >= 0
            Omega_T = acos(NT_i);
        else
            Omega_T = 2*pi - acos(NT_i);
        end
        
        % Visualizzazione dei risultati
        %disp(['Inclinazione i_T (rad): ', num2str(i_T)]);
        %disp(['Ascensione retta del nodo ascendente Omega_T (rad): ', num2str(Omega_T)]);
        
        
        % Creazione del ciclo for per ogni ω_T,n
        for omega_T_n = 0:0.1:2*pi % Definire omega_values con i valori desiderati
            % Definizione delle matrici di rotazione
            R_omega_T_n = [cos(omega_T_n), sin(omega_T_n), 0;
                          -sin(omega_T_n), cos(omega_T_n), 0;
                           0, 0, 1]; % Rotazione intorno all'asse Z
            
            R_i_T = [1, 0, 0;
                     0, cos(i_T), sin(i_T);
                     0, -sin(i_T), cos(i_T)]; % Rotazione intorno all'asse X
            
            R_Omega_T = [cos(Omega_T), sin(Omega_T), 0;
                        -sin(Omega_T), cos(Omega_T), 0;
                         0, 0, 1]; % Rotazione intorno all'asse Z
        
            % Matrice di trasformazione dal sistema inerziale al perifocale
            T_Elio_PF_T_n = R_omega_T_n * R_i_T * R_Omega_T;
        
            % Calcolo della posizione nei sistemi perifocali
            r1_T_n = T_Elio_PF_T_n * (r1);
            r2_T_n = T_Elio_PF_T_n * (r2);
        
            % Calcolo delle anomalie vere
            theta1_T_n = atan2(r1_T_n(2)/norm(r1), r1_T_n(1)/norm(r1));
            theta2_T_n = atan2(r2_T_n(2)/norm(r2), r2_T_n(1)/norm(r2));
            theta1_T_n = PositiveAngles(theta1_T_n);
            theta2_T_n = PositiveAngles(theta2_T_n);

            % Calcolo degli elementi orbitali
            e_T_n = (norm(r2) - norm(r1))/(norm(r1) * r1_T_n(1)/norm(r1) - norm(r2) * r2_T_n(1)/norm(r2));
            if e_T_n <=0 || e_T_n >= 1
                continue;
            end
                e_T_n_vec = [e_T_n_vec; e_T_n];
                %e_T_n = norm(e_T_n);
                a_T_n = norm(r1) * (1 + e_T_n * cos(theta1_T_n))/(1 - e_T_n^2);
            
                % Determinazione delle velocità nei punti di manovra
            
                %p_T_n = a_T_n * (1-e_T_n^2);
                %v1T_n = sqrt(mu_t./p_T_n) * [-sin(theta1_T_n); e_T_n+cos(theta1_T_n); 0];
                %v2T_n = sqrt(mu_t./p_T_n) * [-sin(theta2_T_n); e_T_n+cos(theta2_T_n); 0];
                [~,v1T_n] = par2car(a_T_n,e_T_n,i_T,Omega_T,omega_T_n,theta1_T_n,mu_t);
                [~,v2T_n] = par2car(a_T_n,e_T_n,i_T,Omega_T,omega_T_n,theta2_T_n,mu_t);
            
                % Calcolo del delta-v totale
                delta_v_n = norm(v1T_n - vv_i) + norm(vv_f - v2T_n);

                %calcolo tempo

                deltat = TOF(a_T_n,e_T_n,theta1_T_n,theta2_T_n,mu_t);
            
                if abs(delta_v_ott) > abs(delta_v_n) || delta_v_ott == 0
                    delta_v_ott = delta_v_n;
                    delta_v_ott1= norm(v1T_n - vv_i);
                    delta_v_ott2 = norm(vv_f - v2T_n);
                    th_i_ott = th_i;
                    th_f_ott = th_f;
                    omega_T_ott = omega_T_n;
                    e_T_ott = e_T_n;
                    a_T_ott = a_T_n;
                    theta1_T_ott = theta1_T_n;
                    theta2_T_ott = theta2_T_n;
                    i_T_ott =i_T;
                    Omega_T_ott = Omega_T;
                    delta_t_ott = deltat;
                end
                i = i+1;
                delta_v_vect= [delta_v_vect; delta_v_n];
        end
    end
end

    % Visualizzazione dei risultati
            disp(['Per ω_T,n = ', num2str(omega_T_ott)]);
            disp(['  θ1_T,n (rad): ', num2str(theta1_T_ott)]);
            disp(['  θ2_T,n (rad): ', num2str(theta2_T_ott)]);
            disp(['  e_T,n: ', num2str(e_T_ott)]);
            disp(['  a_T,n: ', num2str(a_T_ott)]);
            disp(['  Δv_n: ', num2str(delta_v_ott)]);
            disp('----------------------------------');
  j = 1:1:i;
   plot(j,delta_v_vect);
   legend('\DeltaV');
   xlabel('Iterazioni')
   ylabel('\DeltaV [Km/s]')
   xlim([1,i])
grid on


%% Plot Orbite con trasferimento Grid Search
figure('Position', [100, 100, 1000, 600]);
hold on
plotOrbit(a_i,e_i,i_i,OM_i,om_i,0,2*pi,dth,mu_t,'r');
plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,dth,mu_t,'b');
plotOrbit( a_T_ott,e_T_ott,i_T_ott,Omega_T_ott,omega_T_ott,theta1_T_ott,theta2_T_ott+2*pi,dth,mu_t,'y');
[rr1,~] = par2car(a_T_ott,e_T_ott,i_T_ott,Omega_T_ott,omega_T_ott,theta1_T_ott,mu_t);
plot3(rr1(1),rr1(2),rr1(3), '.', 'MarkerSize', 20, 'LineWidth', 2);
[rr2] = par2car(a_T_ott,e_T_ott,i_T_ott,Omega_T_ott,omega_T_ott,theta2_T_ott,mu_t);
plot3(rr2(1),rr2(2),rr2(3), '.', 'MarkerSize', 20, 'LineWidth', 2);
Terra3d(Rt);
grid on
xlabel('Asse X [Km]')
ylabel('Asse Y [Km]')
zlabel('Asse Z [Km]')
view(-10,25)
legend('dati', 'Location', 'northeast')
legend('orbita geostazionaria','orbita parcheggio','trasferimento','punto iniziale','punto finale','Terra');