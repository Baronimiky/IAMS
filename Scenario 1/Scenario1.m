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

%[rr_i,vv_i] = par2car(a_i,e_i,i_i,OM_i,om_i,th_i,mu_t);

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

% Conversione delle coordinate cartesiane in parametri orbitali
[a_f, e_f, i_f, OM_f, om_f, th_f] = car2par(rr_f,vv_f,mu_t);

%% Cambio del piano orbitale
[DeltaV_CP, om_2, theta_CP] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu_t);
DeltaT1 = TOF(a_i, e_i, th_i, theta_CP, mu_t); % tempo manovra

% plot punto di partenza manovra
[rr_1] = par2car(a_i,e_i,i_i,OM_i,om_i,theta_CP,mu_t); % punto partenza manovra
plot3(rr_1(1),rr_1(2),rr_1(3),'k+','LineWidth',2);

%% Cambio dell'anomalia del pericentro
[DeltaV_CA, thi_CPA, thf_CPA] = changePericenterArg(a_i, e_i, om_2, om_f, mu_t);

% Seleziono theta opportuni
Deltat2_1 = TOF(a_i,e_i,theta_CP,thi_CPA(1),mu_t); % tempo manovra
Deltat2_2 = TOF(a_i,e_i,theta_CP,thi_CPA(2),mu_t); % tempo manovra

if Deltat2_1 > Deltat2_2
    DeltaT2 = Deltat2_2;
    DeltaT3 = TOF(a_i,e_i,thf_CPA(2),0,mu_t); % Tempo di volo per ulteriore trasferimento
    index = 2;
else
    DeltaT2 = Deltat2_1;
    DeltaT3 = TOF(a_i,e_i,thf_CPA(1),0,mu_t); % Tempo di volo per ulteriore trasferimento
    index = 1;
end

[rr_2] = par2car(a_i,e_i,i_i,OM_f,om_2,thi_CPA(index),mu_t); % punto partenza manovra
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

[DeltaV1_BT, DeltaV2_BT, DeltaT4] = bitangentTransfer(a_i, e_i, a_f, e_f, type, mu_t);

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
% [rr_3] = par2car(a_i,e_i,i_i,OM_f,om_f,thi_bt,mu_t); % punto partenza manovra
% plot3(rr_3(1),rr_3(2),rr_3(3),'k+','LineWidth',2);

%% Trasferimento biellittico bitangente

[DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2] = biellipticTransfer(a_i, e_i, a_f, e_f, Rt, mu_t);

DeltaT6 = TOF(a_f, e_f, pi, th_f, mu_t);

%% PLOT CAMBIO PIANO E PERICENTRO

plotOrbit(a_i,e_i,i_f,OM_f,om_2,0,2*pi,dth,mu_t, 'g');  % cambio piano
plotOrbit(a_i,e_i,i_f,OM_f,om_f,0,2*pi,dth,mu_t, 'c'); % cambio pericentro

%% PLOT LEGENDA E CALCOLO DELTA_V, DELTA_T

plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,dth,mu_t) % 
legend('Orbita iniziale','Orbita finale','Cambio piano','Cambio angolo','Trasferimento bitangente');

Delta_VTot = abs(DeltaV_CP) + abs(DeltaV_CA) + abs(DeltaV1_BT) + abs(DeltaV2_BT) % deltav totale [m/s]
Delta_TTot = DeltaT1 + DeltaT2 + DeltaT3 + DeltaT4 + DeltaT5 % delta Tempo totale [s]


%% CERCO LA COMBINAZIONE PIU' EFFICIENTE DI TRASFERIMENTI

sequenze = {
    {@changeOrbitalPlane, @changePericenterArg, @bitangentTransfer};       % piano → pericentro → bitangente
    {@changeOrbitalPlane, @changePericenterArg, @biellipticTransfer};      % piano → pericentro → biellittico
    {@bitangentTransfer, @changeOrbitalPlane, @changePericenterArg};       % Trasf. → piano → pericentro
    {@biellipticTransfer, @changeOrbitalPlane, @changePericenterArg};      % Biellittico → piano → pericentro
};

bitangent_types = {'pa', 'ap', 'pp', 'aa'};

ra_t_values = linspace(50000, 120000, 10);  % range apocentro intermedio biellittico
theta_vals = linspace(0, 2*pi, 20);         % possibile theta di esecuzione manovre

bestDeltaV = inf;
bestIndex = 0;

for ra_t = ra_t_values
    for s = 1:length(sequenze)
        deltaVtot = 0;
        a = a_i; e = e_i; i = i_i; OM = OM_i; om = om_i; th = th_i;

        for m = 1:length(sequenze{s})
            manovra = sequenze{s}{m};

            if isequal(manovra, @bitangentTransfer)
                bestDeltaVbit = inf;
                for t = 1:length(bitangent_types)
                    type = bitangent_types{t};
                    [dV1, dV2, ~] = bitangentTransfer(a, e, a_f, e_f, type, mu_t);
                    deltaVcandidate = abs(dV1) + abs(dV2);
                    if deltaVcandidate < bestDeltaVbit
                        bestDeltaVbit = deltaVcandidate;
                        bestType = type;
                    end
                end
                deltaVtot = deltaVtot + bestDeltaVbit;
                a = a_f; e = e_f;

            elseif isequal(manovra, @biellipticTransfer)
                [dV1, dV2, dV3, ~, ~] = biellipticTransfer(a, e, a_f, e_f, ra_t, mu_t);
                deltaVtot = deltaVtot + abs(dV1) + abs(dV2) + abs(dV3);
                a = a_f; e = e_f;

            elseif isequal(manovra, @changeOrbitalPlane)
                [dV, ~, ~] = changeOrbitalPlane(a, e, i, OM, om, i_f, OM_f, th, mu_t);
                deltaVtot = deltaVtot + abs(dV);
                i = i_f; OM = OM_f;

            elseif isequal(manovra, @changePericenterArg)
                [dV, ~, ~] = changePericenterArg(a, e, om, om_f, mu_t);
                deltaVtot = deltaVtot + abs(dV);
                om = om_f;
            end
        end

        if deltaVtot < bestDeltaV
            bestDeltaV = deltaVtot;
            bestIndex = s;
            bestConfig.ra_t = ra_t;
        end
    end
end

fprintf("Sequenza ottima: #%d con DeltaV totale = %.6f km/s\n", bestIndex, bestDeltaV);

%%

DeltaV_bench = benchmarkDeltaV(a_i, e_i, i_i, OM_i, om_i, a_f, e_f, i_f, OM_f, om_f, mu_t);
