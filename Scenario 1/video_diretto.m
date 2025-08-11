%% ANIMAZIONE TRASFERIMENTO DIRETTO - CODICE COMPLETO
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

% %% CALCOLO TRAIETTORIA DIRETTA
% [r1,vv_i] = par2car(a_i,e_i,i_i,OM_i,om_i,th_i,mu_t);
% [r2,vv_f] = par2car(a_f,e_f,i_f,OM_f,om_f,th_f,mu_t);
% 
% % Versore normale al piano orbitale di trasferimento
% h_T = cross(r1, r2)/norm(cross(r1, r2));
% 
% % Calcolo dell'inclinazione dell'orbita di trasferimento
% k = [0, 0, 1]; % Versore asse Z
% i_T = acos(dot(h_T, k));
% 
% % Determinazione della linea dei nodi  
% N = cross(k, h_T) / norm(cross(k, h_T));
% 
% % Calcolo della ascensione retta del nodo ascendente (RAAN)
% i_hat = [1, 0, 0]; % Versore asse X
% j_hat = [0, 1, 0]; % Versore asse Y
% 
% NT_i = dot(N, i_hat);
% NT_j = dot(N, j_hat);
% 
% if NT_j >= 0
%     Omega_T = acos(NT_i);
% else
%     Omega_T = 2*pi - acos(NT_i);
% end
% 
% % Inizializzazione variabili per ottimizzazione
% delta_v_ott = 0;
% delta_v_vect = [];
% e_T_n_vec = [];
% i = 0;
% 
% % Loop di ottimizzazione
% for omega_T_n = 0:0.1:2*pi
%     % Matrici di rotazione
%     R_omega_T_n = [cos(omega_T_n), sin(omega_T_n), 0;
%                   -sin(omega_T_n), cos(omega_T_n), 0;
%                    0, 0, 1];
% 
%     R_i_T = [1, 0, 0;
%              0, cos(i_T), sin(i_T);
%              0, -sin(i_T), cos(i_T)];
% 
%     R_Omega_T = [cos(Omega_T), sin(Omega_T), 0;
%                 -sin(Omega_T), cos(Omega_T), 0;
%                  0, 0, 1];
% 
%     % Trasformazione al sistema perifocale
%     T_Elio_PF_T_n = R_omega_T_n * R_i_T * R_Omega_T;
%     r1_T_n = T_Elio_PF_T_n * (r1);
%     r2_T_n = T_Elio_PF_T_n * (r2);
% 
%     % Calcolo anomalie vere
%     theta1_T_n = atan2(r1_T_n(2)/norm(r1), r1_T_n(1)/norm(r1));
%     theta2_T_n = atan2(r2_T_n(2)/norm(r2), r2_T_n(1)/norm(r2));
%     theta1_T_n = mod(theta1_T_n, 2*pi);
%     theta2_T_n = mod(theta2_T_n, 2*pi);
% 
%     % Calcolo elementi orbitali
%     e_T_n = (norm(r2) - norm(r1))/(norm(r1) * r1_T_n(1)/norm(r1) - norm(r2) * r2_T_n(1)/norm(r2));
% 
%     if e_T_n <= 0 || e_T_n >= 1
%         continue;
%     end
% 
%     e_T_n_vec = [e_T_n_vec; e_T_n];
%     a_T_n = norm(r1) * (1 + e_T_n * cos(theta1_T_n))/(1 - e_T_n^2);
% 
%     % Calcolo velocità di trasferimento
%     [~,v1T_n] = par2car(a_T_n,e_T_n,i_T,Omega_T,omega_T_n,theta1_T_n,mu_t);
%     [~,v2T_n] = par2car(a_T_n,e_T_n,i_T,Omega_T,omega_T_n,theta2_T_n,mu_t);
% 
%     % Calcolo delta-v
%     delta_v_n = norm(v1T_n - vv_i) + norm(vv_f - v2T_n);
%     deltat = TOF(a_T_n,e_T_n,theta1_T_n,theta2_T_n,mu_t);
% 
%     % Aggiornamento soluzione ottima
%     if abs(delta_v_ott) > abs(delta_v_n) || delta_v_ott == 0
%         delta_v_ott = delta_v_n;
%         delta_v_ott1 = norm(v1T_n - vv_i);
%         delta_v_ott2 = norm(vv_f - v2T_n);
%         omega_T_ott = omega_T_n;
%         e_T_ott = e_T_n;
%         a_T_ott = a_T_n;
%         theta1_T_ott = theta1_T_n;
%         theta2_T_ott = theta2_T_n;
%         i_T_ott = i_T;
%         Omega_T_ott = Omega_T;
%         delta_t_ott = deltat;
%     end
% 
%     i = i+1;
%     delta_v_vect = [delta_v_vect; delta_v_n];
% end





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








function animateTransfer(a,e,i,OM,om,th_start,th_end,mu,satellite,v,doRotate,view_start,view_end,time_offset,h_text_time,h_text_vel,h_text_dv,delta_v)
    dth_anim = 0.01;
    
    % Calcolo anomalia eccentrica
    E_start = 2*atan(sqrt((1-e)/(1+e)) * tan(th_start/2));
    E_end = 2*atan(sqrt((1-e)/(1+e)) * tan(th_end/2));
    
    if E_end < E_start
        E_end = E_end + 2*pi;
    end
    
    % Calcolo tempo di volo
    n = sqrt(mu/a^3);
    DeltaT_exatto = (E_end - e*sin(E_end) - (E_start - e*sin(E_start))) / n;
    
    % Generazione punti lungo l'orbita
    th_vals = th_start:dth_anim:th_end;
    if th_end < th_start
        th_vals = [th_start:dth_anim:2*pi, 0:dth_anim:th_end];
    end
    N = length(th_vals);
    
    for k = 1:N
        th = th_vals(k);
        [rr, vv] = par2car(a,e,i,OM,om,mod(th,2*pi),mu);
        set(satellite, 'XData', rr(1), 'YData', rr(2), 'ZData', rr(3));
        
        % Aggiornamento informazioni
        v_norm = norm(vv);
        time_passed = time_offset + DeltaT_exatto * (k-1)/N;
        
        set(h_text_time, 'String', sprintf('Tempo: %.1f s', time_passed));
        set(h_text_vel, 'String', sprintf('Velocità: %.2f km/s', v_norm));
        if ~isempty(delta_v)
            set(h_text_dv, 'String', sprintf('ΔV: %.2f km/s', delta_v));
        end
        
        % Rotazione vista
        if doRotate
            az = view_start(1) + (view_end(1) - view_start(1)) * k / N;
            el = view_start(2) + (view_end(2) - view_start(2)) * k / N;
            view(az, el);
        end
        
        drawnow;
        frame = getframe(gcf);
        writeVideo(v, frame);
    end
    
    if doRotate
        view(view_end);
    end
end

%% ANIMAZIONE
fig = figure('Position', [100, 100, 1000, 600]);
hold on;

% Plot orbite
plotOrbit(a_i,e_i,i_i,OM_i,om_i,0,2*pi,dth,mu_t,'r');
plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,dth,mu_t,'b');
plotOrbit(a_T_ott,e_T_ott,i_T_ott,Omega_T_ott,omega_T_ott,theta1_T_ott,theta2_T_ott,dth,mu_t,'y');

% Punti iniziale e finale
[rr1,~] = par2car(a_T_ott,e_T_ott,i_T_ott,Omega_T_ott,omega_T_ott,theta1_T_ott,mu_t);
[rr2,~] = par2car(a_T_ott,e_T_ott,i_T_ott,Omega_T_ott,omega_T_ott,theta2_T_ott,mu_t);
plot3(rr1(1),rr1(2),rr1(3), '.', 'MarkerSize', 20, 'LineWidth', 2, 'Color', 'g');
plot3(rr2(1),rr2(2),rr2(3), '.', 'MarkerSize', 20, 'LineWidth', 2, 'Color', 'm');

% Terra
Terra3d(Rt);

% Legenda
xlabel('Asse X [Km]'); ylabel('Asse Y [Km]'); zlabel('Asse Z [Km]');
legend('Orbita iniziale','Orbita finale','Orbita di trasferimento',...
       'Punto iniziale','Punto finale','Terra', 'Location', 'northeastoutside');
grid on; view(-10,25);

% Textbox per informazioni
pos_time = [0.82 0.13 0.15 0.05];
pos_vel = [0.82 0.07 0.15 0.05];
pos_dv = [0.82 0.01 0.15 0.05];

h_text_time = annotation('textbox', pos_time, 'String', '', ...
    'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 12, 'FontWeight', 'bold');
h_text_vel = annotation('textbox', pos_vel, 'String', '', ...
    'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 12, 'FontWeight', 'bold');
h_text_dv = annotation('textbox', pos_dv, 'String', '', ...
    'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 12, 'FontWeight', 'bold');

% Satellite per animazione
satellite = plot3(NaN, NaN, NaN, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 8, 'HandleVisibility', 'off');

% Video writer
v = VideoWriter('Direct_Transfer.mp4', 'MPEG-4');
v.Quality = 100; v.FrameRate = 60; open(v);

%% ANIMAZIONE DELLA TRAIETTORIA
time_accum = 0;
view_start = [100, 30]; view_end = [0, 30];



% 1. Tratto orbita iniziale
animateTransfer(a_i,e_i,i_i,OM_i,om_i,th_i,theta1_T_ott,mu_t,satellite,v,true,...
               view_start,view_end,time_accum,h_text_time,h_text_vel,h_text_dv,[]);
time_accum = time_accum + TOF(a_i,e_i,th_i,theta1_T_ott,mu_t);

% 2. Trasferimento
animateTransfer(a_T_ott,e_T_ott,i_T_ott,Omega_T_ott,omega_T_ott,theta1_T_ott,theta2_T_ott,mu_t,...
               satellite,v,false,view_end,view_end,time_accum,h_text_time,h_text_vel,h_text_dv,delta_v_ott1);
time_accum = time_accum + delta_t_ott;

% 3. Tratto orbita finale
animateTransfer(a_f,e_f,i_f,OM_f,om_f,theta2_T_ott,th_f,mu_t,satellite,v,false,...
               view_end,view_end,time_accum,h_text_time,h_text_vel,h_text_dv,delta_v_ott);

%% MODIFICA FINALE PRIMA DI CHIUDERE IL VIDEO
set(h_text_dv, 'String', sprintf('ΔV: %.2f km/s', delta_v_ott));
set(h_text_time, 'String', sprintf('Tempo: %.1f s', time_accum));

% 1. Forza il rendering completo
drawnow expose; % Il più completo tra i drawnow

% 2. Cattura multipla con diverse tecniche
for k = 1:5
    % Tecnica 1: getframe standard
    frame = getframe(gcf);
    writeVideo(v, frame);
    
    % Calcola correttamente il rettangolo di cattura
    figPos = get(gcf, 'Position');
    
    pause(0.05); % Breve pausa tra le catture
end

% 3. Chiudi il video
close(v);