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



%% ANIMAZIONE TRASFERIMENTO OTTIMALE SEMPLIFICATA
fig = figure('Position', [100, 100, 1250, 750]);
hold on;

% Plot solo delle parti di orbita rilevanti
plotOrbit(a_i,e_i,i_i,OM_i,om_i,th_i,theta1_T_ott,dth,mu_t,'r');
plotOrbit(a_f,e_f,i_f,OM_f,om_f,theta2_T_ott,th_f,dth,mu_t,'b');
plotOrbit(a_T_ott,e_T_ott,i_T_ott,Omega_T_ott,omega_T_ott,theta1_T_ott,theta2_T_ott,dth,mu_t,'y');

% Punti chiave
[rr_start,~] = par2car(a_i,e_i,i_i,OM_i,om_i,th_i,mu_t);
[rr_man1,~] = par2car(a_i,e_i,i_i,OM_i,om_i,theta1_T_ott,mu_t);
[rr_man2,~] = par2car(a_f,e_f,i_f,OM_f,om_f,theta2_T_ott,mu_t);
[rr_end,~] = par2car(a_f,e_f,i_f,OM_f,om_f,th_f,mu_t);

% Plot punti chiave
plot3(rr_start(1),rr_start(2),rr_start(3), '.', 'MarkerSize', 20, 'LineWidth', 2, 'Color', 'g');
plot3(rr_man1(1),rr_man1(2),rr_man1(3), 'ko', 'MarkerSize', 8, 'LineWidth', 1);
plot3(rr_man2(1),rr_man2(2),rr_man2(3), 'ko', 'MarkerSize', 8, 'LineWidth', 1);
plot3(rr_end(1),rr_end(2),rr_end(3), '.', 'MarkerSize', 20, 'LineWidth', 2, 'Color', 'm');

% Terra
Terra3d(Rt);

% Legenda semplificata
legend('Orbita iniziale','Orbita finale','Trasferimento','Punto iniziale','Punto finale',...
       'Location', 'northeastoutside');
xlabel('Asse X [Km]'); ylabel('Asse Y [Km]'); zlabel('Asse Z [Km]');
grid on; 
view(-10,25);

% Satellite per animazione
satellite = plot3(NaN, NaN, NaN, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 8, ...
                'HandleVisibility', 'off');

% Video writer
v = VideoWriter('Transfer_Animation.mp4', 'MPEG-4');
v.Quality = 100;
v.FrameRate = 60;
open(v);

%% ANIMAZIONE IN 3 FASI SEMPLIFICATA

% 1. Percorso orbita iniziale (da th_i a punto di manovra)
theta_range1 = linspace(th_i, theta1_T_ott, 60);
for theta = theta_range1
    [rr, ~] = par2car(a_i,e_i,i_i,OM_i,om_i,theta,mu_t);
    set(satellite, 'XData', rr(1), 'YData', rr(2), 'ZData', rr(3));
    drawnow;
    writeVideo(v, getframe(gcf));
end

% 2. Trasferimento (da punto di manovra iniziale a finale)
theta_range2 = linspace(theta1_T_ott, theta2_T_ott, 120);
for theta = theta_range2
    [rr, ~] = par2car(a_T_ott,e_T_ott,i_T_ott,Omega_T_ott,omega_T_ott,theta,mu_t);
    set(satellite, 'XData', rr(1), 'YData', rr(2), 'ZData', rr(3));
    drawnow;
    writeVideo(v, getframe(gcf));
end

% 3. Percorso orbita finale (da punto di manovra a th_f)
theta_range3 = linspace(theta2_T_ott, th_f, 60);
for theta = theta_range3
    [rr, ~] = par2car(a_f,e_f,i_f,OM_f,om_f,theta,mu_t);
    set(satellite, 'XData', rr(1), 'YData', rr(2), 'ZData', rr(3));
    drawnow;
    writeVideo(v, getframe(gcf));
end

% % Mantieni l'ultimo frame per 2 secondi
% for k = 1:120
%     writeVideo(v, getframe(gcf));
% end

close(v);
disp('Animazione completata. Video salvato come "Transfer_Animation.mp4"');