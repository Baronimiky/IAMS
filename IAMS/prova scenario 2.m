%%
clear all
clc

a_i = 1.4946e+08;
e_i = 0.016;
i_i = 9.1920e-05;
OM_i = 2.7847;
om_i = 5.2643;

% 2.156790 0.083223 0.23 281.03 91.69 gruppo 2 Giulio
% 2.205422 0.162780 2.50 294.26 1.20 gruppo 5 Matteo
% 1.066293 0.315212 8.51 205.80 172.04 noi

a_f = 1.066293 *1.496*1e8;
e_f = 0.315212;
i_f = 8.51*pi/180;
OM_f = 205.80*pi/180;
om_f = 172.04*pi/180;


mu_s = 132712440018; % costante gravitazionale sole [km^3/s^2]

delta_v_ott = 0;
delta_v_vect = [];
e_T_n_vec = [];
i = 0;

%%
for th_i = 0:0.1:2*pi
    for th_f = 0:0.1:2*pi

        [r1,vv_i] = par2car(a_i,e_i,i_i,OM_i,om_i,th_i,mu_s);
    
        [r2,vv_f] = par2car(a_f,e_f,i_f,OM_f,om_f,th_f,mu_s);
    
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
            e_T_n_vec = [e_T_n_vec; e_T_n];
            %e_T_n = norm(e_T_n);
            a_T_n = norm(r1) * (1 + e_T_n * cos(theta1_T_n))/(1 - e_T_n^2);
            
            % Determinazione delle velocità nei punti di manovra
            
            %p_T_n = a_T_n * (1-e_T_n^2);
            %v1T_n = sqrt(mu_s./p_T_n) * [-sin(theta1_T_n); e_T_n+cos(theta1_T_n); 0];
            %v2T_n = sqrt(mu_s./p_T_n) * [-sin(theta2_T_n); e_T_n+cos(theta2_T_n); 0];
            [~,v1T_n] = par2car(a_T_n,e_T_n,i_T,Omega_T,omega_T_n,theta1_T_n,mu_s);
            [~,v2T_n] = par2car(a_T_n,e_T_n,i_T,Omega_T,omega_T_n,theta2_T_n,mu_s);
            
            % Calcolo del delta-v totale
            delta_v_n = norm(v1T_n - vv_i) + norm(vv_f - v2T_n);
            
            if abs(delta_v_ott) > abs(delta_v_n) || delta_v_ott == 0
                delta_v_ott = delta_v_n;
                th_i_ott = th_i;
                th_f_ott = th_f;
                omega_T_ott = omega_T_n;
                e_T_ott = e_T_n;
                a_T_ott = a_T_n;
                theta1_T_ott = theta1_T_n;
                theta2_T_ott = theta2_T_n;
                i_T_ott =i_T;
                Omega_T_ott = Omega_T;
            end
          i = i+1
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
    dt=TOF(a_T_ott,e_T_ott, theta1_T_ott,theta2_T_ott,mu_s)
  j = 1:1:i;
   plot(j,delta_v_vect);

    
    % [v1T_n] = orbVelocity(a_T_n, e_T_n, i_T, Omega_T, omega_T_n, theta1_T_n, mu_s);
    % [v2T_n] = orbVelocity(a_T_n, e_T_n, i_T, Omega_T, omega_T_n, theta2_T_n, mu_s);
    
    
%%
figure
plotOrbit(a_i,e_i,i_i,OM_i,om_i,0,2*pi,0.01,5,mu_s);
plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,0.01,3,mu_s);
plotOrbit( a_T_ott,e_T_ott,i_T_ott,Omega_T_ott,omega_T_ott,theta1_T_ott,theta2_T_ott,0.0001,3,mu_s);
legend('inizio','fine','trasferimento');
