function [a_T_n,e_T_n,i_T,Omega_T,omega_T_n,theta1_T_n,theta2_T_n] = calcolo_orb(a_i,e_i,i_i,OM_i,om_i,a_f,e_f,i_f,OM_f,om_f,mu_s,x)

omega_T_n = x(3);
 [r1,vv_i] =  par2car(a_i,e_i,i_i,OM_i,om_i,x(1),mu_s);
    
        [r2,vv_f] =  par2car(a_f,e_f,i_f,OM_f,om_f,x(2),mu_s);
    
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
       
            % Definizione delle matrici di rotazione
            R_omega_T_n = [cos(x(3)), sin(x(3)), 0;
                          -sin(x(3)), cos(x(3)), 0;
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
            theta1_T_n=PositiveAngles(theta1_T_n);
            theta2_T_n=PositiveAngles(theta2_T_n);
        
            % Calcolo degli elementi orbitali
            e_T_n = (norm(r2) - norm(r1))/(norm(r1) * r1_T_n(1)/norm(r1) - norm(r2) * r2_T_n(1)/norm(r2)); 
            a_T_n = norm(r1) * (1 + e_T_n * cos(theta1_T_n))/(1 - e_T_n^2);
end

