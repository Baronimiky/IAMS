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

Rt = 6378; % raggio Terra
Rm = 6.955 *1e+5; % raggio Sole


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
            if e_T_n <=0 || e_T_n >= 1
                continue;
            end
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
%exportgraphics(gcf, 'delta_v_for.pdf', 'ContentType', 'vector', 'BackgroundColor','none')
    
    
%% Plot orbite
figure
plotOrbit(a_i,e_i,i_i,OM_i,om_i,0,2*pi,0.01,5,mu_s);
plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,0.01,3,mu_s);
sole(Rm*20);
legend('orbita Terra','orbita Asteroide','Sole');
xlabel('Asse X [Km]')
ylabel('Asse Y [Km]')
zlabel('Asse Z [Km]')
box on
grid on
view(-20, 15)
%exportgraphics(gcf, 'delta_v_for.pdf', 'ContentType', 'vector', 'BackgroundColor','none')

%% Plot Orbite con trasferimento Grid Search
plotOrbit(a_i,e_i,i_i,OM_i,om_i,0,2*pi,0.01,5,mu_s);
plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,0.01,3,mu_s);
plotOrbit( a_T_ott,e_T_ott,i_T_ott,Omega_T_ott,omega_T_ott,theta1_T_ott,theta2_T_ott,0.0001,3,mu_s);
[rr1,~] = par2car(a_T_ott,e_T_ott,i_T_ott,Omega_T_ott,omega_T_ott,theta1_T_ott,mu_s);
plot3(rr1(1),rr1(2),rr1(3), '.', 'MarkerSize', 20, 'LineWidth', 2);
[rr2] = par2car(a_T_ott,e_T_ott,i_T_ott,Omega_T_ott,omega_T_ott,theta2_T_ott,mu_s);
plot3(rr2(1),rr2(2),rr2(3), '.', 'MarkerSize', 20, 'LineWidth', 2);
sole(Rm*20);
Terra3d(Rt*700, rr1);
asteroid(Rt*500,rr2);
legend('orbita terra','orbita asteroide','trasferimento','punto iniziale (Terra)','punto finale (Asteroide)','sole');

%% Ottimizzazione Velocità con Fmincon
clc
% x = [th_i;th_f;omega_T_n];

fun_ott_vel = @(x) fun_vel(a_i,e_i,i_i,OM_i,om_i,x(1),a_f,e_f,i_f,OM_f,om_f,x(2),mu_s,x(3));
x0 = [pi;pi;pi];
A = [];
B = [];
Aeq = [];
Beq = [];
Lb = [0;0;0];
Ub = [2*pi;2*pi;2*pi];

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', ...
    'MaxFunctionEvaluations', 100000, 'MaxIterations', 100000, ...
    'PlotFcn', {@optimplotfval, @optimplotfirstorderopt, @optimplotconstrviolation});

nonlcon = @(x) vinc_nlin(x,a_i,e_i,i_i,OM_i,om_i,a_f,e_f,i_f,OM_f,om_f,mu_s);

x_ott = fmincon(fun_ott_vel, x0, A, B, Aeq, Beq, Lb, Ub, nonlcon, options);
%exportgraphics(gcf, 'fmincon.pdf', 'ContentType', 'vector', 'BackgroundColor','none')
delta_v_ott_fun = fun_ott_vel(x_ott);

%% Ottimizzazione Velocità con Ga
clc

    options = optimoptions('ga', ...
    'UseParallel', false, ...
    'MaxGenerations', 2000, ...
    'PopulationSize', 2000, ...
    'FunctionTolerance', 1e-10, ...
    'PlotFcn', {@gaplotbestf}, ...
    'Display', 'iter');
    
    x_ga = ga(fun_ott_vel,3,A,B,Aeq,Beq,Lb,Ub, nonlcon,options);
    delta_v_ott_fun_ga = fun_ott_vel(x_ga);
    %exportgraphics(gcf, 'ga.pdf', 'ContentType', 'vector', 'BackgroundColor','none')
    [a_T_ga_vec,e_T_ga_vec,i_T_ga_vec,Omega_T_ga_vec,omega_T_ga_vec,theta1_T_ga_vec,theta2_T_ga_vec] = calcolo_orb(a_i,e_i,i_i,OM_i,om_i,a_f,e_f,i_f,OM_f,om_f,mu_s,x_ga);
%% Plot velocità ottimale di Ga

figure('Position', [100, 100, 1000, 600])  
plotOrbit(a_i,e_i,i_i,OM_i,om_i,0,2*pi,0.01,5,mu_s);
plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,0.01,3,mu_s);
plotOrbit(a_T_ga_vec,e_T_ga_vec,i_T_ga_vec,Omega_T_ga_vec,omega_T_ga_vec,theta1_T_ga_vec,theta2_T_ga_vec,0.0001,3,mu_s);
[rr1,~] = par2car(a_T_ga_vec,e_T_ga_vec,i_T_ga_vec,Omega_T_ga_vec,omega_T_ga_vec,theta1_T_ga_vec,mu_s);
plot3(rr1(1),rr1(2),rr1(3), '.', 'MarkerSize', 20, 'LineWidth', 2);
[rr2,~] = par2car(a_T_ga_vec,e_T_ga_vec,i_T_ga_vec,Omega_T_ga_vec,omega_T_ga_vec,theta2_T_ga_vec,mu_s);
plot3(rr2(1),rr2(2),rr2(3), '.', 'MarkerSize', 20, 'LineWidth', 2);
sole(Rm*20);
Terra3d(Rt*700, rr1);
asteroid(Rt*500,rr2);
grid on
xlabel('Asse X [Km]')
ylabel('Asse Y [Km]')
zlabel('Asse Z [Km]')
view(-10,25)
legend('dati', 'Location', 'northeast')
legend('orbita terra','orbita asteroide','trasferimento','punto iniziale (Terra)','punto finale (Asteroide)','sole');
%exportgraphics(gcf, 'orbita_ga.pdf', 'ContentType', 'vector', 'BackgroundColor','none')
%% Ottimizzazione Velocità con Surrogateopt

opts = optimoptions('surrogateopt', ...
    'UseParallel', false, ...         % usa CPU multiple se puoi
    'MaxTime', 600, ...              % tempo massimo in secondi
    'PlotFcn', 'surrogateoptplot','MaxFunctionEvaluations',2000);  % opzionale: mostra progresso

[x_surr, fval_surr] = surrogateopt(fun_ott_vel, Lb, Ub, opts);

%% Ottimizzazione Velocità con Particleswarm

opts = optimoptions('particleswarm', ...
    'SwarmSize', 100, ...
    'UseParallel', false, ...
    'MaxIterations', 200, ...
    'Display', 'iter');

[x_par, fval_par] = particleswarm(fun_ott_vel, 3, Lb,Ub, opts);

%% Ottimizzazione Velocità con Fmincon annidati
clc

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', ...
    'MaxFunctionEvaluations', 100000, 'MaxIterations', 100000, ...
    'PlotFcn', {@optimplotfval, @optimplotfirstorderopt, @optimplotconstrviolation});


problem = createOptimProblem('fmincon', 'objective', fun_ott_vel, 'x0', x0, 'lb', Lb, 'ub', Ub, 'options', options);
ms = MultiStart('UseParallel', true);
[x, fval_multi] = run(ms, problem, 50);

%% Ottimizzazione del tempo di manovra con Fmincon

% x = [th_i;th_f;omega_T_n];

fun_ott_tem = @(x) fun_tem(a_i,e_i,i_i,OM_i,om_i,x(1),a_f,e_f,i_f,OM_f,om_f,x(2),mu_s,x(3));
x0 = [rand*2*pi;rand*2*pi;rand*2*pi];
A = [];
B = [];
Aeq = [];
Beq = [];
Lb = [0;0;0];
Ub = [2*pi;2*pi;2*pi];

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', ...
    'MaxFunctionEvaluations', 100000, 'MaxIterations', 100000, ...
    'PlotFcn', {@optimplotfval, @optimplotfirstorderopt, @optimplotconstrviolation});

nonlcon = @(x) vinc_nlin(x,a_i,e_i,i_i,OM_i,om_i,a_f,e_f,i_f,OM_f,om_f,mu_s);

x_ott_tem = fmincon(fun_ott_tem,x0,A,B,Aeq,Beq,Lb,Ub,nonlcon,options);

delta_t_ott_fun = fun_ott_tem(x_ott_tem);

[a_T_ott_tem,e_T_ott_tem,i_T_ott_tem,Omega_T_ott_tem,omega_T_ott_tem,theta1_T_ott_tem,theta2_T_ott_tem] = calcolo_orb(a_i,e_i,i_i,OM_i,om_i,a_f,e_f,i_f,OM_f,om_f,mu_s,x_ott_tem);

%% Ottimizzazione del tempo di manovra con Ga

    options = optimoptions('ga', ...
    'UseParallel', true, ...
    'MaxGenerations', 2000, ...
    'PopulationSize', 2000, ...
    'FunctionTolerance', 1e-6, ...
    'PlotFcn', {@gaplotbestf}, ...
    'Display', 'iter');
    
    x_ga_tem = ga(fun_ott_tem,3,A,B,Aeq,Beq,Lb,Ub,nonlcon,options);
    delta_t_ott_fun_ga = fun_ott_tem(x_ga_tem);

    [a_T_ga_tem,e_T_ga_tem,i_T_ga_tem,Omega_T_ga_tem,omega_T_ga_tem,theta1_T_ga_tem,theta2_T_ga_tem] = calcolo_orb(a_i,e_i,i_i,OM_i,om_i,a_f,e_f,i_f,OM_f,om_f,mu_s,x_ga_tem);
    
%% Plot ottimizzando il tempo con Ga

figure('Position', [100, 100, 1000, 600])
plotOrbit(a_i,e_i,i_i,OM_i,om_i,0,2*pi,0.01,5,mu_s);
plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,0.01,3,mu_s);
plotOrbit(a_T_ga_tem,e_T_ga_tem,i_T_ga_tem,Omega_T_ga_tem,omega_T_ga_tem,theta1_T_ga_tem,theta2_T_ga_tem,0.0001,3,mu_s);
[rr1,~] = par2car(a_T_ga_tem,e_T_ga_tem,i_T_ga_tem,Omega_T_ga_tem,omega_T_ga_tem,theta1_T_ga_tem,mu_s);
plot3(rr1(1),rr1(2),rr1(3), '.', 'MarkerSize', 20, 'LineWidth', 2);
[rr2] = par2car(a_T_ga_tem,e_T_ga_tem,i_T_ga_tem,Omega_T_ga_tem,omega_T_ga_tem,theta2_T_ga_tem,mu_s);
plot3(rr2(1),rr2(2),rr2(3), '.', 'MarkerSize', 20, 'LineWidth', 2);
sole(Rm*20);
Terra3d(Rt*700, rr1);
asteroid(Rt*500,rr2);
grid on
xlabel('Asse X [Km]')
ylabel('Asse Y [Km]')
zlabel('Asse Z [Km]')
view(-10,25)
legend('dati', 'Location', 'northeast')
legend('orbita terra','orbita asteroide','trasferimento','punto iniziale (Terra)','punto finale (Asteroide)','sole');
exportgraphics(gcf, 'orbita_gatem.pdf', 'ContentType', 'vector', 'BackgroundColor','none')
%% ottimizzazione sia del tempo sia della velocità con Ga/Fmincon

clc
% x = [th_i;th_f;omega_T_n];
delta_v_ott_O = 6.17;
dt_ott_O = 5*1e5;
alpha = 0.4921;
beta = 0.5079;

fun_ott_vec_tem = @(x) fun_vec_tem(a_i,e_i,i_i,OM_i,om_i,x(1),a_f,e_f,i_f,OM_f,om_f,x(2),mu_s,x(3),delta_v_ott_O,dt_ott_O,alpha,beta);
x0 = [rand*2*pi;rand*2*pi;rand*2*pi];
A = [];
B = [];
Aeq = [];
Beq = [];
Lb = [0;0;0];
Ub = [2*pi;2*pi;2*pi];

%options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', ...
%    'MaxFunctionEvaluations', 100000, 'MaxIterations', 100000, ...
%    'PlotFcn', {@optimplotfval, @optimplotfirstorderopt, @optimplotconstrviolation});

options = optimoptions('ga', ...
    'UseParallel', true, ...
    'MaxGenerations', 2000, ...
    'PopulationSize', 2000, ...
    'FunctionTolerance', 1e-6, ...
    'PlotFcn', {@gaplotbestf}, ...
    'Display', 'iter');

nonlcon = @(x) vinc_nlin(x,a_i,e_i,i_i,OM_i,om_i,a_f,e_f,i_f,OM_f,om_f,mu_s);

%x_ott_vec_tem = fmincon(fun_ott_vec_tem,x0,A,B,Aeq,Beq,Lb,Ub,nonlcon,options);
x_ott_vec_tem = ga(fun_ott_vec_tem,3,A,B,Aeq,Beq,Lb,Ub,nonlcon,options);

O_ott = fun_ott_vec_tem(x_ott_vec_tem);

delta_v_ott_vec_tem = fun_vel(a_i,e_i,i_i,OM_i,om_i,x_ott_vec_tem(1),a_f,e_f,i_f,OM_f,om_f,x_ott_vec_tem(2),mu_s,x_ott_vec_tem(3));
dt_ott_vec_tem = fun_tem(a_i,e_i,i_i,OM_i,om_i,x_ott_vec_tem(1),a_f,e_f,i_f,OM_f,om_f,x_ott_vec_tem(2),mu_s,x_ott_vec_tem(3));


[a_T_ott_vec_tem,e_T_ott_vec_tem,i_T_ott_vec_tem,Omega_T_ott_vec_tem,omega_T_ott_vec_tem,theta1_T_ott_vec_tem,theta2_T_ott_vec_tem] = calcolo_orb(a_i,e_i,i_i,OM_i,om_i,a_f,e_f,i_f,OM_f,om_f,mu_s,x_ott_vec_tem);

 %% plot ottimizzando il tempo e la velocità

figure('Position', [100, 100, 1000, 600])
plotOrbit(a_i,e_i,i_i,OM_i,om_i,0,2*pi,0.01,5,mu_s);
plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,0.01,3,mu_s);
plotOrbit(a_T_ott_vec_tem,e_T_ott_vec_tem,i_T_ott_vec_tem,Omega_T_ott_vec_tem,omega_T_ott_vec_tem,theta1_T_ott_vec_tem,theta2_T_ott_vec_tem,0.0001,3,mu_s);
[rr1,~] = par2car(a_T_ott_vec_tem,e_T_ott_vec_tem,i_T_ott_vec_tem,Omega_T_ott_vec_tem,omega_T_ott_vec_tem,theta1_T_ott_vec_tem,mu_s);
plot3(rr1(1),rr1(2),rr1(3), '.', 'MarkerSize', 20, 'LineWidth', 2);
[rr2] = par2car(a_T_ott_vec_tem,e_T_ott_vec_tem,i_T_ott_vec_tem,Omega_T_ott_vec_tem,omega_T_ott_vec_tem,theta2_T_ott_vec_tem,mu_s);
plot3(rr2(1),rr2(2),rr2(3), '.', 'MarkerSize', 20, 'LineWidth', 2);
sole(Rm*20);
Terra3d(Rt*700, rr1);
asteroid(Rt*500,rr2);
grid on
xlabel('Asse X [Km]')
ylabel('Asse Y [Km]')
zlabel('Asse Z [Km]')
view(-10,25)
legend('dati', 'Location', 'northeast')
legend('orbita terra','orbita asteroide','trasferimento','punto iniziale (Terra)','punto finale (Asteroide)','sole');
%exportgraphics(gcf, 'orbita_gavectem.pdf', 'ContentType', 'vector', 'BackgroundColor','none')

%% scelta parametri alpha e beta

% Dati di partenza (sostituisci con i tuoi valori)
alpha_values = linspace(0.1, 10, 10);  % Variazione dei valori di alpha
beta_values = linspace(0.1, 10, 10);   % Variazione dei valori di beta

% Variabili per memorizzare i risultati
delta_v_results = [];  % Risultati delta-v
dt_results = [];       % Risultati dt

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', ...
   'MaxFunctionEvaluations', 100000, 'MaxIterations', 100000, ...
    'PlotFcn', {@optimplotfval, @optimplotfirstorderopt, @optimplotconstrviolation});

% Inizializzazione del ciclo per trovare la curva di Pareto
for alpha = alpha_values
    for beta = beta_values
        % Funzione obiettivo con alpha e beta
        fun_ott_vec_tem = @(x) fun_vec_tem(a_i, e_i, i_i, OM_i, om_i, x(1), a_f, e_f, i_f, OM_f, om_f, x(2), mu_s, x(3), delta_v_ott_O, dt_ott_O, alpha, beta);

        % Esegui l'ottimizzazione con fmincon
        [x_ott_vec_tem, ~] = fmincon(@(x) fun_ott_vec_tem(x), x0, A, B, Aeq, Beq, Lb, Ub, nonlcon, options);

        % Calcola delta_v e dt per la soluzione ottimizzata
        delta_v_ott_vec_tem = fun_vel(a_i, e_i, i_i, OM_i, om_i, x_ott_vec_tem(1), a_f, e_f, i_f, OM_f, om_f, x_ott_vec_tem(2), mu_s, x_ott_vec_tem(3));
        dt_ott_vec_tem = fun_tem(a_i, e_i, i_i, OM_i, om_i, x_ott_vec_tem(1), a_f, e_f, i_f, OM_f, om_f, x_ott_vec_tem(2), mu_s, x_ott_vec_tem(3));

        % Salva i risultati
        delta_v_results = [delta_v_results; delta_v_ott_vec_tem];
        dt_results = [dt_results; dt_ott_vec_tem];
    end
end

%% plot curva di pareto

% Identificazione delle soluzioni non dominate
pareto_front = true(length(delta_v_results), 1);
for i = 1:length(delta_v_results)
    for j = 1:length(delta_v_results)
        if (delta_v_results(j) <= delta_v_results(i) && dt_results(j) < dt_results(i)) || ...
           (delta_v_results(j) < delta_v_results(i) && dt_results(j) <= dt_results(i))
            pareto_front(i) = false;
            break; % Se i risultati i sono dominati da j, esci dal ciclo
        end
    end
end

% Visualizzazione della curva di Pareto
figure;
scatter(delta_v_results, dt_results, 'filled');
xlabel('Delta-v (m/s)');
ylabel('Tempo di volo (s)');
title('Curva di Pareto - Ottimizzazione Delta-v e Tempo di Volo');
grid on;
% Traccia la curva di Pareto
hold on;
plot(delta_v_results(pareto_front), dt_results(pareto_front), 'r-', 'LineWidth', 2); % Traccia la frontiera di Pareto
legend('Soluzioni', 'Frontiera di Pareto');

%% scelta alpha e beta

% Trova il punto più vicino all'origine (come compromesso)
norm_vals = sqrt((delta_v_results/delta_v_ott_O).^2 + (dt_results/dt_ott_O).^2);
[~, best_idx] = min(norm_vals);
hold on;
plot(delta_v_results(best_idx), dt_results(best_idx), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Delta-v (m/s)');
ylabel('Tempo di volo (s)');
title('Curva di Pareto - Ottimizzazione Delta-v e Tempo di Volo');
scatter(delta_v_results, dt_results, 'filled');
legend('Soluzioni', 'Compromesso scelto');
exportgraphics(gcf, 'pareto.pdf', 'ContentType', 'vector', 'BackgroundColor','none')
% --- CALCOLO α E β OTTIMI ---
delta_v_star = delta_v_results(best_idx);
dt_star = dt_results(best_idx);

% Normalizzazione
delta_v_star_norm = delta_v_star / delta_v_ott_O;
dt_star_norm = dt_star / dt_ott_O;

alpha = 1 / delta_v_star_norm;
beta = 1 / dt_star_norm;

% Normalizzazione per somma 1
alpha_beta_sum = alpha + beta;
alpha = alpha / alpha_beta_sum;
beta = beta / alpha_beta_sum;

fprintf('Alpha ottimo: %.4f\n', alpha);
fprintf('Beta ottimo:  %.4f\n', beta);

alpha =  0.4921;
beta = 0.5079;
