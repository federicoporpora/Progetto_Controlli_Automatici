%% MAIN PROJECT FILE - GRUPPO 35
clear all; close all; clc;

% Setup Cartelle per le figure
% Il percorso è relativo: da "code" saliamo a "root" (..) poi in "report/figs"
output_dir = fullfile('..', 'report', 'figs');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Creata cartella per le figure: %s\n', output_dir);
else
    fprintf('Cartella figure esistente: %s\n', output_dir);
end

%% DATI DEL PROGETTO (Traccia)
k1 = 1.63;
k2 = 0.73;
k3 = 1.55;
k4 = 2.5;

a1_eq = 0.0235;
a2_eq = 3.67;

%% PUNTO 1: EQUILIBRIO E LINEARIZZAZIONE

fprintf('\n--- PUNTO 1: Equilibrio ---\n');
check_eq2 = k2*sqrt(a1_eq) - k3*sqrt(a2_eq);
if abs(check_eq2) > 1e-4
    warning('Valori di equilibrio non coerenti.');
end

x1e = a1_eq;
x2e = a2_eq;
xe = [x1e; x2e];

% Calcolo ingresso di equilibrio
ue = (k1 * sqrt(x1e)) / k4;
fprintf('ue calcolato: %.5f\n', ue);

% Calcolo Jacobiane (Matrici A, B, C, D)
J_x = [ -k1/(2*sqrt(x1e)),      0;
         k2/(2*sqrt(x1e)), -k3/(2*sqrt(x2e)) ];

J_u = [ k4;
        0 ];

A = J_x;
B = J_u;
C = [0, 1];
D = 0;

%% PUNTO 2: FUNZIONE DI TRASFERIMENTO G(s)

fprintf('\n--- PUNTO 2: Funzione di Trasferimento ---\n');
sys_ss = ss(A, B, C, D);
poli = eig(A);
p1 = poli(1);
p2 = poli(2);

% Numeratore analitico per struttura a cascata: B(1)*A(2,1)
Num_val = B(1) * A(2,1);
s = tf('s');
G_s = Num_val / ((s - p1) * (s - p2));
mu_G = dcgain(G_s);

fprintf('G(s) calcolata.\n');
fprintf('Guadagno statico mu_G = %.2f (%.2f dB)\n', mu_G, 20*log10(mu_G));

% --- FIGURA 1: BODE G(s) ---
f1 = figure(1);
bode(G_s);
grid on;
title('Diagramma di Bode - Pianta G(s)');
saveas(f1, fullfile(output_dir, 'bode_G_plant.jpg'));
fprintf('Salvato: bode_G_plant.jpg\n');

%% PUNTO 3: SINTESI DEL REGOLATORE

fprintf('\n--- PUNTO 3: Sintesi Regolatore ---\n');

% 3.1 Regolatore Statico
% Specifica errore a regime <= 0.05 con W=3, D=2.5
% Calcolo da LaTeX: Kr >= 50. Scelto Kr = 60.
Kr = 60;
Rs_s = Kr;

% 3.2 Regolatore Dinamico (Rete Anticipatrice)
% Specifica Ta,5% <= 0.05s -> wc >= 60 rad/s.
wc_des = 60; 
% Specifica S% <= 20% -> Mf >= 45.6 -> Target Mf = 60 gradi.
Mf_target = 60; 

% Funzione parziale
L_parz = Kr * G_s;
[mag_wc, phase_wc] = bode(L_parz, wc_des);

% Calcolo parametri rete anticipatrice per inversione
% M* = 1 (0 dB) per avere attraversamento
% Phi* = Mf_target - (180 + fase_attuale)
M_star = 1;
phi_star_deg = Mf_target - (180 + phase_wc);
phi_star_rad = deg2rad(phi_star_deg);

% Formule di inversione
tau = (M_star - cos(phi_star_rad)) / (wc_des * sin(phi_star_rad));
alpha_tau = (cos(phi_star_rad) - 1/M_star) / (wc_des * sin(phi_star_rad));
alpha = alpha_tau / tau;

fprintf('Parametri Rete: tau = %.4f, alpha = %.4f\n', tau, alpha);

Rd_s = (1 + tau*s) / (1 + alpha*tau*s);

% Regolatore Totale
R_s = Rs_s * Rd_s;

% Funzione d'Anello
L_s = R_s * G_s;

%% GENERAZIONE GRAFICI DI PROGETTO (Punto 3)

% --- FIGURA 2: LOOP SHAPING (Zone proibite) ---
f2 = figure(2);
margin(L_s); 
grid on; 
hold on;

% Disegno Zone Proibite (patch grafiche)
% Specifica Disturbi: w <= 2, |L| >= 40dB
patch([0.1