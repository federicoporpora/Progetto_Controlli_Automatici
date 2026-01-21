close all; clear all; clc;

% Setup Cartella Output (come richiesto per il report)
output_dir = fullfile('..', 'report', 'figs');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% PARAMETRI DEL SISTEMA (Serbatoi)
k1 = 1.63;
k2 = 0.73;
k3 = 1.55;
k4 = 2.5;
% Valori di equilibrio (da traccia)
a1_eq = 0.0235;
a2_eq = 3.67;

%% PUNTO 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo dell'equilibrio e linearizzazione

x1e = a1_eq;
x2e = a2_eq;
xe = [x1e; x2e];

% Ingresso di equilibrio
ue = (k1 * sqrt(x1e)) / k4;

% Matrici del sistema linearizzato (Jacobiane)
A = [ -k1/(2*sqrt(x1e)),      0;
       k2/(2*sqrt(x1e)), -k3/(2*sqrt(x2e)) ];
B = [ k4; 0 ];
C = [ 0, 1 ];
D = 0;

fprintf('PUNTO 1: Equilibrio calcolato.\n');
fprintf('ue = %.5f\n', ue);

%% PUNTO 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo Funzione di Trasferimento G(s)
sys_ss = ss(A, B, C, D);
s = tf('s');

% Calcolo G(s) analitico (pulito)
poli = eig(A);
Num_val = B(1) * A(2,1); % Numeratore sistema a cascata
G = Num_val / ((s - poli(1)) * (s - poli(2)));

% Stampa Grafico Bode G(s)
figure(1); clf; set(gcf, 'Color', 'w');
margin(G);
title('Diagramma di Bode di G(s)');
grid on;
saveas(gcf, fullfile(output_dir, 'bode_G_plant.jpg'));

%% PUNTO 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SINTESI DEL REGOLATORE

% Specifiche di Progetto
e_star = 0.05;      % Errore a regime
W_star = 3;         % Ampiezza riferimento
D_star = 2.5;       % Ampiezza disturbo
Mf_target = 60;     % Margine di fase target (sovraelongazione < 20%)
wc_des = 60;        % Pulsazione critica (tempo assestamento < 0.05s)

% -- Regolatore Statico --
% Calcolo guadagno per errore a regime
mu_G = dcgain(G);
Kr_min = (W_star + D_star*mu_G)/e_star; % Formula approssimata lato sicurezza
Kr = 60; % Scelta progettuale (da calcoli precedenti era > 49.95)
R_s = Kr;

G_e = G * R_s; % Sistema esteso

% -- Regolatore Dinamico (Rete Anticipatrice) --
% Calcolo modulo e fase attuali a wc_des
val_Ge = evalfr(G_e, 1j*wc_des);
M_star = 1 / abs(val_Ge); 
phi_star_deg = Mf_target - (180 + rad2deg(angle(val_Ge)));

% Formule di inversione
phi_rad = deg2rad(phi_star_deg);
tau = (M_star - cos(phi_rad)) / (wc_des * sin(phi_rad));
alpha_tau = (cos(phi_rad) - 1/M_star) / (wc_des * sin(phi_rad));
alpha = alpha_tau / tau;

% Costruzione Regolatore
R_d = (1 + tau*s) / (1 + alpha*tau*s);
L = G_e * R_d; % Funzione d'anello

% -- Grafico Loop Shaping --
figure(2); clf; set(gcf, 'Color', 'w');
margin(L);
title('Loop Shaping L(s) con specifiche');
grid on; hold on;
xlim([1e-4 1e9]); ylim([-250 150]);

% Patch Zone Proibite (Stile semplice)
% Disturbo (Bassa frequenza)
patch([1e-4 2 2 1e-4], [40 40 150 150], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r');
text(0.01, 50, 'Disturbi', 'Color', 'r');
% Rumore (Alta frequenza)
patch([1e5 1e9 1e9 1e5], [-63 -63 100 100], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'b');
text(2e5, -40, 'Rumore', 'Color', 'b');

saveas(gcf, fullfile(output_dir, 'loop_shaping.jpg'));

% -- Grafico Bode Comparativo --
figure(3); clf; set(gcf, 'Color', 'w');
bodemag(G, 'b--'); hold on;
bodemag(L, 'r-');
grid on;
legend('G(s)', 'L(s)');
title('Confronto G(s) vs L(s)');
saveas(gcf, fullfile(output_dir, 'bode_comparativo.jpg'));

%% PUNTO 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST E SIMULAZIONI

% Funzioni di Sensitività
S = 1 / (1 + L); % Sensitività (Reiezione disturbi)
F = L / (1 + L); % Sensitività Complementare (Inseguimento / Rumore)

% -- Risposta al Gradino --
figure(4); clf; set(gcf, 'Color', 'w');
step(W_star * F, 0.15);
grid on; hold on;
yline(W_star*1.05, 'k--'); yline(W_star*0.95, 'k--');
title('Risposta al gradino w(t)=3');
saveas(gcf, fullfile(output_dir, 'step_response.jpg'));

% -- Funzioni di Sensitività --
figure(5); clf; set(gcf, 'Color', 'w');
bodemag(S, 'b'); hold on;
bodemag(F, 'r--');
grid on;
legend('|S|', '|F|');
title('Funzioni di Sensitività');
saveas(gcf, fullfile(output_dir, 'sensitivita.jpg'));

fprintf('Tutte le figure salvate in report/figs.\n');