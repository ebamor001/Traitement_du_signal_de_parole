%% Schéma complet et explicatif de la procédure Overlap–Add (avec fenêtres de Hamming)
clear; close all; clc;

L = 40;                   % longueur de trame
R = L/2;                  % recouvrement (50 %)
N = 3*L/2;                  % taille du signal de démonstration
n = 0:N-1;
x = sin(2*pi*n/40);       % signal test (sinus)
w = hamming(L)';          % fenêtre de Hamming

% === Étape 1 : Découpage et fenêtrage ===
x1 = x(1:L).*w;
x2 = x(R+1:R+L).*w;

% === Étape 2 : Reconstruction (somme et pondération) ===
x_recon = zeros(1,N);
weight  = zeros(1,N);

x_recon(1:L) = x_recon(1:L) + x1;
x_recon(R+1:R+L) = x_recon(R+1:R+L) + x2;
weight(1:L) = weight(1:L) + w;
weight(R+1:R+L) = weight(R+1:R+L) + w;

% === Étape 3 : Normalisation ===
x_recon_norm = x_recon ./ weight;

% === Figure unique et explicative ===
figure('Color','w','Position',[100 100 900 750]);

% ---- (a) Signal + fenêtres de Hamming + trames fenêtrées ----
% ---- (a) Signal + fenêtres de Hamming + trames fenêtrées ----
subplot(3,1,1);
hold on;

h1 = plot(n, x, 'k', 'LineWidth', 1.2);              % signal original
h2 = plot(0:L-1, w, 'r--', 'LineWidth', 1.2);        % fenêtre de Hamming 1
plot(R:R+L-1, w, 'r--', 'LineWidth', 1.2);           % fenêtre de Hamming 2
h3 = area(0:L-1, x1, 'FaceColor',[0.3 0.6 1], 'FaceAlpha',0.5, 'EdgeColor','none'); % trame 1 (bleue)
h4 = area(R:R+L-1, x2, 'FaceColor',[1 0.6 0.3], 'FaceAlpha',0.5, 'EdgeColor','none'); % trame 2 (orange)
plot(0:L-1, x1, 'b', 'LineWidth', 1.2);
plot(R:R+L-1, x2, 'Color',[1 0.4 0], 'LineWidth', 1.2);

title('(a) Découpage et fenêtrage du signal');
xlabel('Échantillons n'); ylabel('Amplitude');
legend([h1 h2 h3 h4], ...
       'Signal x[n]', ...
       'Fenêtres de Hamming (rouge pointillé)', ...
       'Trame 1 (bleue)', ...
       'Trame 2 (orange)');
grid on; axis tight;


% ---- (b) Somme des trames et somme des fenêtres ----
subplot(3,1,2);
yyaxis left;
plot(n, x_recon, 'm', 'LineWidth', 1.3);
ylabel('Somme des trames fenêtrées');
yyaxis right;
plot(n, weight, 'g--', 'LineWidth', 1.5);
ylabel('Somme des fenêtres de Hamming');
title('(b) Superposition des trames et pondération (somme des fenêtres)');
xlabel('Échantillons n');
legend('x_{somme}[n]', '\Sigma w[n]', 'Location','best');
grid on; axis tight;

% ---- (c) Reconstruction finale ----
subplot(3,1,3);
plot(n, x, 'k', 'LineWidth', 1.2); hold on;
plot(n, x_recon_norm, '--r', 'LineWidth', 1.2);
plot(n, x - x_recon_norm, 'b:');
legend('x[n] original', 'x_{recon}[n] (normalisé)', 'Erreur', 'Location','best');
title('(c) Signal reconstruit après division par la somme des fenêtres');
xlabel('Échantillons n'); ylabel('Amplitude');
grid on; axis tight;

sgtitle('Procédure complète d''addition–recouvrement (Overlap–Add) avec fenêtres de Hamming');
saveas(gcf, 'schema_overlap_add_complet_hamming.png');
