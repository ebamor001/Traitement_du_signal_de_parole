%% Schéma Overlap–Add avec fenêtres de Hamming (2 trames)
clear; close all; clc;

%% Paramètres
L = 40;              % longueur d'une trame
R = L/2;             % recouvrement 50 %
N = L + R;           % taille minimale pour 2 trames recouvrantes
n = 0:N-1;           % indices
x = sin(2*pi*n/40);  % signal test
w = hamming(L)';     % fenêtre de Hamming

%% === Étape 1 : Découpage et fenêtrage ===
x1 = x(1:L) .* w;            % trame 1 fenêtrée
x2 = x(R+1:R+L) .* w;        % trame 2 fenêtrée

%% === Étape 2 : Reconstruction par Overlap–Add ===
x_recon = zeros(1,N);
weight  = zeros(1,N);

x_recon(1:L) = x_recon(1:L) + x1;
x_recon(R+1:R+L) = x_recon(R+1:R+L) + x2;

weight(1:L) = weight(1:L) + w;
weight(R+1:R+L) = weight(R+1:R+L) + w;

%% === Étape 3 : Normalisation ===
x_recon_norm = x_recon ./ weight;

%% === Figure explicative ===
figure('Color','w','Position',[100 100 900 750]);


%% ----------------------------------------------------
%  (a) Signal + fenêtres + trames fenêtrées
%% ----------------------------------------------------
subplot(3,1,1); hold on;

% Signal entier
h1 = plot(n, x, 'k', 'LineWidth', 1.4);

% Fenêtres de Hamming (trame 1 & 2)
h2 = plot(0:L-1, w, 'r--', 'LineWidth', 1.2); 
plot(R:R+L-1, w, 'r--', 'LineWidth', 1.2);

% Aires colorées pour les trames
h3 = area(0:L-1, x1, 'FaceColor',[0.3 0.6 1], 'FaceAlpha',0.45, 'EdgeColor','none');
h4 = area(R:R+L-1, x2, 'FaceColor',[1.0 0.6 0.3], 'FaceAlpha',0.45, 'EdgeColor','none');

% Contours des trames
plot(0:L-1, x1, 'Color',[0.2 0.45 1], 'LineWidth', 1.3);
plot(R:R+L-1, x2, 'Color',[1 0.4 0], 'LineWidth', 1.3);

title('(a) Découpage et fenêtrage du signal');
xlabel('Échantillons n');
ylabel('Amplitude');

legend([h1 h2 h3 h4], ...
       'Signal x[n]', ...
       'Fenêtres de Hamming', ...
       'Trame 1 (bleue)', ...
       'Trame 2 (orange)', ...
       'Location','best');

grid on; axis tight;


%% ----------------------------------------------------
%  (b) Somme des trames et somme des fenêtres
%% ----------------------------------------------------

subplot(3,1,2); hold on;

%% Axe gauche : somme des trames
yyaxis left
plot(n, x_recon, 'm', 'LineWidth', 1.3);
ylabel('Somme des trames');
ylim([-1.2 1.2]);
yticks([-1 0 1]);   

%% Axe droit : somme des fenêtres
yyaxis right
plot(n, weight, 'g--', 'LineWidth', 1.5);
ylabel('\Sigma w[n]');
ylim([-1.2 1.2]);
yticks([-1 0 1 2]);   


title('(b) Superposition des trames et somme des fenêtres');
xlabel('Échantillons n');

legend('x_{somme}[n]', '\Sigma w[n]', 'Location','best');
grid on; axis tight;


%% ----------------------------------------------------
%  (c) Reconstruction finale
%% ----------------------------------------------------
subplot(3,1,3); hold on;

plot(n, x, 'k', 'LineWidth', 1.4);
plot(n, x_recon_norm, '--r', 'LineWidth', 1.4);
plot(n, x - x_recon_norm, 'b:', 'LineWidth', 1.2);

legend('x[n] original', 'x_{recon}[n] (normalisé)', 'Erreur', 'Location','best');

title('(c) Reconstruction finale après normalisation');
xlabel('Échantillons n');
ylabel('Amplitude');

grid on; axis tight;

%% Titre global
sgtitle('Procédure complète Overlap–Add (Fenêtre de Hamming)', ...
        'FontSize', 14, 'FontWeight', 'bold');

%% Sauvegarde
saveas(gcf, 'schema_overlap_add_complet_hamming.png');
