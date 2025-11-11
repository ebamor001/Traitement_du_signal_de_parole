%% Préliminaire 3 : Filtrage du signal
clear; close all; clc;

load('fcno01fz.mat');   % Chargement du signal
signal = fcno01fz;
Fe = 8000;              % Fréquence d'échantillonnage
k0_values = [1 10 100 300];

%% Figure globale pour les spectrogrammes 
figure('Name', 'Spectrogrammes pour différentes valeurs de k0', 'NumberTitle', 'off');

for i = 1:length(k0_values)
    k0 = k0_values(i);

    % --- Définition du filtre ---
    h = zeros(1, k0+1);
    h(1) = 1; 
    h(end) = 1;  % h(k) = δ(k) + δ(k - k0)

    % --- Application du filtre ---
    y = filter(h, 1, signal);

    % Lecture audio
    disp(['Lecture du signal filtré pour k0 = ', num2str(k0)]);
    sound(y / max(abs(y)), Fe);  % normalisation + lecture du signal filtré
    pause(length(y)/Fe + 1);     % attendre la fin du son avant le suivant


    % --- Affichage du spectrogramme ---
    subplot(2,2,i);
    spectrogram(y, hamming(256), 128, 512, Fe, 'yaxis');
    title(['Spectrogramme - Signal filtré (k_0 = ', num2str(k0), ')']);
    ylim([0 4]); % 0–4 kHz utile
end

sgtitle('Comparaison des spectrogrammes pour différentes valeurs de k_0');

%%  Boucle séparée pour les figures d'analyse du filtre 
for k0 = k0_values
    h = zeros(1, k0+1);
    h(1) = 1; h(end) = 1;

    y = filter(h, 1, signal);
    [H, f] = freqz(h, 1, 512, Fe);

    figure('Name', ['Analyse du filtre pour k0 = ', num2str(k0)], 'NumberTitle', 'off');

    subplot(2,2,1);
    stem(0:k0, h, 'filled');
    title(['Réponse impulsionnelle h(k), k0 = ', num2str(k0)]);
    xlabel('k'); grid on;

    subplot(2,2,2);
    plot(f, abs(H), 'LineWidth', 1.2);
    title('Module de la réponse en fréquence');
    xlabel('Fréquence (Hz)'); ylabel('|H(f)|'); grid on;

    subplot(2,2,3);
    plot(signal); 
    title('Signal original');
    xlabel('Échantillons'); 
    xlim([1 length(signal)]); 
    A = max(abs(signal));
    ylim([-A A]);grid on;


    subplot(2,2,4);
    plot(y); 
    title(['Signal filtré (k0 = ', num2str(k0), ')']);
    xlabel('Échantillons'); 
    xlim([1 length(y)]); 
    A = max(abs(y));
    ylim([-A A]);grid on;
end
