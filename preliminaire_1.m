%% Préliminaire 1 : Analyse d'un bruit blanc gaussien
clear; close all; clc;

N = 2048;          % longueur du signal
sigma = 1;         % écart-type
mu = 0;            % moyenne

%Génération d'un bruit blanc gaussien 
x = sigma * randn(1, N) + mu;

% Autocorrélation théorique
R_theo = sigma^2 * [zeros(1, N-1), 1, zeros(1, N-1)]; % même longueur que xcorr

% Autocorrélation estimée
[R_est_biased, lags] = xcorr(x, "biased");
[R_est_unbiased, ~]  = xcorr(x, "unbiased");

% Affichage
figure;
    
subplot(3,1,1); % Théorique
stem(lags, R_theo, 'r', 'filled');
xlabel('Décalage k');
xlim([min(lags) max(lags)]);
ylabel('R_x(k)');
title('Autocorrélation théorique');
grid on;

subplot(3,1,2); % Estimée biased
plot(lags, R_est_biased, 'b');
xlabel('Décalage k');
xlim([min(lags) max(lags)]);
ylabel('R_x(k)');
title('Autocorrélation estimée (biased)');
grid on;

subplot(3,1,3); % Estimée unbiased
plot(lags, R_est_unbiased, 'g');
xlabel('Décalage k');
xlim([min(lags) max(lags)]);
ylabel('R_x(k)');
title('Autocorrélation estimée (unbiased)');
grid on;

sgtitle('Comparaison des autocorrélations (Théorique / Biased / Unbiased)');

%% Densité spectrale de puissance (DSP) estimée
Fe = 1; Nfft = 1024;
[Pp, f]  = Mon_Periodogramme(x, Nfft, Fe);
[Pd, ~]  = Mon_Daniell(x, Nfft, Fe, 9); %9= taille de lissage
[Pb, ~]  = Mon_Bartlett(x, Nfft, Fe);
[Pw, ~]  = Mon_Welch(x, Nfft, Fe, 0.5);
[Pc, ~]  = Mon_Correlogramme(x, Nfft, Fe);


% DSP théorique (constante = sigma^2)
P_theo = sigma^2 * ones(size(f));
% Conversion en dB
Pp_dB = 10*log10(Pp);
Pb_dB = 10*log10(Pb);
Pw_dB = 10*log10(Pw);
Pd_dB = 10*log10(Pd);
Pc_dB = 10*log10(Pc);
Pth_dB = 10*log10(P_theo);

% Trouver les bornes globales 
ymin = min([Pp_dB Pb_dB Pw_dB Pd_dB Pc_dB]) - 2;
ymax = max([Pp_dB Pb_dB Pw_dB Pd_dB Pc_dB]) + 2;

figure;

subplot(5,1,1);
plot(f, Pp_dB, 'b', 'LineWidth', 1.2); hold on;
plot(f, Pth_dB, 'k--', 'LineWidth', 1);
xlabel('Fréquence (Hz)'); ylabel('Puissance (dB)');
title('Périodogramme simple');
grid on; ylim([ymin ymax]);

subplot(5,1,2);
plot(f, Pd_dB, 'c', 'LineWidth', 1.2); hold on;
plot(f, Pth_dB, 'k--', 'LineWidth', 1);
xlabel('Fréquence (Hz)'); ylabel('Puissance (dB)');
title('Méthode de Daniell (lissage)');
grid on; ylim([ymin ymax]);


subplot(5,1,3);
plot(f, Pb_dB, 'm', 'LineWidth', 1.2); hold on;
plot(f, Pth_dB, 'k--', 'LineWidth', 1);
xlabel('Fréquence (Hz)'); ylabel('Puissance (dB)');
title('Méthode de Bartlett');
grid on; ylim([ymin ymax]);

subplot(5,1,4);
plot(f, Pw_dB, 'g', 'LineWidth', 1.2); hold on;
plot(f, Pth_dB, 'k--', 'LineWidth', 1);
xlabel('Fréquence (Hz)'); ylabel('Puissance (dB)');
title('Méthode de Welch');
grid on; ylim([ymin ymax]);

subplot(5,1,5);
plot(f, Pc_dB, 'r', 'LineWidth', 1.2); hold on;
plot(f, Pth_dB, 'k--', 'LineWidth', 1);
xlabel('Fréquence (Hz)'); ylabel('Puissance (dB)');
title('Méthode du Corrélogramme');
grid on; ylim([ymin ymax]);

sgtitle('Comparaison des estimateurs de la DSP du bruit blanc');


%% platitude
nb_real = 100;                  % nombre de réalisations
flatness = zeros(1, nb_real);

for k = 1:nb_real
    x = randn(1, N);             % bruit blanc gaussien
    [P, ~] = Mon_Welch(x, Nfft, Fe,0.5);  
    flatness(k) = platitude_spectrale(P);
end

moyenne = mean(flatness);
ecart_type = std(flatness);

disp(['Platitude spectrale moyenne = ', num2str(moyenne)]);
disp(['Écart-type = ', num2str(ecart_type)]);