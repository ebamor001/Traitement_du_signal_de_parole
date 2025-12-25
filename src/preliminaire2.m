%% Préliminaire 2 : Bruitage du signal de parole
clear; close all; clc;

% Chargement du signal de parole depuis le .mat
load('fcno01fz.mat');   

signal = fcno01fz;

Fe = 8000; 

% Fonction interne pour bruiter 
function [y, bruit] = addNoise(s, RSBdB)
    Ps = mean(s.^2);
    RSB = 10^(RSBdB/10);
    Pn = Ps / RSB;
    bruit = sqrt(Pn) * randn(size(s));
    y = s + bruit;
end

RSB_values = [5 10 15 20];
for i = 1:length(RSB_values)
    [y, bruit] = addNoise(signal, RSB_values(i));

    t = (0:length(signal)-1) / Fe;
   
    %  Affichage 
    figure;
    subplot(2,2,1);
    plot(t,signal);
    title('Signal de parole original');
    xlabel('Temps (s)');ylabel('Amplitude');
    xlim([1 6.5]); 
    A = max(abs(signal));
    ylim([-A A]);
    grid on;


    subplot(2,2,3);
    spectrogram(signal, hamming(256), 128, 512, Fe, 'yaxis');
    title('Spectrogramme - signal original');

    subplot(2,2,2);
    plot(t,y);
    title(['Signal bruité (RSB = ', num2str(RSB_values(i)), ' dB)']);
    xlabel('Temps (s)');ylabel('Amplitude');
    xlim([1 6.5]); 
    A = max(abs(y));
    ylim([-A A]);
    grid on;


    subplot(2,2,4);
    spectrogram(y, hamming(256), 128, 512, Fe, 'yaxis');
    title('Spectrogramme - signal bruité');
end
