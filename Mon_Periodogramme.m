function [P, f] = Mon_Periodogramme(x, Nfft, Fe)
    %   Nfft : taille de la FFT (une puissance de 2)
    N = length(x);
    X = fft(x, Nfft);
    P = (abs(X).^2)/ (Fe * N);       % densité spectrale de puissance
    P = fftshift(P);                 % centrage en fréquence
    f = linspace(-Fe/2, Fe/2, Nfft);
end
