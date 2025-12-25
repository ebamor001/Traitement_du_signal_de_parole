function [P_bart, f] = Mon_Bartlett(x, Nfft, Fe)
    N = length(x);
    nb_segments = ceil(N / Nfft);             % nombre total de segments (y compris le dernier incomplet)
    P_bart = zeros(1, Nfft);                  % initialisation de la DSP moyenne

    % Zero-padding pour compléter le dernier segment
    x = [x, zeros(1, nb_segments * Nfft - N)];

    % On ne garde que les segments complets
    for i = 1:nb_segments
        seg = x((i-1)*Nfft + 1 : i*Nfft);
        [P, ~] = Mon_Periodogramme(seg, Nfft, Fe);
        P_bart = P_bart + P;
    end

    P_bart=P_bart/ nb_segments;
    f = linspace(-Fe/2, Fe/2, Nfft);
end


    
   