function F = platitude_spectrale(P)
    % Calcule la platitude spectrale d'un spectre P
    P = P(:);                 % vecteur colonne
    P(P <= 0) = eps;          % éviter log(0)
    geo_mean = exp(mean(log(P)));
    ari_mean = mean(P);
    F = geo_mean / ari_mean;
end