clear all
close all
% Importer le son
[y,fe] = audioread("note_basson_filtered_and_downsampled.wav");
N = length(y);
n = linspace(-fe/2, fe/2, N);
%% Trouver les harmoniques du signal Basson et leur magnitude
F = fftshift(fft(y,N));
Fmag = abs(F);
Fphase = angle(F);

f0 = 466/2;                                               % Fréquence fondamentale du son d'origine (La#)
harmoniques = f0*(1:32);                                % Fréquence de chaque harmonique théorique de la fondamentale

index_harmo = zeros(1, length(harmoniques));            %Index des maximums de la fft qui correspondent aux harmoniques théoriques
for k = 1:length(harmoniques)
    % Chercher l'index fréquentiel le plus proche de l'harmonique théorique
    [~, idx] = min(abs(n - harmoniques(k)));
    
    % Chercher localement le maximum dans une petite fenêtre
    range = max(1, idx-1000):min(N, idx+1000);          % Fenétrage pour Éviter les dépassements
    [~, local_max] = max(Fmag(range));                  % Trouver l'index fréquentiel dont la magnitude est la plus haute
    
    % Mettre à jour l'index avec la position du vrai pic
    index_harmo(k) = range(local_max);
end

figure(21);
clf
plot(n, 20*log(Fmag));
hold on;
scatter(n(index_harmo), 20*log(Fmag(index_harmo)), 'ro', 'filled'); 
xlabel('Fréquence Hz');
ylabel('Magnitude (dB)');
title('Détection des harmoniques');
grid on;
legend('FFT', 'Harmoniques détectées');
hold off;
figure(2);
plot(n ,Fphase);
xlabel('Fréquence Hz');
ylabel('Phase');

%% Passe-bas RIF (génération de l'enveloppe)

%Ordre du filtre
Fc = pi/1000;           % Fréquence de coupure
N_lp = 1000;               % Ordre du filtre
m = N_lp*Fc/fe;
K = 2*m+1;

% Génération de la réponse impulsionnelle
k = -N_lp/2:N_lp/2-1; % Indices centrés autour de 0
h = zeros(size(k)); % Initialisation du filtre

% Calcul des coefficients de la réponse impulsionnelle
for i = 1:length(k)
    if k(i) == 0
        h(i) = K / N_lp;
    else
        h(i) = (1/N_lp) * (sin(pi * k(i) * K / N_lp) / sin(pi * k(i) / N_lp));
    end
end

h = hamming(N_lp)'.*h;     % On applique un filtre hamming afin d'éviter l'effet de fuite
y_abs = abs(y);         % On met le signal d'entrée en valeur absolue car l'envloppe temporelle est tjrs au dessus de zéro

y_enveloppe = conv(y_abs, h, 'same');    % On applique le filtre passe bas au signal d'entrée

% Affichage du résultat de l'enveloppe temp
figure(23);
clf
plot(y_abs, 'b'); hold on;
plot(y_enveloppe, 'r'); hold off;
xlabel('Échantillon');
ylabel('Magnitude');
legend('Avant filtrage', 'Après filtrage');
title('Signal et enveloppe temporelle');
grid on;

save("data_basson.mat", "Fmag","Fphase","harmoniques", "index_harmo","y","y_enveloppe","fe","N","n")