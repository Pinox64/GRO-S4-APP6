%% ========================================================================
%  GRO430 – Traitement numérique des signaux (APP6)
%  Script : APP6.m
%  Étudiants : Samuel Hamelin   &   Renaud Gagnon
%  Université de Sherbrooke – Faculté de génie
%  Date : <2025-04-04>
%  ------------------------------------------------------------------------

clear all;
close all;

% Importer l'audio
[y,Fs] = audioread("note_guitare_LAd.wav");
info = audioinfo("note_guitare_LAd.wav");
duree_son = info.Duration;
N = 160000;

w_hamm = hamming(N);
y_f = y .* w_hamm;
Y = fftshift(fft(y_f));

Fmag = abs(Y);
Fphase = angle(Y);

n = linspace(-Fs/2, Fs/2, N);

%% Trouver les harmoniques et leur magnitude
f0 = 466;                                               % Fréquence fondamentale du son d'origine (La#)
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

figure(1);
plot(n, 20*log10(Fmag));
hold on;
scatter(n(index_harmo), 20*log10(Fmag(index_harmo)), 'ro', 'filled'); 
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
N = 1000;               % Ordre du filtre
m = N*Fc/Fs;
K = 2*m+1;

% Génération de la réponse impulsionnelle
k = -N/2:N/2-1; % Indices centrés autour de 0
h = zeros(size(k)); % Initialisation du filtre

% Calcul des coefficients de la réponse impulsionnelle
for i = 1:length(k)
    if k(i) == 0
        h(i) = K / N;
    else
        h(i) = (1/N) * (sin(pi * k(i) * K / N) / sin(pi * k(i) / N));
    end
end

h = hamming(N)'.*h;     % On applique un filtre hamming afin d'éviter l'effet de fuite
y_abs = abs(y);         % On met le signal d'entrée en valeur absolue car l'envloppe temporelle est tjrs au dessus de zéro

y_filtered = conv(y_abs, h, 'same');    % On applique le filtre passe bas au signal d'entrée

% Affichage du résultat de l'enveloppe temp
figure(3);
plot(y_abs, 'b'); hold on;
plot(y_filtered, 'r'); hold off;
xlabel('Échantillon');
ylabel('Magnitude');
legend('Avant filtrage', 'Après filtrage');
title('Signal et enveloppe temporelle');
grid on;

%% Synthétiser les notes nécessaires à la mélodie que l'on veut jouer
% Fréquence de chaque note (en Hz)
freq.DO   = 261.6;
freq.DOd  = 277.2;
freq.RE   = 293.7;
freq.REd  = 311.1;
freq.MI   = 329.6;
freq.FA   = 349.2;
freq.FAd  = 370.0;
freq.SOL  = 392.0;
freq.SOLd = 415.3;
freq.LA   = 440.0;
freq.LAd  = 466.2;
freq.SI   = 493.9;

% Facteur a appliquer a in LAd pour obtenir chaque note
fact.DO   = freq.DO   / freq.LAd;
fact.DOd  = freq.DOd  / freq.LAd;
fact.RE   = freq.RE   / freq.LAd;
fact.REd  = freq.REd  / freq.LAd;
fact.MI   = freq.MI   / freq.LAd;
fact.FA   = freq.FA   / freq.LAd;
fact.FAd  = freq.FAd  / freq.LAd;
fact.SOL  = freq.SOL  / freq.LAd;
fact.SOLd = freq.SOLd / freq.LAd;
fact.LA   = freq.LA   / freq.LAd;
fact.LAd  = freq.LAd  / freq.LAd;  % = 1
fact.SI   = freq.SI   / freq.LAd;

%%Recréer un LA#
freqLad = 466.2;
t = (0:length(y_filtered)-1) / Fs;

sum_sinuses = zeros(1, length(t));
for i = 1:length(index_harmo)
    sum_sinuses = sum_sinuses + Fmag(index_harmo(i)) * cos(2*pi*harmoniques(i)*t+Fphase(index_harmo(i)));
end

synthLA = sum_sinuses' .* y_filtered;
synthLA = synthLA / max(abs(synthLA));

%sound(synthLA, Fs);

%%Recréer un SOL
t = (0:length(y_filtered)-1) / Fs;
SOL_harmoniques = harmoniques.*fact.SOL;
sum_sinuses = zeros(1, length(t));
for i = 1:length(index_harmo)
    sum_sinuses = sum_sinuses + Fmag(index_harmo(i)) * cos(2*pi*SOL_harmoniques(i)*t+Fphase(index_harmo(i)));
end

synthSOL = sum_sinuses' .* y_filtered;
synthSOL = synthSOL / max(abs(synthSOL));

%%Recréer MI
t = (0:length(y_filtered)-1) / Fs;
MI_harmoniques = harmoniques.*fact.MI;
sum_sinuses = zeros(1, length(t));
for i = 1:length(index_harmo)
    sum_sinuses = sum_sinuses + Fmag(index_harmo(i)) * cos(2*pi*MI_harmoniques(i)*t+Fphase(index_harmo(i)));
end

synthMI = sum_sinuses' .* y_filtered;
synthMI = synthMI / max(abs(synthMI));

%%Recréer FA
t = (0:length(y_filtered)-1) / Fs;

FA_harmoniques = harmoniques.*fact.FA;
sum_sinuses = zeros(1, length(t));
for i = 1:length(index_harmo)
    sum_sinuses = sum_sinuses + Fmag(index_harmo(i)) * cos(2*pi*FA_harmoniques(i)*t+Fphase(index_harmo(i)));
end

synthFA = sum_sinuses' .* y_filtered;
synthFA = synthFA / max(abs(synthFA));

%%Recréer un RE
t = (0:length(y_filtered)-1) / Fs;

RE_harmoniques = harmoniques.*fact.RE;
sum_sinuses = zeros(1, length(t));
for i = 1:length(index_harmo)
    sum_sinuses = sum_sinuses + Fmag(index_harmo(i)) * cos(2*pi*RE_harmoniques(i)*t+Fphase(index_harmo(i)));
end

synthRE = sum_sinuses' .* y_filtered;
synthRE = synthRE / max(abs(synthRE));

%% Jouer mélodie complète

T = timer('TimerFcn',@(~,~)disp(''),'StartDelay',0.25);
T2 = timer('TimerFcn',@(~,~)disp(''),'StartDelay',2);
start(T);
wait(T);
stop(T);
sound(synthSOL, Fs);
start(T);
wait(T);
stop(T);
sound(synthSOL, Fs);
start(T);
wait(T);
stop(T);
sound(synthSOL, Fs);
start(T);
wait(T);
stop(T);
sound(synthMI, Fs);
start(T2);
wait(T2);
stop(T2);
%clear sound
sound(synthFA, Fs);
start(T);
wait(T);
stop(T);
sound(synthFA, Fs);
start(T);
wait(T);
stop(T);
sound(synthFA, Fs);
start(T);
wait(T);
stop(T);
sound(synthRE, Fs);

