clear all;
close all;

[y,Fs] = audioread("note_guitare_LAd.wav");     %Lire le son d'origine
info = audioinfo("note_guitare_LAd.wav");       
duree_son = info.Duration;
N = 160000;

w_hamm = hamming(N);
y_f = y .* w_hamm;      %Appliquer un fenétrage hamming sur le signal d'entrée
Y = fftshift(fft(y_f));


Fmag = abs(Y);
Fphase = angle(Y);

n = linspace(-Fs/2, Fs/2, N);

%% Trouver harmoniques
f0 = 466;
harmoniques = f0*(1:32);

index_harmo = zeros(1, length(harmoniques));
for k = 1:length(harmoniques)
    % Chercher l'index le plus proche
    [~, idx] = min(abs(n - harmoniques(k)));
    
    % Chercher localement le maximum dans une petite fenêtre
    range = max(1, idx-1000):min(N, idx+1000); % Éviter les dépassements
    [~, local_max] = max(Fmag(range));
    
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

%% Passe-bas RIF pour l'enveloppe temporelle

%Ordre du filtre
Fc = pi/1000;
N = 1000;
m = N*Fc/Fs;
K = 2*m+1;

% Génération de la réponse impulsionnelle
k = -N/2:N/2-1; % Indices centrés autour de 0
h = zeros(size(k)); % Initialisation du filtre

% Calcul des coefficients
for i = 1:length(k)
    if k(i) == 0
        h(i) = K / N;
    else
        h(i) = (1/N) * (sin(pi * k(i) * K / N) / sin(pi * k(i) / N));
    end
end

h = hamming(N)'.*h;
y_abs = abs(y);

y_filtered = conv(y_abs, h, 'same');

figure(3);
plot(y_abs, 'b'); hold on;
plot(y_filtered, 'r'); hold off;
xlabel('Échantillon');
ylabel('Magnitude');
legend('Avant filtrage', 'Après filtrage');
title('Signal et enveloppe temporelle');
grid on;

%% Recréer LA#

t = 1:160000;
sum_sinuses = zeros(1, length(t));
for i = 1:length(index_harmo)
    sum_sinuses = sum_sinuses + real(Y(index_harmo(i)))*sin(2*pi*n(index_harmo(i))*t-Fphase(index_harmo(i)));
end

synth = sum_sinuses' .* y_filtered;

figure(4);
plot(synth);
hold on;
plot(y);
hold off;

sound(y, Fs);
%audiowrite('note_guitar_LAd_synthetise.wav', synth, Fs);