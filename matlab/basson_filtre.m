clear all
close all
%But: filtrer le bruit et la sinus.
%Pour filtrer le bruit on veut faire un RII un passe bas pour eliminer +15khz
%Pour filtrer la sinus, on veut faire un Coupe bande butterworth pour
%eliminer tout ce qui est entre 980hz à 1020hz
[x,fe] = audioread("note_basson_plus_sinus_1000_Hz_plus_hautes_freqs.wav");
 %% Filtre coupe bande

fc_lp = 40;  % Fréquence de coupure du filtre
N = 1024;   %Ordre du filtre (Donné dans l'énoncé)
m = N*fc_lp/fe
K = m*2+1;
fc_nf = 1000;
w0 = 2*pi*fc_nf/fe; %Fréqence centrale (en Hz) de la bande à élminier (
% A small local function that never computes 0/0 for x=0
function val = h_bas(x, K, N)
    if x == 0
        val = K/N;
    else
        val = (1/N) * sin(pi*K*x/N) / sin(pi*x/N);
    end
end
h_bas2 = @(n) arrayfun(@(x) h_bas(x, K, N), n);
function val = h_bande(x, K, N,w0)
    delta = double(x==0); %Dirac discret
    val = delta - 2*h_bas(x,K,N)*cos(w0*x);
end
% Then define an anonymous "vectorized" function using arrayfun:
h = @(n) arrayfun(@(x) h_bande(x, K, N,w0), n);

Nfft = N;
n_ = -N/2 : N/2-1; 
h_ = h(n_);   % No NaNs, h_ is a vector
hb_ = h_bas2(n_); 
h_fft = fftshift(fft(h_,Nfft));
h_mag = 20*log(abs(h_fft));
h_phase = angle(h_fft);

hb_fft = fftshift(fft(hb_,Nfft));
hb_mag = 20*log(abs(hb_fft));
hb_phase = angle(hb_fft);

f_ = linspace(-fe/2, fe/2, Nfft);   % axe des fréquences en Hz

figure("Name","reponse freq du passe bas")
clf
%freqz(h_,1024);
%hold on;
%freqz(h_b,1024);
subplot(2,1,1);
plot(f_,hb_mag);
title("mag");
subplot(2,1,2);
plot(f_,hb_phase);
title("phase");
xlabel("Frequence");
ylabel("Amplitude (db)");
%legend("Coupe bande", "Passe bas");

figure("Name","reponse freq du coupe bande")
clf
%freqz(h_,1024);
%hold on;
%freqz(h_b,1024);
subplot(2,1,1);
plot(f_,h_mag);
title("mag");
subplot(2,1,2);
plot(f_,h_phase);
title("phase");
legend("Coupe bande", "Passe bas");

% Affichage des poles et des zero du coupe band
figure ("Name", "Pôles et zéros du filtre RIF coupe bande");
zplane(h_,1);
title("Pôles et zéros du filtre RIF coupe bande");

figure ("Name","Coupe bande");
clf
plot(n_,h_);
title("Réponse a une impulsion du coupe bande")
%h_trunc = h_(end/2, end)
h_mogus = hamming(Nfft)'.*h_;
h_reconstructed = real(ifft(ifftshift(h_)));
filteredsinSound = conv(conv(x,h_mogus, "same"),h_mogus,"same");
filteredsinSound = filteredsinSound(length(h_mogus):length(filteredsinSound)-length(h_mogus));
audiowrite("note_basson_plus_sinus_1000_Hz_plus_hautes_freqs_NoSin.wav",filteredsinSound,fe)

%% Filtre RII (butterworth) pour filtrer ce qui est plus haut que 15KHz

%fc = 12000;
fp = 8000; % Fréquence de passage (on coupe moins de 1db)
fs = 15000; % Fréquence qu'on veut complètement couper
ny = fe/2; % Fréquence de nyquist du signal
Wp = fp/ny; % Fréquence de passage normalisée selon nyquist de fe
Ws = fs/ny; % Fréquence de stop normalisée selon nyquist de fe

Rp = 1; % Ripple minimal dans la fréquence de passage (en db)
Rs = 40; % Atténuation minimale à la fréquence de stop (en db)

[N, Wn] = buttord(Wp,Ws,Rp,Rs)
[b,a] = butter(N,Wn,"low");
 fullyFiltered = filter(b,a,filteredsinSound);
 audiowrite("note_basson_plus_sinus_1000_Hz_plus_hautes_freqs_hf_couper.wav",fullyFiltered,fe);
 %figure (1);
 %t = 1:length(x);
 %plot(t,x);
 %hold on;
 %plot(t,fullyFiltered);
 %legend("in","out");

 % Affichage des poles et des zero du filtre valant du beurre
 figure ("Name", "Pôles et zéros du filtre RII (butterworth)");
 zplane(b,a);
 title("Pôles et zéros du filtre RII (butterworth)");

y = fullyFiltered;
w_hamm = hamming(length(y));
y_f = y .* w_hamm;      %Appliquer un fenétrage hamming sur le signal d'entrée
Y = fftshift(fft(y_f));


%% Rééchantillonnement du signal Basson
y_downsampled = downsample(fullyFiltered,2);
fe2 = fe/2
N = length(y_downsampled);
n = linspace(-fe2/2, fe2/2, N);

sound(y_downsampled, fe2)
audiowrite("note_basson_filtered_and_downsampled.wav", y_downsampled, fe2)