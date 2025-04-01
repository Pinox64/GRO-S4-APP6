close all
clear all

%But: filtrer le bruit et la sinus.
%Pour filtrer le bruit on veut faire un RII un passe bas pour eliminer +15khz
%Pour filtrer la sinus, on veut faire un Coupe bande butterworth pour
%eliminer tout ce qui est entre 980hz à 1020hz

%% Filtre RII (butterworth) pour filtrer ce qui est plus haut que 15KHz
[x,fe] = audioread("note_basson_plus_sinus_1000_Hz_plus_hautes_freqs.wav");
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
 y = filter(b,a,x);
 
 audiowrite("note_basson_plus_sinus_1000_Hz_plus_hautes_freqs_hf_couper.wav",y,fe);
 figure (1);
 t = 1:length(x);
 plot(t,x);
 hold on;
 plot(t,y);
 legend("in","out");

 %% Filtre coupe bande

fc = 1000;  % Fréquence de coupure du filtre
N = 1024; %Ordre du filtre (Donné dans l'énoncé)
m = N*fc/fe;
K = m*2+1
w0 = 2*pi*fc/fe; %Fréqence centrale (en Hz) de la bande à élminier (
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
    delta = double(x==0) %Dirac discret
    val = delta - 2*h_bas(x,K,N)*cos(w0*x);
end
% Then define an anonymous "vectorized" function using arrayfun:
h = @(n) arrayfun(@(x) h_bande(x, K, N,w0), n);

n_ = -N/2 : N/2-1; 
h_ = h(n_);   % No NaNs, h_ is a vector
h_fft = fftshift(fft(h_,1024));
h_mag = 20*log(abs(h_fft));
h_phase = angle(h_fft);
freqz(h_);
figure ("Name","Coupe bande");
plot(n_,h_);
title("Réponse fréquentielle du coupe bande")

filteredsinSound = conv(y,h_);
audiowrite("note_basson_plus_sinus_1000_Hz_plus_hautes_freqs_fully_filtered.wav",filteredsinSound,fe)
