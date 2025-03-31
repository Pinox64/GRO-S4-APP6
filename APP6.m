[y,Fs] = audioread("note_guitare_LAd.wav");
N = 160000;
F = fftshift(fft(Fs, N));

Fmag = abs(F);
Fphase = angle(F);

figure(1);
plot(N ,20*log(Fmag));
figure(2);
plot(N ,Fphase);