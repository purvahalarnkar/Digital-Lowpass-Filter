clc
clear all
close all

Fs=48000;
Fpass=1000;
Fstop=8000;
Rp=1;
Rs=60;
Wp=Fpass/(Fs/2);
Ws=Fstop/(Fs/2);

% butterworth
[n_butter,Wn_butter]=buttord(Wp,Ws,Rp,Rs);
[num_butter,den_butter]=butter(n_butter,Wn_butter);
Fc_butter=Wn_butter*(Fs/2); 
[h_butter,w_butter]=freqz(num_butter,den_butter,10000);
figure(1)
plot(w_butter/pi*Fs/2,20*log10(abs(h_butter))); 
axis([0 2e3 -5 1])
grid on
title('Frequency response: Butterworth filter')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
figure(2)
zplane(num_butter,den_butter)
title('Pole-zero plot: Butterworth filter')

% chebyshev type 1
[n_cheb1,Wn_cheb1] = cheb1ord(Wp,Ws,Rp,Rs);
[num_cheb1,den_cheb1] = cheby1(n_cheb1,Rp,Wp);
Fc_cheb1=Wn_cheb1*(Fs/2); 
[h_cheb1,w_cheb1]=freqz(num_cheb1,den_cheb1);
figure(3)
plot(w_cheb1/pi*Fs/2,20*log10(abs(h_cheb1)));
axis([0 1.5e3 -5 5]) 
grid on
title('Frequency response: Chebyshev type 1 filter')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
figure(4)
zplane(num_cheb1,den_cheb1)
title('Pole-zero plot: Chebyshev type 1 filter')

%chebyshev type 2
[n_cheb2,Wn_cheb2] = cheb2ord(Wp,Ws,Rp,Rs);
[num_cheb2,den_cheb2] = cheby2(n_cheb2,Rs,Ws);
Fc_cheb2=Wn_cheb2*(Fs/2); 
[h_cheb2,w_cheb2]=freqz(num_cheb2,den_cheb2,10000);
figure(5)
plot(w_cheb2/pi*Fs/2,20*log10(abs(h_cheb2)));
axis([0 18e3 -100 5])
grid on
title('Frequency response: Chebyshev type 2 filter')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
figure(6)
zplane(num_cheb2,den_cheb2)
title('Pole-zero plot: Chebyshev type 2 filter')

%elliptic
[n_ellip,Wn_ellip] = ellipord(Wp,Ws,Rp,Rs);
[num_ellip,den_ellip] = ellip(n_ellip,Rp,Rs,Wp);
Fc_ellip=Wn_ellip*(Fs/2); 
[h_ellip,w_ellip]=freqz(num_ellip,den_ellip,100000);
figure(7)
plot(w_ellip/pi*Fs/2,20*log10(abs(h_ellip)));
axis([0 10e3 -100 5])
grid on
title('Frequency response: Elliptic filter')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
figure(8)
zplane(num_ellip,den_ellip)
title('Pole-zero plot: Elliptic filter')

%linear phase FIR
dev = [(10^(Rp/20)-1)/(10^(Rp/20)+1)  10^(-Rs/20)];
[n_fir,fo,ao,w]=firpmord([Fpass Fstop],[1 0],dev,Fs);
den_fir=1;

%blackman window
num_firB=fir1(n_fir,[Wp Ws],blackman(n_fir+1));
[h_firB,w_firB]=freqz(num_firB,den_fir,10000);
figure(9)
plot(w_firB/pi*Fs/2,20*log10(abs(h_firB)));
grid on
title('Frequency response: Linear-phase FIR filter, Blackmann window')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
figure(10)
zplane(num_firB,den_fir)
title('Pole-zero plot: Linear-phase FIR filter, Blackmann window')

%blackmanharris window
num_firBH=fir1(n_fir,[Wp Ws],blackmanharris(n_fir+1));
[h_firBH,w_firBH]=freqz(num_firBH,den_fir,10000);
figure(11)
plot(w_firBH/pi*Fs/2,20*log10(abs(h_firBH)));
grid on
title('Frequency response: Linear-phase FIR filter, Blackmann Harris window')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
figure(12)
zplane(num_firBH,den_fir)
title('Pole-zero plot: Linear-phase FIR filter, Blackmann Harris window')

%chebyshev window
num_firC=fir1(n_fir,[Wp Ws],chebwin(n_fir+1,50));
[h_firC,w_firC]=freqz(num_firC,den_fir,10000);
figure(13)
plot(w_firC/pi*Fs/2,20*log10(abs(h_firC)));
grid on
title('Frequency response: Linear-phase FIR filter, Chebyshev window')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
figure(14)
zplane(num_firC,den_fir)
title('Pole-zero plot: Linear-phase FIR filter, Chebyshev window')

%Flat top window
num_firF=fir1(n_fir,[Wp Ws],flattopwin(n_fir+1));
[h_firF,w_firF]=freqz(num_firF,den_fir,10000);
figure(15)
plot(w_firF/pi*Fs/2,20*log10(abs(h_firF)));
grid on
title('Frequency response: Linear-phase FIR filter, Flattop window')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
figure(16)
zplane(num_firF,den_fir);
title('Pole-zero plot: Linear-phase FIR filter, Flattop window')

%gaussian window
num_firG=fir1(n_fir,[Wp Ws],gausswin(n_fir+1));
[h_firG,w_firG]=freqz(num_firG,den_fir,10000);
figure(17)
plot(w_firG/pi*Fs/2,20*log10(abs(h_firG)));
grid on
title('Frequency response: Linear-phase FIR filter, Gaussian window')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
figure(18)
zplane(num_firB,den_fir);
title('Pole-zero plot: Linear-phase FIR filter, Gaussian window')

%hanning window
num_firHn=fir1(n_fir,[Wp Ws],hann(n_fir+1));
[h_firHn,w_firHn]=freqz(num_firHn,den_fir,10000);
figure(19)
plot(w_firHn/pi*Fs/2,20*log10(abs(h_firHn)));
grid on
title('Frequency response: Linear-phase FIR filter, Hanning window')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
figure(20)
zplane(num_firHn,den_fir)
title('Pole-zero plot: Linear-phase FIR filter, Hanning window')

%hamming window
num_firHm=fir1(n_fir,[Wp Ws],hamming(n_fir+1));
[h_firHm,w_firHm]=freqz(num_firHm,den_fir,10000);
figure(21)
plot(w_firHm/pi*Fs/2,20*log10(abs(h_firHm)));
grid on
title('Frequency response: Linear-phase FIR filter, Hamming window')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
figure(22)
zplane(num_firHm,den_fir)
title('Pole-zero plot: Linear-phase FIR filter, Hamming window')

%kaiser window
num_firK=fir1(n_fir,[Wp Ws],kaiser(n_fir+1));
[h_firK,w_firK]=freqz(num_firK,den_fir,10000);
figure(23)
plot(w_firK/pi*Fs/2,20*log10(abs(h_firK)));
grid on
title('Frequency response: Linear-phase FIR filter, Kaiser window')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
figure(24)
zplane(num_firK,den_fir)
title('Pole-zero plot: Linear-phase FIR filter, Kaiser window')

%rectangular window
num_firR=fir1(n_fir,[Wp Ws],rectwin(n_fir+1));
[h_firR,w_firR]=freqz(num_firR,den_fir,10000);
figure(25)
plot(w_firR/pi*Fs/2,20*log10(abs(h_firR)));
grid on
title('Frequency response: Linear-phase FIR filter, Rectangular window')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
figure(26)
zplane(num_firR,den_fir)
title('Pole-zero plot: Linear-phase FIR filter, Rectangular window')

%tukey window
num_firT=fir1(n_fir,[Wp Ws],tukeywin(n_fir+1));
[h_firT,w_firT]=freqz(num_firT,den_fir,10000);
figure(27)
plot(w_firT/pi*Fs/2,20*log10(abs(h_firT)));
grid on
title('Frequency response: Linear-phase FIR filter, Tukey window')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
figure(28)
zplane(num_firT,den_fir)
title('Pole-zero plot: Linear-phase FIR filter, Tukey window')

%triangular window
num_firTr=fir1(n_fir,[Wp Ws],triang(n_fir+1));
[h_firTr,w_firTr]=freqz(num_firTr,den_fir,10000);
figure(29)
plot(w_firTr/pi*Fs/2,20*log10(abs(h_firTr)));
grid on
title('Frequency response: Linear-phase FIR filter, Triangular window')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
figure(30)
zplane(num_firTr,den_fir)
title('Pole-zero plot: Linear-phase FIR filter, Triangular window')
