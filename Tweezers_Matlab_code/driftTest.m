clear all;
close all;
clc;
load('driftCorrect2.mat');
Fs = 500;
fNyq    =   Fs / 2;
delta_t =   1 / Fs;
X = data.scombined{1}-2;
X = X-mean(X);
Y = data.force{1}-55;
figure;
plot(decimate(X,5),decimate(Y,5));
Tmax = max(data.time)-min(data.time);
N = length(data.time);
dF = Fs/N;
f = (-Fs/2:dF:Fs/2-dF)';
%P = fft(X).*conj(fft(X));
%f = ([0 : length(X)-1] / (max(data.time)+delta_t))';

BPF = 1-((-0.005 < abs(f)) & (abs(f) < 0.005));
figure;
      plot(f,BPF);

      spectrum = fftshift(fft(X))/N;
      spectrumY = fftshift(fft(Y))/N;
      figure;
      subplot(2,1,1);
      plot(f,abs(spectrum));
      
      spectrum = BPF.*spectrum;
      spectrumY = fftshift(fft(Y))/N;
     subplot(2,1,2);
     plot(f,abs(spectrum));
     
     X=ifft(ifftshift(spectrum))*N; %inverse ifft
     %Y=ifft(ifftshift(spectrumY))*N; %inverse ifft
     
     figure;
     
     plot(decimate(X,5),decimate(Y,5));
% ii = find(f<0.01);
% Xfft = fft(X);
% Xfft(ii) = 0;
% X = ifft(Xfft);
%figure;
%loglog(f,P);

% figure;
% plot(decimate(X,50),decimate(Y,50));
