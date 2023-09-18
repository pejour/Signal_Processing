% TP 2
% Exercice 1
N =100;
Fe = 10000;
Te = 1/Fe;
f1 = 1000;
f2 = 3000;
intervalle = [0:Te:(N-1)*Te];
X = cos(2*pi*f1*intervalle);
Y = cos(2*pi*f2*intervalle);
Z= X+Y;
plot (intervalle,X);



frequ = [0:(Fe/N):Fe-(Fe/N)];
TF = abs(fft(Z));
figure;
subplot(211);
semilogy(frequ,TF);
xlabel('frequence (Hz)');
ylabel('log(TF(|X(t)|))');
subplot(212);
plot(frequ,TF);
xlabel('frequence(Hz)');
ylabel('TF(|X(t)|)');


% Exercice 2:

fc = 2000;
intervalle2 = [-5*Te:Te:5*Te];
B = 2*fc*sinc(2*fc*intervalle2);
s = filter(B,1,X);
figure;
plot(intervalle,s);