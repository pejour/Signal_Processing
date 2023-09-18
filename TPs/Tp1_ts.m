% TP 1
% Exercice1:
% q1&2
fe1 = 10000;
Te = 1/fe1;
f0 = 1100;
N = 90;
subplot(211);
intervalle = [0:Te:(N-1)*Te];
X = cos(2*pi*f0*intervalle);
plot(intervalle,X);
xlabel('Temps (en s)');
ylabel('cosinus');

% q3&4
subplot(212);
Te2 = 1/1000;
intervalle2 = [0:Te2:(N-1)*Te2];
Y = cos(2*pi*f0*intervalle2);
plot(intervalle2,Y);
xlabel('Temps (en s)');
ylabel('cosinus 2');



% Exercice 2:
% q1
frequ = [0:(fe1/N):fe1-(fe1/N)];
TF = abs(fft(X));
figure;
subplot(211);
semilogy(frequ,TF);
xlabel('frequence (Hz)');
ylabel('log(TF(|X(t)|))');
subplot(212);
plot(frequ,TF);
xlabel('frequence(Hz)');
ylabel('TF(|X(t)|)');
