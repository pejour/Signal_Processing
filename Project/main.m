clear;
close all;


% Récupération des données à transmettre.

full_data = load ("DonneesBinome4.mat").bits;
data_1000 = full_data(1:1000);

DATA = data_1000;


% Initialisation des constantes

Fs = 300;   % bits/s
Ts = 1 / Fs;
Fe = 48000;
Te = 1 / Fe;
Ns = floor (Ts / Te);


nb_data = length (DATA);
duree = Ts * nb_data;
nb_echantillon = floor (duree / Te);


F0 = 6000;
F1 = 2000;


SNR = 2;


K = 57;



%% Signal NRZ

% Récupération du signal NRZ

Temps = 0:Te:duree;

NRZ = zeros (1, nb_echantillon);
for i = 1:length(Temps)
    t = Temps(i);
    data_indice = floor (t / Ts) + 1;
    if data_indice <= nb_data
        NRZ(i) = DATA(data_indice);
    else
        NRZ(i) = DATA(nb_data);
    end
end

Ne_ = 2^nextpow2(length(NRZ));

figure();

plot (Temps(1:floor(Ns*10)), NRZ(1:floor(Ns*10)));
xlabel ("Temps en secondes");
ylabel ("NRZ (t)");
title ("Signal NRZ en fonction du temps");



% Calcule de la densité spectrale de puissance

Frequences = 0:Fe/Ne_:(Ne_-1)*Fe/Ne_;

Snrz_reelle = 1/nb_echantillon * (abs(fft(NRZ, Ne_))).^2;
Snrz_theorique = (Ts * sinc (pi * Frequences * Ts).^2 + dirac (Frequences)) / 4;
Snrz_theorique = Snrz_theorique + flip (Snrz_theorique);


figure();

semilogy (Frequences, Snrz_reelle, Frequences, Snrz_theorique);
xlabel ("Fréquences en hertz");
ylabel ("Snrz (f)");
title ("Densité spectrale de puissance du signal NRZ");


%% Signal modulé x (t)

% Génération du signal modulé

Cos0 = cos (2 * pi * F0 * Temps);
Cos1 = cos (2 * pi * F1 * Temps);
X = (1 - NRZ) .* Cos0 + NRZ .* Cos1;


figure();

plot (Temps(1:floor(Ns*10)), X(1:floor(Ns*10)));
xlabel ("Temps en secondes");
ylabel ("x (t)");
title ("Signal modulé en fréquence");



%% Densité spectrale de puissance de x (t) théorique

Rx = zeros (1, nb_echantillon);
for i = 1:nb_echantillon
    Rx(i) = mean (X(1+i:length(X)) .* X(1:length(X)-i));
end
Sx_theorique = fft(Rx, Ne_).^2;


% Calculer la densité spectrale de puissance de x (t) réelle

Sx_reelle = 1/nb_echantillon * (abs(fft(X, Ne_))).^2;


figure ();

subplot(211);
semilogy (Frequences, Sx_reelle);
xlabel ("Fréquences en hertz");
ylabel ("Sx (f)");
title ("Densité spectale de puissance de x (t) (réelle)");

subplot(212);
semilogy (Frequences, Sx_theorique);
xlabel ("Fréquences en hertz");
ylabel ("Sx (f)");
title ("Densité spectale de puissance de x (t) (théorique)");


%% Generation du bruit blanc

% Génération du bruit blanc

Px = mean (abs(X).^2);
sigma = Px * 10^(-SNR / 10);
bruit = sigma * randn (1, length(X));

X_bruite = X + bruit;

Sx_bruite = 1/nb_echantillon * (abs(fft(X_bruite, Ne_))).^2;

figure ();

plot (Temps(1:floor(Ns*10)), X_bruite(1:floor(Ns*10)));
xlabel ("Temps en secondes");
ylabel ("x (t)");
title ("Signal x (t) bruité");


%% Filtrages

% Filtrage passe-bas

fc = (F0 + F1) / 2;
fc_neg = (2 * Fe - F0 - F1) / 2;
largeur = 2 * fc;
ordre = 11;
intervalle = -(ordre-1) * Te / 2 :Te: (ordre-1) * Te / 2;
B = largeur * sinc(largeur * intervalle);

x_filtre_pb = filter(B, 1, [X_bruite zeros(1, (ordre-1) / 2)]);
x_filtre_pb_padding = x_filtre_pb((ordre-1)/2:length(x_filtre_pb));

Sx_pb = 1/nb_echantillon * (abs(fft(x_filtre_pb_padding, Ne_))).^2;


figure ();

plot (Temps(1:floor(Ns*10)), x_filtre_pb_padding(1:floor(Ns*10)));
xlabel ("Temps en secondes");
ylabel ("x~ (t)");
title ("Signal x (t) filtré (filtre passe-bas)");


% filtrage passe-haut

indice_fc = floor (fc * Ne_ / Fe);
indice_fc_neg = floor (fc_neg * Ne_ / Fe);
HIpb = [ones(1, indice_fc) zeros(1, indice_fc_neg - indice_fc) ones(1, length(Frequences) - indice_fc_neg)];

HIph = (1 - HIpb);

H = ifft (HIph);

x_filtre_ph = filter (H, 1, [X_bruite zeros(1, (ordre-1) / 2)]);
x_filtre_ph_padding = x_filtre_ph((ordre-1)/2:length(x_filtre_ph));

Sx_ph = 1/nb_echantillon * (abs(fft(x_filtre_ph_padding, Ne_))).^2;



figure ();

plot (Temps(1:floor(Ns*10)), x_filtre_ph_padding(1:floor(Ns*10)));
xlabel ("Temps en secondes");
ylabel ("x~ (t)");
title ("Signal x (t) filtré (filtre passe-haut)");



%% Réponse en fréquence

figure ();

subplot (211);
plot (Temps, H(1:length(Temps)));
xlabel ("Temps en secondes");
ylabel ("Réponse impulsionnelle");
title ("Réponse impulsionnelle du filtre passe-haut");

subplot (212);
plot (Frequences, HIph);
xlabel ("Fréquences en Hz");
ylabel ("Réponse en fréquence");
title ("Réponse en fréquence du filtre pass-haut");


figure ();

subplot (211);
plot (Frequences, 10000*HIph, Frequences, Sx_bruite);  % Pour plus de visibilité, nous avons multiplié la réponse en fréquence par 10 000.
xlabel ("Fréquences en Hz");
ylabel ("Densité spectrale de puissance");
title ("Densité spectrale en entrée du filtre");

subplot (212);
plot (Frequences, Sx_ph);
xlabel ("Fréquences en Hz");
ylabel ("Densité spectrale de puissance");
title ("Densité spectrale en sortie du filtre");




%% Détection d'énergie

suite_binaire_reconstruite =  zeros (nb_data);
energie = zeros (nb_data);
nb_erronnes = 0;

for i = 1:nb_data
    energie(i) = sum (x_filtre_ph_padding ((i-1)*Ns+1:i*Ns) .^ 2);
    if energie(i) < K
        % detection d'énergie
         suite_binaire_reconstruite(i) = 1;
        
        % rapport bit erronnés
        if DATA(i) ~= 1
            nb_erronnes = nb_erronnes + 1;
        end
    elseif DATA(i) ~= 0
        nb_erronnes = nb_erronnes + 1;
    end
end

taux_erreur_binaire = nb_erronnes / nb_data * 100


% Affichage de l'image reconstituée

%pcode reconstitution_image;
%reconstitution_image (suite_binaire_reconstruite);
%which reconstitution_image;


%% Démodulateur FSK

% Sans gestion de synchronisation

x_filtre_resized = x_filtre_ph_padding(2:length(x_filtre_ph_padding));

prod_F0 = x_filtre_resized .* cos (2*pi*F0*Temps);
prod_F1 = x_filtre_resized .* cos (2*pi*F1*Temps);

suite_binaire_reconstruite_FSK =  zeros (nb_data);
integrale = zeros (nb_data);
nb_erronnes_FSK = 0;

for i = 1:nb_data
    integrale(i) = sum (prod_F1((i-1)*Ns+1:i*Ns) - prod_F0((i-1)*Ns+1:i*Ns));
    if integrale(i) < 0
        % detection intégrale positive
         suite_binaire_reconstruite_FSK(i) = 1;
        
        % rapport bit erronnés
        if DATA(i) ~= 1
            nb_erronnes_FSK = nb_erronnes_FSK + 1;
        end
    elseif DATA(i) ~= 0
        nb_erronnes_FSK = nb_erronnes_FSK + 1;
    end
end

taux_erreur_binaire_FSK = nb_erronnes_FSK / nb_data * 100


% Gestion d'une erreur de synchronisation de phase de porteuse

prod_cosF0 = x_filtre_resized .* cos (2*pi*F0*Temps);
prod_cosF1 = x_filtre_resized .* cos (2*pi*F1*Temps);
prod_sinF0 = x_filtre_resized .* sin (2*pi*F0*Temps);
prod_sinF1 = x_filtre_resized .* sin (2*pi*F1*Temps);

suite_binaire_reconstruite_FSK_synchro =  zeros (nb_data);
integrale_synchro = zeros (nb_data);
nb_erronnes_FSK_synchro = 0;

for i = 1:nb_data
    integrale_synchro(i) = sum (prod_cosF1((i-1)*Ns+1:i*Ns))^2 + sum (prod_sinF1((i-1)*Ns+1:i*Ns))^2 - sum (prod_cosF0((i-1)*Ns+1:i*Ns))^2 - sum (prod_sinF0((i-1)*Ns+1:i*Ns))^2;
    if integrale_synchro(i) > 0
        % detection intégrale positive
         suite_binaire_reconstruite_FSK_synchro(i) = 1;
        
        % rapport bit erronnés
        if DATA(i) ~= 1
            nb_erronnes_FSK_synchro = nb_erronnes_FSK_synchro + 1;
        end
    elseif DATA(i) ~= 0
        nb_erronnes_FSK_synchro = nb_erronnes_FSK_synchro + 1;
    end
end

taux_erreur_binaire_FSK_synchro = nb_erronnes_FSK_synchro / nb_data * 100