clear all
close all
clc

SNR_in = 10; % pas en dB, en lineaire

%Chargement signaux
load 1z88153a8.mat
signal=sigd;

load white8.mat 
noise=sigd(1:length(signal));

clear sigd

Esig = signal' * signal;
%Esig = rms(signal);
Enoise = noise' * noise;
%Enoise = rms(noise);

signal_bruite = signal + noise*sqrt(Esig/(SNR_in*Enoise));

noise = noise*sqrt(Esig/(SNR_in*Enoise));

%%Normalisation signaux
signal = signal/max(abs(signal_bruite));

noise = noise/max(abs(signal_bruite));

signal_bruite = signal_bruite/max(abs(signal_bruite));

sound(signal_bruite)

%Variables de traitement

Fs = 8000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(signal);   % Length of signal
t = (0:L-1)*T;        % Time vector
LSample = 256;        % Nombre d'echantillon des trames echantillonées (temps des frames = 256/8000 = 32ms)
FFT_Ordre = 256;      % Ordre de la FFT, ne doit pas etre inferieur a LSample
Beta=0.02;            % Beta
Alpha=5;              % Alpha
Gamma = 1;         % Exposant utilise pour la puissance spectrale
G = 1;                %Facteur de Normalisation
LSample_ov = 228;     % nombre d'échantillons par trame overlappées (228 présent, 14 passés et 14 futurs)
index_Bruit_SNR = 4;

%Affiche le Signal, le Bruit et le Signal Bruite
figure
subplot(3,1,1)
plot(t,signal); title('Signal original'); 
subplot(3,1,2)
plot(t,noise); title('Signal de bruit'); 
subplot(3,1,3)
plot(t,signal_bruite); title('Signal bruité'); xlabel('Temps (s)');

% Transformee de Fourrier du signal initial
P2 = fftshift(fft(signal)); 
f2 = Fs*(-L/2:L/2-1)/L;

figure
subplot(2,1,1)
plot(f2,abs(P2).^2);
title('Puissance du signal intial');
xlabel('Fréquence (en Hz)');
subplot(2,1,2)
plot(f2,abs(P2));
title('Représentation Fréquentielle du signal intial'); xlabel('Fréquence (en Hz)');

%Creation du bruit moyen
Bruit=0;
j=0;
for i=1:7
    E_Bruit = signal_bruite((i-1)*LSample+1:i*LSample);
    Bruit = Bruit + abs(fftshift(fft(E_Bruit,FFT_Ordre))).^2;
    j=j+1;
end
PS_Bruit = Bruit/i;


%F_Bruit = fft(Bruit);
%PS_Bruit = abs(F_Bruit).^2;
f3 = Fs*(-LSample/2:LSample/2-1)/LSample;
%Affiche le Bruit
figure;
% subplot(211)
plot(f3,fftshift(Bruit));
title('Représentation de la puissance moyenne du Bruit'); xlabel('Fréquence (Hz)');
% subplot(212)
% plot(signal_bruite((i-1)*LSample+1:i*LSample));
% title('Dernier échantillon utilisé'); xlabel('Temps (s)'); 

index=15;
%Exemple des valeurs des différentes composantes lors d'un calcul de trames
Echantillon = signal_bruite((index-1)*LSample+1:index*LSample);
F_Echantillon = fftshift(fft(Echantillon));
PS_Echantillon = abs(F_Echantillon(1:LSample)).^2;
Echantillon_NB = signal((index-1)*LSample+1:index*LSample);
F_Echantillon_NB = fftshift(fft(Echantillon_NB));
PS_Echantillon_NB = abs(F_Echantillon_NB(1:LSample)).^2;

figure
plot(f3, PS_Echantillon - Alpha*PS_Bruit, 'Color', 'm', 'DisplayName','PS Echantillon - PS Bruit'); title('Visualisation des différentes puissances spectrales manipulées avec Alpha = 5 dans une trame de parole')
legend('-DynamicLegend');
hold on
plot(f3, PS_Echantillon, 'Color', 'g', 'DisplayName','PS Echantillon');
plot(f3, PS_Bruit, 'Color', 'r', 'DisplayName','PS Bruit');
plot(f3, PS_Echantillon_NB,'b', 'DisplayName','PS Echantillon Initial');xlabel('Fréquence (en Hz)');
hold off

%Calcul du nouveau signal Debruitage simple
max_index = floor(length(signal)/LSample);
Signal_Debruite = 0;
for index=1:max_index-1
    Echantillon = signal_bruite((index-1)*LSample+1:index*LSample); %Identification de l'echantillon
    F_Echantillon = fftshift(fft(Echantillon)); 
    Argument_E_Signal = angle(F_Echantillon); %Recuperation de la phase
    PS_Echantillon = abs(F_Echantillon).^2;
    PS_Echantillon_Debruite = PS_Echantillon - PS_Bruit; %Operation de Debruitage
    for j=1:LSample
        if PS_Echantillon_Debruite(j) < 0 %suppression des valeurs resultantes negatives (tel que mentionne dans l'article)
            PS_Echantillon_Debruite(j) = 0.0000000; 
        end
    end
    Spectre_E_Debruite = sqrt(PS_Echantillon_Debruite);
    Echantillon_Debruite = complex(Spectre_E_Debruite.*cos(Argument_E_Signal),Spectre_E_Debruite.*sin(Argument_E_Signal)); % Reintegration des phases avec les valeurs des modules calculees
    Signal_Debruite = cat(1,Signal_Debruite,ifft(ifftshift(Echantillon_Debruite))); % Construction du signal final
end
pause(1.5);
sound(real(Signal_Debruite)); 

%Calcul du nouveau signal avec la solution proposee
Signal_DebruiteArticle = 0;

% Détermination de alpha0
alpha0=0;
for index=1:7 %Pour chaque trame du début qui ne comporte que du bruit.
                % On utilise ces trames car dans ces trames Signal=Bruit.
                % On a donc SNR = log (Bruit/Bruit) = log(1) = 0.
    Echantillon = signal_bruite((index-1)*LSample+1:index*LSample); %Identification de l'échantillon
    PS_Echantillon = abs(fftshift(fft(Echantillon))).^2;
    for alpha=1:1:50 %On cherche le alpha qui permet de supprimer le bruit, ce qui equivaut a ce que tous les echantillons passe sous le seuil definit par Beta
        PS_Echantillon_Debruite=PS_Echantillon - alpha*PS_Bruit;
        for j=1:LSample
            if PS_Echantillon_Debruite(j) < Beta*PS_Bruit(j)
                PS_Echantillon_Debruite(j) = Beta*PS_Bruit(j); 
            end
        end
        if mean(PS_Echantillon_Debruite) == mean(Beta*PS_Bruit)
            break;
        end
    end
    alpha0=alpha0+alpha;
end
%Alpha0 correspond a la moyenne des alphas qui pour chaque trame du début supprime le signal (qui correspond au bruit pour ces trames).
alpha0=alpha0/index;
disp(['Alpha0 = ', num2str(alpha0)]);

%Calcul de la pente 's'
s=1/((1-alpha0)/20);

%Echantillon 1 necessaire pour l'overlapping.
Echantillon = signal_bruite(1:LSample);
F_Echantillon = fftshift(fft(Echantillon));
Argument_E_Signal = angle(F_Echantillon);
PS_Echantillon = abs(F_Echantillon).^2;

%Calcul Alpha
SNR_Trame = 10*log((rms(Echantillon)^2)/(rms(signal_bruite(index_Bruit_SNR*LSample:(index_Bruit_SNR+1)*LSample-1))^2));
if SNR_Trame < -5 
    SNR_Trame = -5;
end
Alpha = alpha0 - SNR_Trame/s;
PS_Echantillon_Debruite = PS_Echantillon - Alpha*PS_Bruit;
    for j=1:LSample
        if PS_Echantillon_Debruite(j) < Beta*PS_Bruit(j)
            PS_Echantillon_Debruite(j) = Beta*PS_Bruit(j); 
        end
    end
Spectre_E_Debruite = sqrt(PS_Echantillon_Debruite);
Echantillon_Debruite = complex(Spectre_E_Debruite.*cos(Argument_E_Signal),Spectre_E_Debruite.*sin(Argument_E_Signal));
Signal_DebruiteArticle=tukeywin(256,0.1).*(real(ifft(ifftshift(Echantillon_Debruite))));

%Autres echantillons
max_index = floor(length(signal)/(LSample_ov+14));
for index=2:max_index
    Echantillon = signal_bruite(((index-1)*LSample_ov)-14+(14*index):(index*LSample_ov)+15+(14*index));  % 14 anciens + 228 nouveaux + 14 futurs = 256
    F_Echantillon = fftshift(fft(Echantillon, FFT_Ordre));
    Argument_E_Signal = angle(F_Echantillon);
    PS_Echantillon = abs(F_Echantillon).^2;
    SNR_Trame = 10*log((rms(Echantillon)^2)/(rms(signal_bruite(index_Bruit_SNR*LSample:(index_Bruit_SNR+1)*LSample-1))^2));
    Alpha = alpha0 - SNR_Trame/s;
    PS_Echantillon_Debruite = G*(PS_Echantillon.^Gamma - Alpha*(PS_Bruit.^Gamma));
    for j=1:LSample
        if PS_Echantillon_Debruite(j)^(1/Gamma) < Beta*PS_Bruit(j)
            PS_Echantillon_Debruite(j) = Beta*PS_Bruit(j);
        else 
            PS_Echantillon_Debruite(j) = PS_Echantillon_Debruite(j)^(1/Gamma);
        end
    end
    Spectre_E_Debruite = sqrt(PS_Echantillon_Debruite);
    Echantillon_Debruite = complex(Spectre_E_Debruite.*cos(Argument_E_Signal),Spectre_E_Debruite.*sin(Argument_E_Signal));
    Signal_Debruite_Tuk=tukeywin(256,0.1).*(real(ifft(ifftshift(Echantillon_Debruite),FFT_Ordre))); %tukey 14 ech de part et d autres
    Signal_DebruiteArticle(length(Signal_DebruiteArticle)-13:length(Signal_DebruiteArticle)) = Signal_Debruite_Tuk(1:14) + Signal_DebruiteArticle(length(Signal_DebruiteArticle)-13:length(Signal_DebruiteArticle));
    Signal_DebruiteArticle=cat(1,Signal_DebruiteArticle, Signal_Debruite_Tuk(15:256)); % /!\ ne pas ajouter les ech de la trame suivante
end
pause(2.5);
sound(real(Signal_DebruiteArticle));

% Affichage des signaux finaux
t2 = (0:length(Signal_Debruite)-1)*(1/8000);
t22 = (0:length(Signal_DebruiteArticle)-1)*(1/8000);
figure;
subplot(411)
plot(t,signal);
title('Signal initial'); 
subplot(412)
plot(t, signal_bruite/max(abs(signal_bruite)));
title('Signal bruité'); 
subplot(413)
plot(t2, real(Signal_Debruite)/max(abs(Signal_Debruite)));
title('Signal débruité avec la solution classique'); 
subplot(414)
plot(t22,Signal_DebruiteArticle/max(abs(Signal_DebruiteArticle)));
title('Signal débruité avec la solution proposée'); xlabel('Temps (s)');

% figure
% plot(t,signal);
% title('Signal initial'); xlabel('Temps (s)');
% hold on
% plot(t22,Signal_DebruiteArticle/max(abs(Signal_DebruiteArticle)));
% hold on
% plot(t2, Signal_Debruite/max(abs(Signal_Debruite)));


%Affichage des signaux finaux en fréquence
f3 = Fs*(-length(Signal_Debruite)/2:length(Signal_Debruite)/2-1)/length(Signal_Debruite);
f33 = Fs*(-length(Signal_DebruiteArticle)/2:length(Signal_DebruiteArticle)/2-1)/length(Signal_DebruiteArticle);
figure
subplot(411)
plot(f2, fftshift(fft(signal)));
title('Signal initial');
subplot(412)
plot(f2, fftshift(fft(signal_bruite)));
title('Signal bruité');
subplot(413)
plot(f3, fftshift(fft(Signal_Debruite)));
title('Signal débruité avec la solution classique');
subplot(414)
plot(f33, fftshift(fft(Signal_DebruiteArticle)));
title('Signal débruité avec la solution proposée');xlabel('Fréquence (en Hz)');


%Calcul des SNR résultants des différentes méthodes 

SNR_base = 10*log((rms(signal_bruite(20*LSample:21*LSample-1))^2)/(rms(signal_bruite(index_Bruit_SNR*LSample:(index_Bruit_SNR+1)*LSample-1))^2))

SNR_classique = 10*log((rms(real(Signal_Debruite(20*LSample:21*LSample-1)))^2)/(rms(real(Signal_Debruite(index_Bruit_SNR*LSample:(index_Bruit_SNR+1)*LSample-1)))^2))

SNR_article = 10*log((rms(Signal_DebruiteArticle(20*LSample:21*LSample-1))^2)/(rms(Signal_DebruiteArticle(index_Bruit_SNR*LSample:(index_Bruit_SNR+1)*LSample-1))^2))



% a=ifft(fft(signal));
% b=ifft(complex(abs(fft(signal)).*cos(angle(fft(signal))),abs(fft(signal)).*sin(angle(fft(signal)))));
% figure;
% subplot(311)
% plot(t, a);
% title('Signal fft and ifft without abs');
% subplot(312)
% plot(t, b);
% title('Signal fft and ifft with abs');
% subplot(313)
% plot(t,a-b);
% title('Erreur lié à abs');