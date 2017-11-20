% simulation parameters
saturation = 0.5;
noiseVariance = 0.1;

% input signal
u = 0.5*prbs(10,1);
PERIOD_LEN = length(u);
Te = 0.2; % sample time
N = length(u);
sim_time = N*Te;

% simulation
simin = struct();
simin.signals = struct('values', u);
simin.time = linspace(0,N*Te, N);
sim('ce1_1_sim')
y = simout;

% true system
G = tf([4],[1 1 4]);
Z = c2d(G, Te, 'zoh');


%% Spectral analysis method
model = spectral_analysis(y,u,Te,'biased');

% Windowing
hann = @(M) 0.5+0.5*cos(pi*[0:M-1]'/(M-1));
hamming = @(M) 0.54+0.46*cos(pi*[0:M-1]'/(M-1));

window = hann(40);
model_hann = spectral_analysis(y,u,Te,'biased',window);

% Bode plot
figure
hold on
bode(model)
bode(Z,model.Frequency)
title('Spectral analysis: whole data (no window)')
hold off
printpdf(gcf, 'ce2_2_no_window.pdf');

figure
hold on
bode(model_hann)
bode(Z,model_hann.Frequency)
title('Spectral analysis: Hann window')
hold off
printpdf(gcf, 'ce2_2_hann_window.pdf');

%% Averaging
N_AVG = 8;

% averaging without window
model = spectral_analysis_avg(y,u,Te,N_AVG,'biased');

% averaging with Hann window
window = hann(40);
model_hann = spectral_analysis_avg(y,u,Te,N_AVG,'biased',window);

figure
subplot(1,2,1)
hold on
bode(model)
bode(Z,model.Frequency)
title('Spectral analysis: Averaging')
hold off
subplot(1,2,2)
hold on
bode(model_hann)
bode(Z,model.Frequency)
title('Spectral analysis: Averaging with Hann window')
hold off

printpdf(gcf, 'ce2_2_averaging.pdf',1.5,1);

%% unbiased plots
figure
subplot(2,2,1)
hold on
model_biased = spectral_analysis(y,u,Te,'biased');
model_unbiased = spectral_analysis(y,u,Te,'unbiased');
bode(model_biased)
bode(model_unbiased)
bode(Z,model_biased.Frequency)
title('Spectral analysis: whole data')
legend('biased','unbiased','true model')
hold off

window = hann(40);
subplot(2,2,2)
hold on
model_biased = spectral_analysis(y,u,Te,'biased',window);
model_unbiased = spectral_analysis(y,u,Te,'unbiased',window);
bode(model_biased)
bode(model_unbiased)
bode(Z,model_biased.Frequency)
title('Spectral analysis: Windowing')
legend('biased','unbiased','true model')
hold off

subplot(2,2,3)
hold on
model_biased = spectral_analysis_avg(y,u,Te,N_AVG,'biased');
model_unbiased = spectral_analysis_avg(y,u,Te,N_AVG,'unbiased');
bode(model_biased)
bode(model_unbiased)
bode(Z,model_biased.Frequency)
title('Spectral analysis: Averaging')
legend('biased','unbiased','true model')
hold off

printpdf(gcf, 'ce2_2_unbiased.pdf',1.5,2);

