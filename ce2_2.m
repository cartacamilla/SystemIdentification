SCALEOPT = 'biased';

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
title('Bode Diagram: whole data (no window)')
hold off
printpdf(gcf, 'ce2_2_no_window.pdf');

figure
hold on
bode(model_hann)
bode(Z,model_hann.Frequency)
title('Bode Diagram: Hann window')
hold off
printpdf(gcf, 'ce2_2_hann_window.pdf');

%% Averaging
N_AVG = 10;

% averaging without window
model = spectral_analysis_avg(y,u,Te,N_AVG,'biased');

% averaging with Hann window
window = hann(40);
model_hann = spectral_analysis_avg(y,u,Te,N_AVG,'biased',window);

figure
hold on
bode(model)
bode(Z,model.Frequency)
title('Bode Diagram: Averaging')
hold off

figure
hold on
bode(model_hann)
bode(Z,model_hann.Frequency)
title('Bode Diagram: Averaging with Hann window')
hold off

printpdf(gcf, 'ce2_2_averaging.pdf');

