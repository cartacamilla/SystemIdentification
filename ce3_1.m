laserbeamdata = load('laserbeamdataN.mat');
y = laserbeamdata.y;
u = laserbeamdata.u;
Te = 1e-3; % sampling time

%% Spectral analysis
% we have a random non-periodic input signal u and the output y
SCALEOPT = 'biased';

model = spectral_analysis(y, u, Te, SCALEOPT);

% Windowing
WINDOW_SIZE = 100;
hann = @(M) 0.5+0.5*cos(pi*[0:M-1]'/(M-1));
window = hann(WINDOW_SIZE);

model_hann = spectral_analysis(y, u, Te, SCALEOPT, window);

N_AVG = 4;
model_avg = spectral_analysis_avg(y,u,Te,N_AVG,SCALEOPT);
model_avg_hann = spectral_analysis_avg(y,u,Te,N_AVG,SCALEOPT, window);

% Bode plot
figure
hold on
bode(model)
%bode(model_avg)
bode(model_hann, model_hann.Frequency)
hold off
