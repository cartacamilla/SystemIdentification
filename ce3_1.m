laserbeamdata = load('laserbeamdataN.mat');
y = laserbeamdata.y;
u = laserbeamdata.u;
Te = 1e-3; % sampling time

%% Spectral analysis
% we have a random non-periodic input signal u and the output y
SCALEOPT = 'biased';

Ryu = xcorr(y,u, SCALEOPT);
Ruu = xcorr(u,u, SCALEOPT);

N = length(u);
Ryu = Ryu(N:end);
Ruu = Ruu(N:end);

Gr = fft(Ryu)./fft(Ruu);

omega_s = 2*pi/Te;
freq = 0:omega_s/N:(N-1)/N*omega_s;

NYQUIST_INDEX = round(N/2);
freq = freq(1:NYQUIST_INDEX);
Gr = Gr(1:NYQUIST_INDEX);
model = frd(Gr, freq, Te);

% Windowing
WINDOW_SIZE = 40;

hann = @(M) 0.5+0.5*cos(pi*[0:M-1]'/(M-1));
padding = zeros(length(Ryu) - WINDOW_SIZE, 1);
window = [hann(WINDOW_SIZE); padding];
Ryu_hann = Ryu.*window;
Gr = fft(Ryu_hann)./fft(Ruu);
Gr = Gr(1:NYQUIST_INDEX);
model_hann = frd(Gr, freq, Te);

% Bode plot
figure
hold on
bode(model)
bode(model_hann)
legend('whole data','windowing')
hold off
