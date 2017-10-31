% input signal
saturation = 0.5;
noiseVariance = 0.1;
%noiseVariance = 0;

N_PERIODS = 1;
u = 0.5*prbs(10,N_PERIODS);
PERIOD_LEN = length(u)/N_PERIODS;
Te = 0.2; % sample time
N = length(u);
sim_time = N*Te;

% simulation
simin = struct();
simin.signals = struct('values', u);
simin.time = linspace(0,N*Te, N);
sim('ce1_1_sim')
y = simout;


%% Spectral analysis method
Ryu = xcorr(y,u);
Ruu = xcorr(u,u);

Ryu = Ryu(N:end);
Ruu = Ruu(N:end);

Gr = fft(Ryu)./fft(Ruu);

omega_s = 2*pi/Te;
freq = 0:omega_s/N:(N-1)/N*omega_s;

NYQUIST_INDEX = round(PERIOD_LEN/2);
freq = freq(1:NYQUIST_INDEX);
Gr = Gr(1:NYQUIST_INDEX);
model = frd(Gr, freq, Te);
bode(model)

% compare with true system
hold on
G = tf([4],[1 1 4]);
Z = c2d(G, Te, 'zoh');
bode(Z,freq)
hold off
