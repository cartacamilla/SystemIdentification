% input signal
saturation = 0.5;
noiseVariance = 0.1;

N_PERIODS = 1;
PERIOD_LEN = 2^10-1;
u = 0.5*prbs(10,N_PERIODS);
Te = 0.2; % sample time
N = length(u);
sim_time = N*Te;

% simulation
simin = struct();
simin.signals = struct('values', u);
simin.time = linspace(0,N*Te, N);
sim('ce1_1_sim')
y = simout;

omega_s = 2*pi/Te;
freq = [];
for i = 0:PERIOD_LEN-1
   freq = [freq; i*omega_s/127];
end
NYQUIST_INDEX = round(PERIOD_LEN/2);
%% Spectral analysis method
Ryu = intcor(y,u);
Ruu = intcor(u,u);

Gr = fft(Ryu)./fft(Ruu);

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
