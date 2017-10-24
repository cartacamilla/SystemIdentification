% input signal
saturation = 0.5;
noiseVariance = 0.1;

N_PERIODS = 16;
PERIOD_LEN = 2^7-1;
u = 0.5*prbs(7,N_PERIODS);
Te = 0.2; % sample time
N = length(u);
sim_time = N*Te;

% simulation
simin = struct();
simin.signals = struct('values', u);
simin.time = linspace(0,N*Te, N);
sim('ce1_1_sim')

%plot(simout)
%% Fourier transform
omega_s = 2*pi/Te;
avg = zeros(PERIOD_LEN,1);
for i = 1:PERIOD_LEN:N
    sig = simout(i:i+127-1);
    avg = avg + fft(sig);
end
freq = [];
for i = 0:PERIOD_LEN-1
   freq = [freq; i*omega_s/127];
end

NYQUIST_INDEX = round(PERIOD_LEN/2);

Y = avg / N_PERIODS;
plot(freq(1:NYQUIST_INDEX), abs(Y(1:NYQUIST_INDEX)))

figure
U = fft(u(1:127));
plot(freq(1:NYQUIST_INDEX), abs(U(1:NYQUIST_INDEX)))
%% Reconstruction

Gr = Y ./ U;

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
