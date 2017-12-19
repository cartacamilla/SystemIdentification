clear all
close all
clc

FPdata = load('FPdata.mat');
Te =  0.04;

%% Impulse response using correlation and deconvolution approach
N = length(FPdata.u)/6;
k=250;
u = FPdata.u(end-N+1:end);
y = FPdata.y(end-N+1:end);
Ruu = intcor(u, u); % use incor because u periodic
Ryu = intcor(y, u);
U = toeplitz(Ruu(N:N+k));
theta = pinv(U'*U)*U'*Ryu(N:N+k);

figure
plot(1/Te*theta)
hold on

% deconvolution
U = toeplitz(u, [u(1);zeros(N-1,1)]);
Uk = U(:,1:k); % truncate U
theta = pinv(Uk)*y/Te;
plot(theta)
hold off

%% Initialisation
u = FPdata.u;
y = FPdata.y;
N = length(y);
data = iddata(y,u,Te);

%% spectral analysis
model = spectral_analysis(y,u,Te,'biased',hamming(10));
bode(model)

%% fourier analyis
N_PERIODS = 6;
PERIOD_LEN = N/N_PERIODS;
omega_s = 2*pi/Te;
avg = zeros(PERIOD_LEN,1);
for i = 1:PERIOD_LEN:N
    sig = y(i:i+PERIOD_LEN-1);
    avg = avg + fft(sig);
end
freq = [];
for i = 0:PERIOD_LEN-1
   freq = [freq; i*omega_s/127];
end

Y = avg / N_PERIODS;
U = fft(u(1:PERIOD_LEN));

% Reconstruction
Gr = Y ./ U;

NYQUIST_INDEX = round(PERIOD_LEN/2);
freq = freq(1:NYQUIST_INDEX);
Gr = Gr(1:NYQUIST_INDEX);

model = frd(Gr, freq, Te);
hold on
bode(model)
