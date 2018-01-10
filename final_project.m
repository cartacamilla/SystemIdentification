clear all
close all
clc

FPdata = load('FPdata.mat');
Te =  0.04;

% Windowing
hann_window = @(M) 0.5+0.5*cos(pi*[0:M-1]'/(M-1));


%% Impulse response using correlation and deconvolution approach
N = length(FPdata.u)/6;
k=400;
u = FPdata.u(end-N+1:end);
y = FPdata.y(end-N+1:end);
Rur = xcorr(u, u, 'unbiased'); % use incor because u periodic
Ryr = xcorr(y, u, 'unbiased');
U = toeplitz(Rur(N:N+k));
theta = pinv(U'*U)*U'*Ryr(N:N+k);

figure
plot(1/Te*theta)
title('impulse response')
hold on

% deconvolution
U = toeplitz(u, [u(1);zeros(N-1,1)]);
Uk = U(:,1:k); % truncate U
theta = pinv(Uk)*y/Te;
plot(theta)
legend('correlation','deconvolution')
hold off

hold on;

window = hann_window(250);
plot(window)
hold off

%% spectral analyis
u = FPdata.u;
y = FPdata.y;
r = FPdata.r;

N_AVG = 1;
N = floor(length(u)/N_AVG);

SCALEOPT = 'biased';
window = ones(N,1);
%window = hann(200);
window = window(length(window)/2 +1:end);
padding = zeros(N - length(window), 1);
window = [window; padding];

fft_Ryr = zeros(N,1);
fft_Rur = zeros(N,1);
for i = 0:N_AVG-1
    % split into chunks
    yc = y(i*N + 1:(i+1)*N);
    uc = u(i*N + 1:(i+1)*N);
    rc = r(i*N + 1:(i+1)*N);

    % correlation
    Ryr = xcorr(yc,rc, SCALEOPT);
    Rur = xcorr(uc,rc, SCALEOPT);

    Ryr = Ryr(N:end);
    Rur = Rur(N:end);

    % windowing
    Ryr = Ryr .* window;

    % averaging
    fft_Ryr = fft_Ryr+fft(Ryr);
    fft_Rur = fft_Rur+fft(Rur);
end

% Reconstruction
Gr = fft_Ryr./fft_Rur;

omega_s = 2*pi/Te;
freq = 0:omega_s/N:(N-1)/N*omega_s;
freq = (0:N-1).*(omega_s/N);

NYQUIST_INDEX = round(N/2);
Gr = Gr(1:NYQUIST_INDEX);
freq = freq(1:NYQUIST_INDEX);

model = frd(Gr, freq, Te);
figure
bode(model)
title('spectral analysis')

%%% fourier analyis
u = FPdata.u;
y = FPdata.y;

N_PERIODS = 6;
N = length(FPdata.u)/N_PERIODS;
omega_s = 2*pi/Te;
avg_y = zeros(N,1);
avg_u = zeros(N,1);
for i = 1:N:N*N_PERIODS
    in = y(i:i+N-1);
    avg_y = avg_y + fft(in);
    out = u(i:i+N-1);
    avg_u = avg_u + fft(out);
end
freq = [];
for i = 0:N-1
   freq = [freq; i*omega_s/N];
end
freq = (0:N-1).*(omega_s/N);

Y = avg_y / N_PERIODS;
U = avg_u / N_PERIODS;

% Reconstruction
Gr = Y ./ U;

NYQUIST_INDEX = round(N/2);
freq = freq(1:NYQUIST_INDEX);
Gr = Gr(1:NYQUIST_INDEX);

model = frd(Gr, freq, Te);
%figure
hold on
bode(model)
%title('fourier analysis')
title('frequency analysis')
legend('spectral analysis', 'fourier analysis')
hold off

% we see 2 resonance frequencies at 14.1 rad/s and 26.6 rad/s
%% Analysis of physical system
syms x1 x2 x3 J1 J2 J3 k1 k2 k3 T s
A=[
    (J1*s^2+k1), -k1, 0;
    -k1, (J2*s^2+k1+k2), -k2;
    0, -k2, (J3*s^2+k2)
];
B=[T;0;0];
x=[x1;x2;x3];
res = solve(A*x-B, x);
G = res.x3/T;

% G has the structure K * 1/s^2 * 1/((s/w1)^2+1) * 1/((s/w2)^2+1);
% J1 = 1; J2 = 1; J3 = 1; k1 = 1; k2 = 1; K3 = 1;
% [n,d] = numden(G);
% n = sym2poly(subs(n));
% d = sym2poly(subs(d));
% G = tf(n,d);
%bode(G)

%%
w1 = 14.1;
w2 = 26.6;
k = 10;
s = tf('s');
G = k * 1/s^2 * 1/(((s/w1)^2+1)*((s/w2)^2+1));
hold on;
bode(G, freq)
hold off;
%% state space model
syms x1 x2 x3 dx1 dx2 dx3 J1 J2 J3 k1 k2 k3 km u
A=[zeros(3), eye(3);
    [
        -k1/J1, k1/J1, 0;
        k1/J2, -(k1+k2)/J2, k2/J2;
        0, k2/J3, -k2/J3
    ], zeros(3)
];
B = [0 0 0 km 0 0]';
C = [0 0 1 0 0 0];
D = zeros(1,1);
x = [x1; x2; x3; dx1; dx2; dx3];
%% spectral analyis without averaging
u = FPdata.u;
y = FPdata.y;
r = FPdata.r;
N = length(u);
omega_s = 2*pi/Te;

SCALEOPT = 'biased';
window = ones(N,1);
window = hann(200);
window = window(length(window)/2 +1:end);
padding = zeros(N - length(window), 1);
window = [window; padding];

% correlation
Ryr = xcorr(y,r, SCALEOPT);
Rur = xcorr(u,r, SCALEOPT);
Ryr = Ryr(N:end);
Rur = Rur(N:end);
% windowing
Ryr = Ryr .* window;
% averaging
fft_Ryr = fft(Ryr);
fft_Rur = fft(Rur);
% Reconstruction
Gr = fft_Ryr./fft_Rur;
freq = (0:N-1).*(omega_s/N);
model = frd(Gr, freq, Te);
figure
bode(model)
title('spectral analysis')
