
% Windowing
hann_window = @(M) 0.5+0.5*cos(pi*[0:M-1]'/(M-1));


%% Impulse response using correlation and deconvolution approach
N = length(FPdata.u)/6;
k=400; % max length
u = FPdata.u(end-N+1:end);
y = FPdata.y(end-N+1:end);
Ruu = intcor(u, u); % use intcor because u periodic
Ryu = intcor(y, u);
U = toeplitz(Ruu(N:N+k));
theta = pinv(U'*U)*U'*Ryu(N:N+k);

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

%% spectral analyis
u = FPdata.u;
y = FPdata.y;
r = FPdata.r;

N_AVG = 1;
N = floor(length(u)/N_AVG);

SCALEOPT = 'biased';
%window = ones(N,1);
window = hann_window(250);
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
    Rur = Rur .* window;
    % averaging
    fft_Ryr = fft_Ryr+fft(Ryr);
    fft_Rur = fft_Rur+fft(Rur);
end

% Reconstruction
Gr = fft_Ryr./fft_Rur;

omega_s = 2*pi/Te;
freq = (0:N-1).*(omega_s/N);

NYQUIST_INDEX = round(N/2);
Gr = Gr(1:NYQUIST_INDEX);
freq = freq(1:NYQUIST_INDEX);

spectral_analysis_model = frd(Gr, freq, Te);
figure
hold on;
bode(spectral_analysis_model); hold off;
title('spectral analysis')

%% fourier analyis
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

fourier_analysis_model = frd(Gr, freq, Te);
%figure
hold on
bode(fourier_analysis_model)
%title('fourier analysis')
title('frequency analysis')
legend('spectral analysis', 'fourier analysis')
hold off

% we see 2 resonance frequencies at 14.1 rad/s and 26.6 rad/s
