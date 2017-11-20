function model = spectral_analysis_avg(y,u,Te,N_AVG,SCALEOPT,window)

N = floor(length(u)/N_AVG);

if nargin < 5
    SCALEOPT = 'biased';
end
if nargin < 6
    window = ones(N,1);
end
padding = zeros(N - length(window), 1);
window = [window; padding];

fft_Ryu = zeros(N,1);
fft_Ruu = zeros(N,1);
for i = 0:N_AVG-1
    % split into chunks
    yc = y(i*N + 1:(i+1)*N);
    uc = u(i*N + 1:(i+1)*N);

    % correlation
    Ryu = xcorr(yc,uc, SCALEOPT);
    Ruu = xcorr(uc,uc, SCALEOPT);

    Ryu = Ryu(N:end);
    Ruu = Ruu(N:end);

    % windowing
    Ryu = Ryu .* window;

    % averaging
    fft_Ryu = fft_Ryu+fft(Ryu);
    fft_Ruu = fft_Ruu+fft(Ruu);
end

% Reconstruction
Gr = fft_Ryu./fft_Ruu;

omega_s = 2*pi/Te;
freq = 0:omega_s/N:(N-1)/N*omega_s;

NYQUIST_INDEX = round(N/2);
Gr = Gr(1:NYQUIST_INDEX);
freq = freq(1:NYQUIST_INDEX);

model = frd(Gr, freq, Te);
