function model = spectral_analysis(y,u,Te,SCALEOPT,window)

N = length(u);

if nargin < 4
    SCALEOPT = 'biased';
end
if nargin < 5
    window = ones(N,1);
end

% correlation
Ryu = xcorr(y,u, SCALEOPT);
Ruu = xcorr(u,u, SCALEOPT);

Ryu = Ryu(N:end);
Ruu = Ruu(N:end);

% Windowing
padding = zeros(N - length(window), 1);
window = [window; padding];
Ryu = Ryu.*window;

% Reconstruction
Gr = fft(Ryu)./fft(Ruu);

omega_s = 2*pi/Te;
freq = 0:omega_s/N:(N-1)/N*omega_s;

NYQUIST_INDEX = round(N/2);
Gr = Gr(1:NYQUIST_INDEX);
freq = freq(1:NYQUIST_INDEX);

model = frd(Gr, freq, Te);
