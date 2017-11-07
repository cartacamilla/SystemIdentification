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

% Windowing
WINDOW_SIZE = 40;

hann = @(M) 0.5+0.5*cos(pi*[0:M-1]'/(M-1));
hamming = @(M) 0.54+0.46*cos(pi*[0:M-1]'/(M-1));

padding = zeros(length(Ryu) - WINDOW_SIZE, 1);
window = [hamming(WINDOW_SIZE); padding];
Ryu_hamming = Ryu.*window;

Gr = fft(Ryu_hamming)./fft(Ruu);
Gr = Gr(1:NYQUIST_INDEX);
model_hamming = frd(Gr, freq, Te);

window = [hann(WINDOW_SIZE); padding];
Ryu_hann = Ryu.*[hann(WINDOW_SIZE); padding];
Gr = fft(Ryu_hamming)./fft(Ruu);
Gr = Gr(1:NYQUIST_INDEX);
model_hann = frd(Gr, freq, Te);

bode(model)

% compare with true system
hold on
bode(model)
bode(model_hann)
bode(model_hamming)
G = tf([4],[1 1 4]);
Z = c2d(G, Te, 'zoh');
bode(Z,freq)
%legend('no window','Hann','Hamming','System')
hold off

%% averaging
M_DIV = 10;
CHUNK_LEN = floor(PERIOD_LEN/M_DIV);

WINDOW_SIZE = 40;
padding = zeros(CHUNK_LEN - WINDOW_SIZE, 1);
window = [hamming(WINDOW_SIZE); padding];

fft_Ryu = zeros(CHUNK_LEN,1);
fft_Ruu = zeros(CHUNK_LEN,1);
for i = 0:9
    yc = y(i*CHUNK_LEN + 1:(i+1)*CHUNK_LEN);
    uc = u(i*CHUNK_LEN + 1:(i+1)*CHUNK_LEN);
    Ryu = xcorr(yc,uc);
    Ruu = xcorr(uc,uc);
    
    Ryu = Ryu(CHUNK_LEN:end);
    Ruu = Ruu(CHUNK_LEN:end);
    
    % windowing
    %Ryu = Ryu .* window;
    
    fft_Ryu = fft_Ryu+fft(Ryu);
    fft_Ruu = fft_Ruu+fft(Ruu);
end

NYQUIST_INDEX = round(CHUNK_LEN/2);
Gr = fft_Ryu./fft_Ruu;
Gr = Gr(1:NYQUIST_INDEX);

omega_s = 2*pi/Te;
freq = 0:omega_s/CHUNK_LEN:(CHUNK_LEN-1)/CHUNK_LEN*omega_s;
freq = freq(1:NYQUIST_INDEX);

model = frd(Gr, freq, Te);

figure
hold on
bode(model)

G = tf([4],[1 1 4]);
Z = c2d(G, Te, 'zoh');
bode(Z,freq)
hold off