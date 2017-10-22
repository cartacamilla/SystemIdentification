noiseVariance = 0.1;

% input signal
Uprbs = 0.5* prbs(7,4);
Te = 0.2; % sample time
N = length(Uprbs);
sim_time = N*Te;
saturation = 0.5;

% simulation
simin = struct();
simin.signals = struct('values', Uprbs);
simin.time = linspace(0,sim_time, N);
sim('ce1_1_sim')

% extract one period of the signals
N = length(Uprbs)/4;
Uprbs = Uprbs(1:N);
Y = simout(end+1-N:end)/Te;
time = simin.time(1:N);

% Correlation approach using intcor
Ruu = intcor(Uprbs, Uprbs);
Ryu = intcor(Y, Uprbs);
U = toeplitz(Ruu);
theta_intcor = pinv(U)*Ryu;

% Correlation approach using xcorr
Ryu = xcorr(Y, Uprbs);
Ruu = xcorr(Uprbs, Uprbs);
U = toeplitz(Ruu);
theta_xcorr = pinv(U)*Ryu;

G = tf([4],[1 1 4], 'InputDelay', Te);
Z = c2d(G, Te, 'zoh');
[G_impulse, G_time] = impulse(Z, time(end));

% error
err_intcor = norm(theta_intcor(1:127) - G_impulse)
err_xcorr = norm(theta_xcorr(1:127) - G_impulse)

close all;
hold on;
stairs(time, theta_intcor(1:N));
stairs(time, theta_xcorr(1:N));
plot(G_time, G_impulse);

