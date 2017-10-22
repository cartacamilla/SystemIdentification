% input signal
Uprbs = 0.5* prbs(7,4);
Te = 0.2; % sample time
N = length(Uprbs);
sim_time = N*Te;

% simulation
simin = struct();
simin.signals = struct('values', Uprbs);
simin.time = linspace(0,sim_time, N);
sim('ce1_1_sim')

% extract one period of the signals
N = length(Uprbs)/4;
Uprbs = Uprbs(1:N);
Y = simout(end+1-N:end);

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

close all;
hold on;
plot(theta_intcor);
plot(theta_xcorr);
