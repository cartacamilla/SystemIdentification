close all;
hold off;

% input signal
Uprbs = 0.5* prbs(7,4);

% sample time
Te = 0.2;

N = length(Uprbs);
sim_time = N*Te; % seconds

simin = struct();
simin.signals = struct('values', Uprbs);
simin.time = linspace(0,sim_time, N);

sim('ce1_1_sim')
plot(simout)

Ruu = intcor(Uprbs, Uprbs);
Ruu_x = xcorr(Uprbs, Uprbs);

Y = simout;

Ryu = intcor(Y, Uprbs);
Ryu_x = xcorr(Y, Uprbs)

U = toeplitz(Ruu(1:N,1));

theta = pinv(U)*Ryu(3*N:length(Ryu),1);

plot(theta)