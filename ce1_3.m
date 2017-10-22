noiseVariance = 0.1;

Te = 0.2;
sim_time = 70;
N = sim_time/Te;
saturation = 0.5;

% random input signal
a = -saturation; % scaling to avoid saturation
b = saturation;
u = a + (b-a).*rand(N,1);

% simulation
simin = struct();
simin.signals = struct('values', u);
simin.time = linspace(0,N*Te, N);
sim('ce1_1_sim')

% deconvolution
U = toeplitz(u, [u(1);zeros(N-1,1)]);
k = 60;
Uk = U(:,1:k); % truncate U

theta = pinv(Uk)*simout/Te;

G = tf([4],[1 1 4]);
H = c2d(G, Te, 'zoh');

% error
[H_impulse, H_time] = impulse(H, simin.time(k));
error = norm(H_impulse - theta)

close all;
hold on;
plot(simin.time(1:k), theta);
plot(H_time, H_impulse);
