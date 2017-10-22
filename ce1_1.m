close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.1

Te = 0.2;
maxTime = 20;

%Shannon
w0=2;
T = 1/(2*w0);
sim_time = 50;
%The sampling period is sufficiently small
sample_length = sim_time/Te;
saturation = 0.5;

G = tf([4],[1,1,4], 'InputDelay',1);
G_discrete = c2d(G, Te, 'zoh');
[G_step, G_time] = step(G_discrete, sim_time-Te);

[G_impulse, G_time] = impulse(G_discrete, sim_time-Te);


input = [zeros(1/Te,1);ones(sample_length - 1/Te,1)];

simin = struct();
simin.signals = struct('values', input);
simin.time = linspace(0,sample_length*Te, sample_length);

%The noise is masking the signal but we can still recognize a second order
%function

noiseVariance = 0;
sim('ce1_1_sim')

figure
stairs(simin.time(1:maxTime/Te), simout(1:maxTime/Te)/saturation); hold on;

noiseVariance = 0.01;
sim('ce1_1_sim')

stairs(simin.time(1:maxTime/Te), simout(1:maxTime/Te)/saturation);hold on; 

plot(G_time(1:maxTime/Te), G_step(1:maxTime/Te)); 



title('Step Response','Interpreter','latex')
legend('Simulated step response without noise','Simulated step response with noise', 'True step response')
xlabel('Time(s)','Interpreter','latex')
ylabel('Amplitude','Interpreter','latex')


input = [zeros(1/Te,1);1;zeros(sample_length - 1/Te-1,1)];

simin = struct();
simin.signals = struct('values', input);
simin.time = linspace(0,sample_length*Te, sample_length);

%The noise is masking the signal but we can still recognize a second order
%function

noiseVariance = 0;
sim('ce1_1_sim')

figure
stairs(simin.time(1:maxTime/Te), simout(1:maxTime/Te)/(saturation*Te)); hold on;

noiseVariance = 0.01;
sim('ce1_1_sim')

stairs(simin.time(1:maxTime/Te), simout(1:maxTime/Te)/(saturation*Te));hold on; 

plot(G_time(1:maxTime/Te), G_impulse(1:maxTime/Te)); 



title('Impulse Response','Interpreter','latex')
legend('Simulated impulse response without noise','Simulated impulse response with noise', 'True impulse response')
xlabel('Time(s)','Interpreter','latex')
ylabel('Amplitude','Interpreter','latex')

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.2
close all;
clc;

x = prbs(5,4);
[R,h] = intcor(x,x);

figure 
grid on
stem(h,R, 'o')

title('Autocorrelation for prbs(5,4) with intcor()','Interpreter','latex')
legend('PRBS with n = 5 and p = 4')

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.3

close all;
clc;

noiseVariance = 0.01;
a = -0.3;
b = -a;
sim_time = 100; % seconds
N = sim_time/Te;

u = a + (b-a).*rand(N,1);

simin = struct();
simin.signals = struct('values', u);
simin.time = linspace(0,N*Te, N);

sim('ce1_1_sim')

U = toeplitz(u, [u(1);zeros(N-1,1)]);
Y = simout;

k = 300;
Uk = U(:,1:k);

theta = pinv(Uk)*Y;
theta = theta;

%theta = U\Y;

figure
stairs(simin.time(1:maxTime/Te), theta(1:maxTime/Te));hold on;

G = tf([4],[1,1,4]);
G_discrete = c2d(G, Te, 'zoh');
[G_impulse, G_time] = impulse(Te*G_discrete, sim_time-Te);
plot(G_time(1:maxTime/Te), G_impulse(1:maxTime/Te)); 


title('Impulse Response using Numerical Deconvolution','Interpreter','latex')
legend('Estimated impulse response','True impulse response')
xlabel('Time(s)','Interpreter','latex')
ylabel('Amplitude','Interpreter','latex')
%%

