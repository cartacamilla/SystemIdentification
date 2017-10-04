Te = 0.2;

%Shannon
w0=2;
T = 1/(2*w0);

%The sampling period is sufficiently small

sample_length = 50/Te;


input = 0.5*[zeros(1/Te,1);ones(sample_length - 1/Te,1)];

simin = struct();
simin.signals = struct('values', input);
simin.time = linspace(0,sample_length*Te, sample_length);

%sim('ce1_1_sim')

%The noise is masking the signal but we can still recognize a second order
%function

simin.signals.values = [1; zeros(sample_length-1,1)];

sim('ce1_1_sim')

plot(simout)