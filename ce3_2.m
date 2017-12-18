clear all
clc
close all

laserbeamdata = load('laserbeamdataN.mat');
y = laserbeamdata.y;
u = laserbeamdata.u;
Te = 1e-3; % sampling time
%%
N = length(y);
m = 50;
Phi = toeplitz(u, [u(1);zeros(N-1,1)]);

% truncate
Phi = Phi(:,1:m);

% least squares
% theta = pinv(Phi'*Phi)*Phi'*y;
theta = Phi\y
figure
stem(theta)
xlabel('m')
ylabel('Value')

%%
% reconstruct
yhat = Phi*theta;

x_ax = Te*[1:length(yhat)];
% compare
figure
plot(x_ax,y);
hold on
plot(x_ax,yhat);
hold off

xlabel('Time[s]')
ylabel('Amplitude')
legend('Real output','Reconstruction')

% loss function
J = sum((y - yhat).^2)/N

%%

% covariance of parameters theta

var_hat = J/(N-m);

cov = var_hat*pinv(Phi'*Phi);
two_sigma = 2*sqrt(diag(cov));

figure
errorbar(theta, two_sigma)

xlabel('Time')
ylabel('Value')
