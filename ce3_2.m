laserbeamdata = load('laserbeamdataN.mat');
y = laserbeamdata.y;
u = laserbeamdata.u;
Te = 1e-3; % sampling time

N = length(y);
m = 50;
Phi = toeplitz(u, [u(1);zeros(N-1,1)]);

% truncate
Phi = Phi(:,1:m);

% least squares
theta = pinv(Phi'*Phi)*Phi'*y;

% reconstruct
yhat = Phi*theta;

% compare
plot(y);
hold on
plot(yhat);
hold off

% loss function
J = sum((y - yhat).^2);


% covariance of parameters theta

var_hat = J/(N-m);

cov = var_hat*pinv(Phi'*Phi);
two_sigma = 2*sqrt(diag(cov));

errorbar(theta, two_sigma)
