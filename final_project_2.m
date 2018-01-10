clear all
close all
clc
%%%
FPdata = load('FPdata.mat');
Te =  0.04;

u = FPdata.u;
y = FPdata.y;
r = FPdata.r;
N = length(y);

data = detrend(iddata(y,u,Te));
estim_data = detrend(iddata(y(1:N/2),u(1:N/2),Te));
valid_data = detrend(iddata(y(N/2+1:end),u(N/2+1:end),Te));


% ARX models for different orders n
%% na estimate
order = 2:10;
loss_func = [];
for na = order
    SYS = arx(data, [na na 0]);
    loss_func = [loss_func, SYS.EstimationInfo.LossFcn];
end

scaling_factor = 0.001;
penalty_BIC = scaling_factor*order*log(N)/N;
global_crit = loss_func + penalty_BIC;
[~,i] = min(global_crit);
na_est = order(i);

figure; hold on;
plot(order, loss_func, 'b');
plot(order, penalty_BIC, '--k')
plot(order, global_crit, '-*r')
hold off;
title('ARX loss function');
xlabel('Model order na');
ylabel('Loss function');
legend('Loss Function', 'Penalty term BIC', 'Global criterion');
na = na_est

%% nk estimate
nb = na; nk = 0;
SYS = arx(data, [na nb nk]);
figure;
% plot B with a confidence interval of 2 sigma
errorbar(SYS.b, SYS.db*2)

lower = SYS.b + 2 * SYS.db;
higher = SYS.b - 2 * SYS.db;
nk_est = 0;
for i = 1:length(SYS.b)
    if(lower(i).*higher(i) <= 0)
        nk_est = nk_est + 1;
    else
        break; % if 0 in 2sigma confidence interval
    end
end
nk = nk_est

%% nb estimate
order = 1:10;
loss_func = [];
for nb = order
    orders = [na nb nk];
    SYS = arx(data, orders);
    loss_func = [loss_func, SYS.EstimationInfo.LossFcn];
end

scaling_factor = 1e-4;
penalty_BIC = scaling_factor*order*log(N)/N;
global_crit = loss_func + penalty_BIC;
[~,i] = min(global_crit);
nb_est = order(i);

figure; hold on;
plot(order, loss_func, 'b');
plot(order, penalty_BIC, '--k')
plot(order, global_crit, '-*r');
hold off;
title('ARX loss function');
xlabel('Model order nb'); ylabel('Loss function');
legend('Loss Function', 'Penalty term BIC', 'Global criterion');
nb = nb_est
n = max(na, nb+nk)

%% Validate order with ARMAX model by observing pole-zero-cancellation
for n=4:8
    SYS = armax(data, [n n n 0]);
    figure
    h = iopzplot(SYS);
    showConfidence(h, 2);
    title('ARMAX n = '+string(n))
end

%% compare with selstruc
%NN = struc(5:7,5:7,0:2);
%V = arxstruc(estim_data, valid_data, NN);
%order = selstruc(V)
%%
n = max(na_est, nb_est + nk_est)
na = n; nb = n; nk = nk_est;
nc = n; nd = n; nf = n; nx = n;

models = struct('sys',{},'name',{});
models(1).sys = arx(estim_data, [na nb nk]);
models(1).name = 'ARX';
models(2).sys = iv4(estim_data, [na nb nk]);
models(2).name = 'Instrumental Variables';
models(3).sys = armax(estim_data, [na nb nc nk]);
models(3).name = 'ARMAX';
models(4).sys = oe(estim_data, [nb nf nk]);
models(4).name = 'Output-Error';
models(5).sys = bj(estim_data, [nb nc nd nf nk]);
models(5).name = 'Box-Jenkins';
models(6).sys = n4sid(estim_data, nx);
models(6).name = 'State-space';

%% time domain validation
fitval = [];
names = {};
for i = 1:length(models)
    [YT, FIT, X0] = compare(valid_data, models(i).sys);
    fitval = [fitval; FIT];
    names{i} = models(i).name;
end
figure
bar(fitval)
xticklabels(names);
xtickangle(45)
%% Time domain plot
figure
for i = 1:length(models)
    subplot(3,2,i)
    compare(valid_data, models(i).sys)
    title(models(i).name)
end
%% Statistical validation
figure
for i = 1:3
    subplot(3,1,i)
    resid(valid_data, models(i).sys)
    title(models(i).name)
end
figure
for i = 4:6
    subplot(3,1,i-3)
    resid(valid_data, models(i).sys)
    title(models(i).name)
end
%% Frequency domain validation

