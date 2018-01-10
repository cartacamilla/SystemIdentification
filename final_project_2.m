clear all
close all
clc

FPdata = load('FPdata.mat');
Te =  0.04;

u = FPdata.u;
y = FPdata.y;
r = FPdata.r;
N = length(y);

data = detrend(iddata(y,u,Te));
estim_data = detrend(iddata(y(1:N/2),u(1:N/2),Te));
valid_data = detrend(iddata(y(N/2+1:end),u(N/2+1:end),Te));

%% parametric identification

% ARX models for different orders n, plot loss function
% -> model_arx = arx(data, [n n 0])
% -> model_arx.EstimationInfo.LossFcn

LossFcn = [];
for n = 1:10
    model_arx = arx(data, [n n 0]);
    LossFcn = [LossFcn; model_arx.EstimationInfo.LossFcn];
end
bar(1:10,LossFcn)
title('ARX loss function')
xlabel('Order')
ylabel('EstimationInfo.LossFcn')
% We observe that after order 3 or 4 the loss function becomes flat
n = 4;

model_arx = arx(data, [n n 0]);
figure;
% plot B with a confidence interval of 2 sigma
errorbar(model_arx.b, model_arx.db*2)
% we observe none of the first coefficients of B are zero -> delay = 0
nk = 0;

% Validate order with ARMAX model
% plot pole-zero-cancellation iopzplot(), showConfidence(h, 2)
for n=1:10
    model_armax = armax(data, [n n n 0]);
    figure
    h = iopzplot(model_armax);
    showConfidence(h, 2);
    title('ARMAX n = '+string(n))
end
% for orders above 6 zeros and poles start cancelling each other out
n = 6;

figure
m = armax(data, [n,n,n,0]);
errorbar(m.b, 2*m.db)
% confirms delay = 0
nk = 0;

%% compare with selstruc
NN = struc(1:10,1:10,0:10);
V = arxstruc(estim_data, valid_data, NN);
order = selstruc(V)
%%
na = 4;
nb = 1;
nk = 4;
nc = na;
nd = na;
nf = na;
nx = max(nb + nk, na);


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

fitval = [];
names = {};
for i = 1:length(models)
    [YT, FIT, X0] = compare(valid_data, models(i).sys);
    fitval = [fitval; FIT];
    names{i} = models(i).name;
end
bar(fitval)
xticklabels(names);
xtickangle(45)

%%
figure
title('Validation')
for i = 1:3
    subplot(3,1,i)
    resid(valid_data, models(i).sys)
    title(models(i).name)
end
figure
title('Validation')
for i = 4:6
    subplot(3,1,i)
    resid(valid_data, models(i).sys)
    title(models(i).name)
end
