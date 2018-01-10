estim_data = detrend(iddata(y(1:N/2),u(1:N/2),Te));
valid_data = detrend(iddata(y(N/2+1:end),u(N/2+1:end),Te));

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

