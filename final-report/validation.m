u = FPdata.u;
y = FPdata.y;
N = length(y);
train_data = detrend(iddata(y(1:N/2),u(1:N/2),Te));
valid_data = detrend(iddata(y(N/2+1:end),u(N/2+1:end),Te));

n = max(na_est, nb_est)
na = n;
nb = n;
nk = nk_est;
nc = n; nd = n; nf = n; nx = n;

models = struct('sys',{},'name',{});
models(1).sys = arx(train_data, [na nb nk]);
models(1).name = 'ARX';
models(2).sys = iv4(train_data, [na nb nk]);
models(2).name = 'Instrumental Variables';
models(3).sys = armax(train_data, [na nb nc nk]);
models(3).name = 'ARMAX';
models(4).sys = oe(train_data, [nb nf nk]);
models(4).name = 'Output-Error';
models(5).sys = bj(train_data, [nb nc nd nf nk]);
models(5).name = 'Box-Jenkins';
models(6).sys = n4sid(train_data, nx);
models(6).name = 'State-space';

%% Time domain validation
figure
for i = 1:length(models)
    subplot(3,2,i)
    compare(valid_data, models(i).sys)
    title(models(i).name)
end
printpdf(gcf, 'final-report/images/4_time_domain_valid.pdf', 1, 1.8)
%% time domain compare
fitval = [];
names = {};
for i = 1:length(models)
    [YT, FIT, X0] = compare(valid_data, models(i).sys);
    fitval = [fitval; FIT];
    names{i} = models(i).name;
end
figure
bar(fitval)
ylim([0 Inf])
xticklabels(names);
xtickangle(45)
ylabel('Fit (%)')
title('Time-domain comparison')
saveas(gcf, 'final-report/images/4_time_domain_comp', 'png');
%% Statistical validation
figure
for i = 1:3
    subplot(3,1,i)
    resid(valid_data, models(i).sys)
    title(models(i).name)
end
printpdf(gcf, 'final-report/images/4_resid_1.pdf', 1, 1.8)
figure
for i = 4:6
    subplot(3,1,i-3)
    resid(valid_data, models(i).sys)
    title(models(i).name)
end
printpdf(gcf, 'final-report/images/4_resid_2.pdf', 1, 1.8)
%% Frequency-domain validation
figure
for i = 1:length(models)
    subplot(3,2,i)
    compare(spectral_analysis_model, models(i).sys)
    title(models(i).name)
end
printpdf(gcf, 'final-report/images/4_freq_domain_valid.pdf', 1.2, 2.4)
%% Frequency-domain compare
fitval = [];
names = {};
for i = 1:length(models)
    [YT, FIT, X0] = compare(spectral_analysis_model, models(i).sys);
    fitval = [fitval; FIT];
    names{i} = models(i).name;
end
figure
bar(fitval)
ylim([0 Inf])
xticklabels(names);
xtickangle(45)
ylabel('Fit (%)')
title('Frequency-domain comparison')
saveas(gcf, 'final-report/images/4_freq_domain_comp', 'png');
%% Frequency-domain manual comparsion
% define valid frequency spectrum, cut out noisy part where phase explodes.
freq = spectral_analysis_model.Frequency(2:1450);
figure
for i = 1:6
    subplot(3,2,i)
    bode(models(i).sys); hold on;
    bode(spectral_analysis_model, freq); hold off;
    %legend(models(i).name, 'SA');
    title(models(i).name)
end
printpdf(gcf, 'final-report/images/4_freq_visual_comp.pdf', 1, 1.8)