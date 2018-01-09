clear all
close all
clc

FPdata = load('FPdata.mat');
Te =  0.04;

%% parametric identification

u = FPdata.u;
y = FPdata.y;
r = FPdata.r;
N = length(y);


data = detrend(iddata(y,u,Te));
testing = detrend(iddata(y(1:N/2),u(1:N/2),Te));
validation = detrend(iddata(y(N/2+1:end),u(N/2+1:end),Te));


% ARX models for different orders n, plot loss function
% -> model_arx = arx(data, [n n 0])
% -> model_arx.EstimationInfo.LossFcn
% ARMAX models, plot pole-zero-cancellation iopzplot(), showConfidence(h, 2)

% find delay nk or d respectively
% errorbar(model_arx.b, model_arx.db, 'x')

% NN = struc(1:10, 1:10, 1:10);
% V = arxstruc(testing, validation, NN)
% order = selstruc(V)

% na=, nb=, nk=
% model_arx = arx(testing, [na nb nk]) armax() oe() bj() n4sid()
% [YT, FIT, X0] = compare(validation, model_arx)
% fit_compare = [fit_compare; FIT]
% models = categorical({'ARX','ARMAX','OE','BJ',SS});
% prices = [1.23 0.99 2.3];
% bar(models, fit_compare)

% resid(validation, )