clear all
close all
clc

FPdata = load('FPdata.mat');
Te =  0.04;

u = FPdata.u;
y = FPdata.y;

%% Impulse response using correlation approach
%imp_resp_corr_approach(y,u,Te);

%% spectral and fourier analyis
freq_reconstruction_FA = fourier_analysis(y,u,Te,true);

freq_reconstruction_SA = spectral_analysis(y,u,Te,6,true,'hann');

% we see 2 resonance frequencies at 114 rad/s and 216 rad/s

%% parametric identification

order_estimation(y,u,Te,freq_reconstruction_SA);

%% parametric models

parametric_models(y,u,Te,6,4,2);




