clear all
close all
clc

FPdata = load('FPdata.mat');
Te =  0.04;

%%
run('frequency_analysis.m')

%%
run('physical_model.m')

%%
run('estimation.m')

%%
run('validation.m')
