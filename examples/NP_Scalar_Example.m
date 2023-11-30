% 
% Neyman-Pearson Detector - Scalar Example
%
% Threshold and Detection probability for NP detector w/ MC Simulation
%
% DETECTION - SiSy - 14/11/2023 - Jules GOMEL ☺
% AY 2023/2024 - Prof. S. Bidon

%% Constants
clear all
close all
clc

% False alarm probability Pfa
Pfa=0.01;

% Number of antennas N, dimension
N=1;
% Noise power
sigma2=1;
% a > 0 : value for alternative hyp H1
a=5;

% Number of Monte-Carlo samples K
K=200/Pfa;

%% Pfa-Threshold

% Generate samples following H0
x_H0=sqrt(sigma2)*randn(N,K);

% Derive test statistics (in our case, t=x)
t_H0=sort(x_H0);

% Compute theoretical threshold
eta_th=norminv(1-Pfa,0,sqrt(sigma2));

% Find exp eta so that Nbr_fa=MC*Pfa
eta_exp=t_H0(floor(K*(1-Pfa)));
%% Detection Probability

% Generate samples following H1
x_H1=a+sqrt(sigma2)*randn(N,K);

% Derive test statistics (in our case, t=x)
t_H1=x_H1;

Pd_th=1-normcdf(eta_th,a,sqrt(sigma2));
Pd_exp=length(find(t_H1>eta_exp))/K;

