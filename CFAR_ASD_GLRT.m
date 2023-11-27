% 
% CFAR Adaptive Subspace Detector - GLRT - Research Project
%
% Work on Shawn Kraut & Louis L. Scharf 1998 article for a detection class project
%
% Ref : S. Kraut and L. L. Scharf, "The CFAR adaptive subspace detector is a scale-invariant GLRT," 
% in IEEE Transactions on Signal Processing vol. 47, no. 9, pp. 2538-2541,
%  Sept. 1999, doi: 10.1109/78.782198. 
%
% DETECTION CLASS - SiSy - 21/11/2023 - ISAE-Supaero - Jules GOMEL ☺
% Academic Year 2023/2024 - Prof. S. Bidon

%% Constants
clear all
close all
clc

% False alarm probability Pfa
Pfa = 1E-1;

% Number of antennas/dimension N
N = 10;

% Number of measurements K
K = 1;

% Number of Monte-Carlo samples M
M = 200 / Pfa;

% Deterministic and known steering vector psi: value for alternative hyp H1
psi = 2 * ones(N, 1) + 1i * 2 * ones(N, 1);

%% Generate K synthetic measurements

% Mu depending on the case we want to study
% Change H ☺
H = 0;
switch H
    case 0
        mu = 0 * ones(1, N);
    case 1
        mu = 2 * ones(1, N);
end

% Noise power 
sigma2 = 1E0;

% Random covariance matrix
randR = randn(N);
symR = (randR + randR') / 2;
R = (symR' * symR);

% Generate synthetic measurements!
W = sqrt(1/2) * (randn(N, K, 1) + 1i * randn(N, K, 1));
y = psi + sqrtm(sigma2 * R) * W;

%% Generate M training vectors in N dimensions

% Generate white complex noise 
W = sqrt(1/2) * (randn(N, M, 1) + 1i * randn(N, M, 1));

% Colored noise
X = sqrtm(R) * W;

%% 1) Sample covariance matrix S

S = abs(X * (X') / M);
invS = inv(S);

%% 2) Test statistic of ASD: estimation of cos^2

cos2_hat = 1 / (psi' * invS * psi) * (abs(psi' * invS * y)).^2 ./ (diag(y' * invS * y)');

%% 3) Threshold computation

%% 4) Monte Carlo simulation
