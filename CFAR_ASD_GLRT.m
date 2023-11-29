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

% Number of measurements K : Number of Monte-Carlo samples
K = 200 / Pfa;

% Number of training vectors M
M = 5000;

% Number for MC-mean
MC=100;

% Deterministic and known steering vector psi: value for alternative hyp H1
psi = 2 * ones(N, 1) + 1i * 2 * ones(N, 1);

% Noise power 
sigma2 = 1E0;

% Random covariance matrix
randR = randn(N);
symR = (randR + randR') / 2;
R = (symR' * symR);

% Generate M training vectors in N-dimension

% Generate white complex noise 
W = sqrt(1/2) * (randn(N, M, 1) + 1i * randn(N, M, 1));

% Colored noise
X = sqrtm(R) * W;

S = X * (X') / M;
invS = inv(S);

%% Part 4.1) Eta as a function of Pfa
Pfa_vec=1E-3:1E-3:9E-1;
eta_exp_array=zeros(length(Pfa_vec),MC);

for p=1:length(Pfa_vec)
    for m=1:MC
        % False alarm rate
        Pfa=Pfa_vec(p);
    
        % Number of measurements K : Number of Monte-Carlo samples
        K = floor(200 / Pfa);
    
        % Generate synthetic measurements under H0! (mu=0)
        W = sqrt(1/2) * (randn(N, K, 1) + 1i * randn(N, K, 1));
        y = sqrtm(sigma2 * R) * W;
        
        % Real because sometimes complex with null imag part 
        cos2_hat = real(1 / (psi' * invS * psi) * (abs(psi' * invS * y)).^2 ./ denom_cos2(y,invS,K));
    
        % Sort the test statistic 
        cos2_hat=sort(cos2_hat);
        
        % Find exp eta so that Nbr_fa=MC*Pfa
        eta_exp_array(p,m)=cos2_hat(floor(K*(1-Pfa)));
    end
end

%% Display eta=f(Pfa)

%eta_exp_array=mean(eta_exp_array,2);

figure
hold on
xlabel('P_{fa}')
ylabel('\eta')
title('Threshold as a function of Probability of False alarm')
plot(Pfa_vec,eta_exp_array,color='b',LineWidth=1)
hold off

%% 3)b) Pd as a function of SNR for a fixed Pfa
Pfa_vec=[1E-1,5E-2,1E-2];

sigma2_vec=1E-1:1E0:1E2;
Pd_exp_vec=zeros(length(Pfa_vec),length(sigma2_vec),MC);

for p=1:length(Pfa_vec)
    for i=1:length(sigma2_vec)
        for m=1:MC
            % Extract Pfa from vector
            Pfa=Pfa_vec(p);
            sigma2=sigma2_vec(i);
            
            % Number of measurements K : Number of Monte-Carlo samples
            K = floor(200 / Pfa);
    
            % Generate synthetic measurements under H0
            W = sqrt(1/2) * (randn(N, K, 1) + 1i * randn(N, K, 1));
            y_H0 = sqrtm(sigma2 * R) * W;
    
            % Real because sometimes complex with null imag part 
            cos2_hat_H0 = real(1 / (psi' * invS * psi) * (abs(psi' * invS * y_H0)).^2 ./ denom_cos2(y_H0,invS,K));
    
            % Sort the test statistic 
            cos2_hat_H0=sort(cos2_hat_H0);
    
            % Find exp eta so that Nbr_fa=MC*Pfa
            eta_exp=cos2_hat_H0(floor(K*(1-Pfa)));
            
            % Generate synthetic measurements under H1
            W = sqrt(1/2) * (randn(N, K, 1) + 1i * randn(N, K, 1));
            y_H1 = psi + sqrtm(sigma2 * R) * W;
            
            % Real because sometimes complex with null imag part 
            cos2_hat_H1 = real(1 / (psi' * invS * psi) * (abs(psi' * invS * y_H1)).^2 ./ denom_cos2(y_H1,invS,K));
    
            
            Pd_exp_vec(p,i,m)=length(find(cos2_hat_H1>eta_exp))/K;
    
        end
    end
end


SNR=psi'*inv(R)*psi./sigma2_vec;
%%  Display

figure
hold on

plot(SNR, Pd_exp_vec(1, :), 'r')
plot(SNR, Pd_exp_vec(2, :), 'g')
plot(SNR, Pd_exp_vec(3, :), 'b')
xlim([0 100])
xlabel('SNR')
ylabel('P_d')
title('P_d as a function of SNR, for a given P_{fa}')
legend('P_{fa}=1E-1', 'P_{fa}=5E-2', 'P_{fa}=1E-2', 'P_{fa}=1E-3')

hold off

%% ROC - Pd = f(Pfa)

Pfa_vec=1E-3:1E-3:1E-1;
sigma2=20;
Pd_ROC_vec=zeros(length(Pfa_vec),MC);

for p=1:length(Pfa_vec)
    for m=1:MC
        % Extract Pfa from vector
        Pfa=Pfa_vec(p);
        % Number of measurements K : Number of Monte-Carlo samples
        K = floor(200 / Pfa);

        % Generate synthetic measurements under H0
        W = sqrt(1/2) * (randn(N, K, 1) + 1i * randn(N, K, 1));
        y_H0 = sqrtm(sigma2 * R) * W;

        % Real because sometimes complex with null imag part 
        cos2_hat_H0 = real(1 / (psi' * invS * psi) * (abs(psi' * invS * y_H0)).^2 ./ denom_cos2(y_H0,invS,K));

        % Sort the test statistic 
        cos2_hat_H0=sort(cos2_hat_H0);

        % Find exp eta so that Nbr_fa=MC*Pfa
        eta_exp=cos2_hat_H0(floor(K*(1-Pfa)));
        
        % Generate synthetic measurements under H1
        W = sqrt(1/2) * (randn(N, K, 1) + 1i * randn(N, K, 1));
        y_H1 = psi + sqrtm(sigma2 * R) * W;
        
        % Real because sometimes complex with null imag part 
        cos2_hat_H1 = real(1 / (psi' * invS * psi) * (abs(psi' * invS * y_H1)).^2 ./ denom_cos2(y_H1,invS,K));

        
        Pd_ROC_vec(p,m)=length(find(cos2_hat_H1>eta_exp))/K;
    end
end

%% Display ROC

Pd_ROC_vec=mean(Pd_ROC_vec,2);

figure
hold on
xlabel('P_{fa}')
ylabel('P_d')
title('P_d as function of P_{fa}, fixed SNR')
plot(Pfa_vec,Pd_ROC_vec(:),color='b',LineWidth=1)
%xlim([0 0.01])
legend('P_d')
hold off

%% Functions

function den=denom_cos2(y,invS,K)
    den=zeros(1,K);
    for k=1:K
        %k-th column
        y_k=y(:,k);
        % Direct computation of the diag of y'S-1y to avoid huge arrays
        den(k)=y_k'*invS*y_k;
    end
end



