% 
% Neyman-Pearson Detector - Vectorial Example
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
Pfa=5E-3;

% Number of antennas N, dimension
N=10;
% Noise power
sigma2=10;
% a deterministic and known : value for alternative hyp H1
a=2*ones(N,1);

% Number of Monte-Carlo samples K
K=200/Pfa;

%% Pfa-Threshold

% Generate samples following H0
x_H0=sqrt(sigma2/2)*(randn(N,K)+1i*randn(N,K));

% Derive test statistics (in our case, t=x)
t_H0=sort(real(a'*x_H0)/(sqrt(.5*sigma2*(a'*a))));

% Compute theoretical threshold
eta_th=norminv(1-Pfa);

% Find exp eta so that Nbr_fa=MC*Pfa
eta_exp=t_H0(floor(K*(1-Pfa)));
%% Detection Probability

% Generate samples following H1
x_H1=a+sqrt(sigma2/2)*(randn(N,K)+1i*randn(N,K));

% Derive test statistics (in our case, t=x)
t_H1=sort(real(a'*x_H1)/(sqrt(.5*sigma2*(a'*a))));

Pd_th=1-normcdf(eta_th-sqrt(2*(a'*a)/sigma2));
Pd_exp=length(find(t_H1>eta_exp))/K;

%% Signal-To-Noise Ratio!
Pfa_vec=[1E-1,5E-2,1E-2,1E-3];

%sigma2_vec=[1E-1,1,2,4,1E1,20,1E2,1E3,5E3,1E4];
sigma2_vec=1E-1:1:5E2;
Pd_th_vec=zeros(length(Pfa_vec),length(sigma2_vec));

for p=1:length(Pfa_vec)
    for i=1:length(sigma2_vec)
        % Extract Pfa from vector
        Pfa=Pfa_vec(p);
        % Generate samples following H0
        x_H0=sqrt(sigma2_vec(i)/2)*(randn(N,K)+1i*randn(N,K));

        % Derive test statistics (in our case, t=x)
        t_H0=sort(real(a'*x_H0)/(sqrt(.5*sigma2_vec(i)*(a'*a))));

        % Compute theoretical threshold
        eta_th=norminv(1-Pfa);

        % Find exp eta so that Nbr_fa=MC*Pfa
        eta_exp=t_H0(floor(K*(1-Pfa)));

        % Generate samples following H1
        x_H1=a+sqrt(sigma2_vec(i)/2)*(randn(N,K)+1i*randn(N,K));

        % Derive test statistics (in our case, t=x)
        t_H1=sort(real(a'*x_H1)/(sqrt(.5*sigma2_vec(i)*(a'*a))));

        Pd_th_vec(p,i)=1-normcdf(eta_th-sqrt(2*(a'*a)/sigma2_vec(i)));
        %Pd_exp=length(find(t_H1>eta_exp))/K;
    end
end
SNR_opt=10*log10(((a'*a)*ones(size(sigma2_vec)))./sigma2_vec);

%% ROC ! Receiver Operating Curve

Pfa_vec=1E-3:1E-3:1E-1;
sigma2=5;
Pd_ROC_vec=zeros(1,length(Pfa_vec));

for p=1:length(Pfa_vec)
        % Extract Pfa from vector
        Pfa=Pfa_vec(p);
        % Generate samples following H0
        x_H0=sqrt(sigma2/2)*(randn(N,K)+1i*randn(N,K));

        % Derive test statistics (in our case, t=x)
        t_H0=sort(real(a'*x_H0)/(sqrt(.5*sigma2*(a'*a))));

        % Compute theoretical threshold
        eta_th=norminv(1-Pfa);

        % Find exp eta so that Nbr_fa=MC*Pfa
        eta_exp=t_H0(floor(K*(1-Pfa)));

        % Generate samples following H1
        x_H1=a+sqrt(sigma2/2)*(randn(N,K)+1i*randn(N,K));

        % Derive test statistics (in our case, t=x)
        t_H1=sort(real(a'*x_H1)/(sqrt(.5*sigma2*(a'*a))));

        Pd_ROC_vec(p)=1-normcdf(eta_th-sqrt(2*(a'*a)/sigma2));
        %Pd_exp=length(find(t_H1>eta_exp))/K;
end

%% Display

close all
clc


figure
hold on
xlabel('SNR_{opt} (dB)')
ylabel('P_d')
title('P_d as function of the optimal SNR, for a given P_{fa}=1E-1')
plot(SNR_opt,Pd_th_vec(1,:),color='b',marker='*')
legend('P_d')
hold off

figure
hold on
xlabel('SNR_{opt} (dB)')
ylabel('P_d')
title('P_d as function of the optimal SNR, for a given P_{fa}')
plot(SNR_opt,Pd_th_vec(1,:),color='r',marker='*')
plot(SNR_opt,Pd_th_vec(2,:),color='g',marker='*')
plot(SNR_opt,Pd_th_vec(3,:),color='b',marker='*')
plot(SNR_opt,Pd_th_vec(4,:),color='k',marker='*')
legend('P_{fa}=1E-1','P_{fa}=5E-2','P_{fa}=1E-2','P_{fa}=1E-3')
hold off

figure
hold on
xlabel('P_{fa}')
ylabel('P_d')
title('P_d as function of P_{fa}, fixed SNR')
plot(Pfa_vec,Pd_ROC_vec(:),color='b',marker='*')
legend('P_d')
hold off