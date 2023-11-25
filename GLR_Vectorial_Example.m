% 
% GLR Detector - Vectorial Example
%
% Threshold and Detection probability for GLRT detector w/ MC Simulation
%
% DETECTION - SiSy - 14/11/2023 - Jules GOMEL ☺
% AY 2023/2024 - Prof. S. Bidon

%% Constants
clear all
close all
clc

% False alarm probability Pfa
Pfa=1E-2;

% Number of antennas N, dimension
N=10;
% Noise power + Cov matrix of noise
sigma2=1E-1;

sampleR=randn(N);
R=sigma2*(sampleR'*sampleR);
invR=inv(R);
mu=zeros(1,N);
% a deterministic and known : value for alternative hyp H1
a=2*ones(N,1)+1i*2*ones(N,1);

% Number of Monte-Carlo samples K
K=200/Pfa;


%% Pfa-Threshold

% Generate samples following H0
x_H0=sqrt(1/2)*(mvnrnd(mu,R,K)'+1i*mvnrnd(mu,R,K)');

% Derive test statistics (in our case, t=x)
t_H0=sort(abs(a'*invR*x_H0).^2/(a'*(invR*a)));

% Compute theoretical threshold
eta_th=icdf('Exponential', 1-Pfa, 1);

% Find exp eta so that Nbr_fa=MC*Pfa
eta_exp=t_H0(floor(K*(1-Pfa)));

% Check eta_exp=eta_th or close to
clc
disp('Ecart entre eta théorique et expérimental) : ')
display(abs(eta_th-eta_exp))
%% Eta as a function of Pfa
Pfa_vec=1E-4:1E-4:1E-1;
eta_th_array=zeros(1,length(Pfa_vec));

for p=1:length(Pfa_vec)
    Pfa=Pfa_vec(p);
    
    % Compute theoretical threshold
    eta_th_array(p)=icdf('Exponential', 1-Pfa, 1);
end

figure
hold on
xlabel('P_{fa}')
ylabel('\eta')
title('Threshold as a function of Probability of False alarm')
plot(Pfa_vec,eta_th_array,color='b',marker='pentagram')
hold off

%% Probability of detection

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

%% Signal-To-Noise Ratio!

SNR_opt=20*log10(((a'*a)*ones(size(sigma2_vec)))./sigma2_vec);

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