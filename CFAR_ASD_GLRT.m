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
Pfa=1E-2;

% Number of antennas/dimension N
N=10;

% Number of Monte-Carlo samples M
M=200/Pfa;

% Noise power + Cov matrix of noise
sigma2=1E-1;

sampleR=randn(N);
R=sigma2*(sampleR'*sampleR);
invR=inv(R);

% Zero mean
mu=zeros(1,N);
% psi deterministic and known : value for alternative hyp H1
psi=2*ones(N,1)+1i*2*ones(N,1);

% For synthetic measurement data : fi



%% Generate M training vectors, N dim

% Generate white complex noise 
W=sqrt(1/2)*(randn(N,M,1)+1i*randn(N,M,1));

% Colored noise
X=sqrtm(R)*W;

%% Training vector sample covariance matrix S

S=X*(X')/M;

% S is not perfectly real

re_im_ratio=real(S)./imag(S);

% But Re>>Im so we can take abs

S=abs(S);


%% Estimation of cos^2


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