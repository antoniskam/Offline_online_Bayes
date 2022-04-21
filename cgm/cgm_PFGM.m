%{
 On-line solution of the filtering problem for the fatigue crack growth model
 with process noise with the PF and GMM resampling.
%}

clear; clc; close all;
set(0,'DefaultTextInterpreter','latex')

load('cgm_stat_rng1.mat')
%rng(1)
rng('shuffle')

%% setup

% functions to evaluate the state and likelihood
state_fun_x=@(x,N)state_fun(x,N,bound);
state_fun_u=@(u,N) state_fun_x(u2x(u),N);
lh_a=@(a,meas)likelihood_eval(a,meas,mu_meas,sig_meas);

%% filter

nGM=8;

np=5000;           % number of particles
N_T=0.1*np;         % resampling treshold
np_sm=5000;         % number of smoothing particles
np_pred=5000;       % number of prediction particles

res_i=[];           % resampling 'register'
nw_i=[];            % number of used mixtures

% parameter samples
part_u=randn(np,4);
part_x=u2x(part_u);
w=1/np*ones(np,1);

% pre-allocation of parameter estimate statistics in physical space
x_p_mu=zeros(nn,4);
x_p_mu(1,:)=mean(part_x);
x_p_std=zeros(nn,4);
x_p_std(1,:)=std(part_x);
x_p_cov=zeros(nn,4,4);
x_p_cov(1,:,:)=cov(part_x);
x_p_corr=zeros(nn,4,4);
x_p_corr(1,:,:)=corr(part_x);
x_p_median=zeros(nn,4);
x_p_median(1,:)=median(part_x);
x_p_95_l=zeros(nn,4);
x_p_95_l(1,:)=quantile(part_x,0.025);
x_p_95_u=zeros(nn,4);
x_p_95_u(1,:)=quantile(part_x,0.975);
x_p_75_l=zeros(nn,4);
x_p_75_l(1,:)=quantile(part_x,0.125);
x_p_75_u=zeros(nn,4);
x_p_75_u(1,:)=quantile(part_x,0.875);

% pre-allocation of parameter estimate statistics in standard normal space
u_p_mu=zeros(nn,4);
u_p_mu(1,:)=mean(part_u);
u_p_std=zeros(nn,4);
u_p_std(1,:)=std(part_u);
u_p_cov=zeros(nn,4,4);
u_p_cov(1,:,:)=cov(part_u);
u_p_corr=zeros(nn,4,4);
u_p_corr(1,:,:)=corr(part_u);
u_p_median=zeros(nn,4);
u_p_median(1,:)=median(part_u);
u_p_95_l=zeros(nn,4);
u_p_95_l(1,:)=quantile(part_u,0.025);
u_p_95_u=zeros(nn,4);
u_p_95_u(1,:)=quantile(part_u,0.975);
u_p_75_l=zeros(nn,4);
u_p_75_l(1,:)=quantile(part_u,0.125);
u_p_75_u=zeros(nn,4);
u_p_75_u(1,:)=quantile(part_u,0.875);

% pre-allocation of state estimate statistics (physical & log space)
part_a=[part_x(:,1),log(part_x(:,1))];
a_p_mu=zeros(nn,2);
a_p_mu(1,:)=mean(part_a);
a_p_std=zeros(nn,2);
a_p_std(1,:)=std(part_a);
a_p_median=zeros(nn,2);
a_p_median(1,:)=median(part_a);
a_p_95_l=zeros(nn,2);
a_p_95_l(1,:)=quantile(part_a,0.025);
a_p_95_u=zeros(nn,2);
a_p_95_u(1,:)=quantile(part_a,0.975);
a_p_75_l=zeros(nn,2);
a_p_75_l(1,:)=quantile(part_a,0.125);
a_p_75_u=zeros(nn,2);
a_p_75_u(1,:)=quantile(part_a,0.875);

N_eff_seq=zeros(nn,1);  % effective sample size (based on weights)
N_eff_seq(1)=np;
N_part=N_eff_seq;       % number of different particles

logcE=zeros(nn-1,1);     % model evidence
H=[];

tic
for i=2:nn
    
    % likelihood evaluation
    a_cur=state_fun_x(part_x,n(i));
    part_a=[a_cur,log(a_cur)];
    lh_cur=lh_a(a_cur,a_meas(i-1));
    logcE(i-1)=log(w'*lh_cur);
    w=w.*lh_cur;
    w=w/sum(w);
    
    % effective sample size
    N_eff=1/sum(w.^2);
    N_eff_seq(i)=N_eff;
    np_temp=size(part_x,1);
    N_part(i)=length(unique(part_x(:,1)));
    
    % parameter estimate statistics in physical space
    [mu_x,cov_x,std_x,corr_x]=w_stat(part_x,w);
    x_p_mu(i,:)= mu_x;
    x_p_std(i,:)= std_x;
    x_p_cov(i,:,:)= cov_x;
    x_p_corr(i,:,:)= corr_x;
    x_p_median(i,:)=w_qtile(part_x,w,0.5);
    x_p_95_l(i,:)=w_qtile(part_x,w,0.025);
    x_p_95_u(i,:)=w_qtile(part_x,w,0.975);
    x_p_75_l(i,:)=w_qtile(part_x,w,0.125);
    x_p_75_u(i,:)=w_qtile(part_x,w,0.875);
    
    % parameter estimate statistics in standard normal space
    [mu_u,cov_u,std_u,corr_u]=w_stat(part_u,w);
    u_p_mu(i,:)= mu_u;
    u_p_std(i,:)= std_u;
    u_p_cov(i,:,:)= cov_u;
    u_p_corr(i,:,:)= corr_u;
    u_p_median(i,:)=w_qtile(part_u,w,0.5);
    u_p_95_l(i,:)=w_qtile(part_u,w,0.025);
    u_p_95_u(i,:)=w_qtile(part_u,w,0.975);
    u_p_75_l(i,:)=w_qtile(part_u,w,0.125);
    u_p_75_u(i,:)=w_qtile(part_u,w,0.875);
    
    % state estimate statistics (physical & log space)
    [mu_a,~,std_a,~]=w_stat(part_a,w);
    a_p_mu(i,:)= mu_a;
    a_p_std(i,:)= std_a;
    a_p_median(i,:)=w_qtile(part_a,w,0.5);
    a_p_95_l(i,:)=w_qtile(part_a,w,0.025);
    a_p_95_u(i,:)=w_qtile(part_a,w,0.975);
    a_p_75_l(i,:)=w_qtile(part_a,w,0.125);
    a_p_75_u(i,:)=w_qtile(part_a,w,0.875);
    
    if i==n_pred
        part_x_pred=part_x;
        part_u_pred=part_u;
        w_pred=w;
        N_eff_pred=1/sum(w_pred.^2);
    end
    
    % resampling
    if N_eff<N_T && i~=nn
        H1=-w'*log(mvnpdf(part_u,mu_u,cov_u)+realmin);
        res_i=[res_i,i]; %#ok<AGROW>
        
        [mu, si, p] = EMGM(part_u',w,nGM);
        nw=length(p);
        nw_i=[nw_i,nw]; %#ok<AGROW>
        for j=1:np
            indw  = randsample(nw, 1, true,p);
            part_u(j,:) = mvnrnd(mu(:,indw), si(:,:,indw));
        end
        part_x = u2x(part_u);
        w=1/np*ones(np,1); % reweighting
        H2=-w'*log(mvnpdf(part_u,mu_u,cov_u)+realmin);
        H=[H;H1,H2];
        
    end
end

%% plots

%% sample degeneracy

f1=figure('units','centimeters','Position', [5 5 14 7]);
tiledlayout(1,2);

nexttile;
hold on
box on
title('$\;\;\;\;\;$ Sample Degeneracy')
ylabel('$N_{eff}$')
xlabel('$n$')
plot(n,N_eff_seq,'k')
axis([0 n_e 0 1.2*np])

nexttile;
hold on
box on
title('$\;\;\;\;\;$ Unique Parameters')
ylabel('number of unique parameters')
xlabel('$n$')
plot(n,N_part,'k')
axis([0 n_e 0 1.2*np])

%% model evidence

logcE=cumsum(logcE);
f1=figure('units','centimeters','Position', [5 5 7.5 7]);
hold on
box on
plot(n(2:end),logcE,'k')
title('Model Evidence $ln(cE)$')
xlabel('$n$')
ylabel('$ln(cE)$')

%% filtered state

f2=figure('units','centimeters','Position', [5 5 18.5 7]);
tl=tiledlayout(1,2);

nexttile;
hold on
box on
title('Filtered State $\ln(a(n))$')
scatter(n(2:end),lna_meas,10,'g','filled')
plot(n,lna_true,'k','LineWidth',1)
plot(n,a_p_mu(:,2),'c')
plot(n,a_p_95_l(:,2),'r')
plot(n,a_p_75_l(:,2),'b')
plot(n,a_p_95_u(:,2),'r')
plot(n,a_p_75_u(:,2),'b')
axis([0 n_e 0 5])
ylabel('$\ln(a(n))$')
xlabel('$n$')

nexttile;
hold on
box on
title('Filtered State $a(n)$')
scatter(n(2:end),a_meas,10,'g','filled')
plot(n,a_true,'k','LineWidth',1)
plot(n,a_p_mu(:,1),'c')
plot(n,a_p_95_l(:,1),'r')
plot(n,a_p_75_l(:,1),'b')
plot(n,a_p_95_u(:,1),'r')
plot(n,a_p_75_u(:,1),'b')
axis([0 n_e 0 60])
ylabel('$a(n)$')
xlabel('$n$')

lg=legend('measured state','true state','mean state estimate','$95\%$ interval','$75\%$ interval','interpreter','latex','Orientation','Horizontal');
lg.Layout.Tile = 'south';

%% filtered parameters

f3=figure('units','centimeters','Position', [5 5 18.5 7]);
%f3=figure('units','centimeters','Position', [5 5 14 10]);
%f3=figure('units','centimeters','Position', [5 5 14 7]);
tiledlayout(1,4);

nexttile;
hold on
box on
title('$a_0$')
ylabel('$a_0$')
xlabel('$n$')
plot(n2,a0*ones(2,1),'k','LineWidth',1)
plot(n,x_p_mu(:,1),'c')
plot(n,x_p_95_l(:,1),'r')
plot(n,x_p_75_l(:,1),'b')
plot(n,x_p_95_u(:,1),'r')
plot(n,x_p_75_u(:,1),'b')

nexttile;
hold on
box on
title('$\Delta S$')
ylabel('$\Delta S$')
xlabel('$n$')
plot(n2,dS*ones(2,1),'k','LineWidth',1)
plot(n,x_p_mu(:,2),'c')
plot(n,x_p_95_l(:,2),'r')
plot(n,x_p_75_l(:,2),'b')
plot(n,x_p_95_u(:,2),'r')
plot(n,x_p_75_u(:,2),'b')

nexttile;
hold on
box on
title('$C_{ln}$')
ylabel('$C_{ln}$')
xlabel('$n$')
plot(n2,lnC*ones(2,1),'k','LineWidth',1)
plot(n,x_p_mu(:,3),'c')
plot(n,x_p_95_l(:,3),'r')
plot(n,x_p_75_l(:,3),'b')
plot(n,x_p_95_u(:,3),'r')
plot(n,x_p_75_u(:,3),'b')

nexttile;
hold on
box on
title('$m$')
ylabel('$m$')
xlabel('$n$')
plot(n2,m*ones(2,1),'k','LineWidth',1)
plot(n,x_p_mu(:,4),'c')
plot(n,x_p_95_l(:,4),'r')
plot(n,x_p_75_l(:,4),'b')
plot(n,x_p_95_u(:,4),'r')
plot(n,x_p_75_u(:,4),'b')

lg=legend('true value','mean parameter estimate','$95\%$ interval','$75\%$ interval','interpreter','latex','Orientation','Horizontal');
lg.Layout.Tile = 'south';

%% final particles

%f4=figure('units','centimeters','Position', [5 5 14 12]);
f4=figure('units','centimeters','Position', [5 5 15 12]);
%f4=figure('units','centimeters','Position', [5 5 14 7]);

ind = randsample(1:np_temp,np,true,w);              % bootstrap sampling
part_post=part_x(ind,:);
hold on
[S,AX,BigAx,H,HAx] = plotmatrix(part_post);
H(1).NumBins=40;
H(2).NumBins=40;
H(3).NumBins=40;
H(4).NumBins=40;

labels = {'$a_0$' '$\Delta S$' '$C_{ln}$' '$m$'};

for i = 1:4                                       % label the plots
  xlabel(AX(4,i), labels{i})
  ylabel(AX(i,1), labels{i})
end

title('Joint Posterior ($n=10^{7}$)')

%% parameters, prior-posterior

a0_mu_post=x_p_mu(end,1);   a0_sig_post=x_p_std(end,1);
dS_mu_post=x_p_mu(end,2);   dS_sig_post=x_p_std(end,2);
lnC_mu_post=x_p_mu(end,3);   lnC_sig_post=x_p_std(end,3);
m_mu_post=x_p_mu(end,4);   m_sig_post=x_p_std(end,4);

a00=linspace(0,3,100);
a0post=linspace(a0_mu_post-5*a0_sig_post,a0_mu_post+5*a0_sig_post,100);
fa0=w_kde(part_x(:,1),w,a0post);

dS0=linspace(dS_mu-5*dS_sigma,dS_mu+5*dS_sigma,100);
dSpost=linspace(dS_mu_post-5*dS_sig_post,dS_mu_post+5*dS_sig_post,100);
fdS=w_kde(part_x(:,2),w,dSpost);

lnC0=linspace(lnC_mu-5*lnC_sigma,lnC_mu+5*lnC_sigma,100);
lnCpost=linspace(lnC_mu_post-5*lnC_sig_post,lnC_mu_post+5*lnC_sig_post,100);
flnC=w_kde(part_x(:,3),w,lnCpost);

m0=linspace(m_mu-5*m_sigma,m_mu+5*m_sigma,100);
mpost=linspace(m_mu_post-5*m_sig_post,m_mu_post+5*m_sig_post,100);
fm=w_kde(part_x(:,4),w,mpost);

f5=figure('units','centimeters','Position', [5 5 14 12]);
%f5=figure('units','centimeters','Position', [5 5 14 7]);
tiledlayout(2,2);
%sgtitle('Final Marginal Parameter Estimates','interpreter','latex')

nexttile;
hold on
box on
title('$a_0$')
xlabel('$a_0$')
ylabel('$f_{a_0}(a_0)$')
plot(a0post,fa0,'r')
plot(a00,a0_dist.pdf(a00),'b')
scatter(a0,0,20,'k','filled')

nexttile;
hold on
box on
title('$\Delta S$')
xlabel('$\Delta S$')
ylabel('$f_{\Delta S}(\Delta S)$')
plot(dSpost,fdS,'r')
plot(dS0,dS_dist.pdf(dS0),'b')
scatter(dS,0,20,'k','filled')

nexttile;
hold on
box on
title('$C_{ln}$')
xlabel('$C_{ln}$')
ylabel('$f_{C_{ln}}(C_{ln})$')
plot(lnCpost,flnC,'r')
plot(lnC0,lnC_dist.pdf(lnC0),'b')
scatter(lnC,0,20,'k','filled')

nexttile;
hold on
box on
title('$m$')
xlabel('$m$')
ylabel('$f_{m}(m)$')
plot(mpost,fm,'r')
plot(m0,m_dist.pdf(m0),'b')
scatter(m,0,20,'k','filled')

lg=legend('smoothed posterior','prior','true value','interpreter','latex','Orientation','Horizontal','NumColumns',3);
lg.Layout.Tile = 'south';

%% standard deviations

%f6=figure('units','centimeters','Position', [5 5 14 12]);
f6=figure('units','centimeters','Position', [5 5 18.5 7]);
tl = tiledlayout(1,4);
%sgtitle('Filtered Standard Deviations $\sigma$','interpreter','latex')

nexttile;
hold on
box on
plot(n,x_p_std(:,1))
title('$a_0$')
xlabel('$n$')
ylabel('$\sigma_{a_0}$')

nexttile;
hold on
box on
plot(n,x_p_std(:,2))
title('$\Delta S$')
xlabel('$n$')
ylabel('$\sigma_{\Delta S}$')

nexttile;
hold on
box on
plot(n,x_p_std(:,3))
title('$C_{ln}$')
xlabel('$n$')
ylabel('$\sigma_{C_{ln}}$')

nexttile;
hold on
box on
plot(n,x_p_std(:,4))
title('$m$')
xlabel('$n$')
ylabel('$\sigma_m$')

%% correlations

f7=figure('units','centimeters','Position', [5 5 7.5 7]);
hold on
box on
plot(n,x_p_corr(:,1,2))
plot(n,x_p_corr(:,1,3))
plot(n,x_p_corr(:,1,4))
plot(n,x_p_corr(:,2,3))
plot(n,x_p_corr(:,2,4))
plot(n,x_p_corr(:,3,4))
axis([n2 -1 1])
title('Filtered Correlations $\rho$')
xlabel('$n$')
ylabel('$\rho$')

lg=legend('$\rho_{a_0, \Delta S}$','$\rho_{a_0, C_{ln}}$','$\rho_{a_0, m}$','$\rho_{\Delta S, C_{ln}}$','$\rho_{\Delta S, m}$','$\rho_{C_{ln}, m}$','interpreter','latex','location','north','Orientation','Horizontal','NumColumns',3);
%lg.Layout.Tile = 'south';

%% smoothing

% state computation
state_fun_x=@(x,N)state_fun(x,N,bound);
a_sm=state_fun_x(part_x,n);
w_sm=w;

[a_sm_mu,a_sm_std]=w_mom(a_sm,w_sm);
a_sm_median=w_qtile(a_sm,w_sm,0.5);
a_sm_95_l=w_qtile(a_sm,w_sm,0.025);
a_sm_95_u=w_qtile(a_sm,w_sm,0.975);
a_sm_75_l=w_qtile(a_sm,w_sm,0.125);
a_sm_75_u=w_qtile(a_sm,w_sm,0.875);

lna_sm=log(a_sm);
[lna_sm_mu,lna_sm_std]=w_mom(lna_sm,w_sm);
lna_sm_median=w_qtile(lna_sm,w_sm,0.5);
lna_sm_95_l=w_qtile(lna_sm,w_sm,0.025);
lna_sm_95_u=w_qtile(lna_sm,w_sm,0.975);
lna_sm_75_l=w_qtile(lna_sm,w_sm,0.125);
lna_sm_75_u=w_qtile(lna_sm,w_sm,0.875);


%%%%%%%%%%%%%%%%%

f8=figure('units','centimeters','Position', [5 5 18.5 7]);
tl = tiledlayout(1,2);

% f8=figure('units','centimeters','Position', [5 5 14 7]);
% tl = tiledlayout(1,2);

nexttile;
box on
hold on
title('Smoothed State $\ln(a(n))$')
scatter(n(2:end),lna_meas,10,'g','filled')
plot(n,lna_true,'k','LineWidth',1)
plot(n,lna_sm_mu,'c')
plot(n,lna_sm_95_l,'r')
plot(n,lna_sm_75_l,'b')
plot(n,lna_sm_95_u,'r')
plot(n,lna_sm_75_u,'b')
axis([0 n_e 0 max([lna_meas(end),lna_true(end)])*1.05])
ylabel('$ln(a(n))$')
xlabel('$n$')

nexttile;
box on
hold on
title('Smoothed State $a(n)$')
scatter(n(2:end),a_meas,10,'g','filled')
plot(n,a_true,'k','LineWidth',1)
plot(n,a_sm_mu,'c')
plot(n,a_sm_95_l,'r')
plot(n,a_sm_75_l,'b')
plot(n,a_sm_95_u,'r')
plot(n,a_sm_75_u,'b')
axis([0 n_e 0 max([a_meas(end),a_true(end)])])
ylabel('$a(n)$')
xlabel('$n$')

lg=legend('measured state','true state','mean state estimate','$95\%$ interval','$75\%$ interval','interpreter','latex','Orientation','Horizontal');
lg.Layout.Tile = 'south';

%% prediction

% state computation
state_fun_x=@(x,N)state_fun(x,N,bound);
a_sm=state_fun_x(part_x_pred,n);

[a_sm_mu,a_sm_std]=w_mom(a_sm,w_pred);
a_sm_median=w_qtile(a_sm,w_pred,0.5);
a_sm_95_l=w_qtile(a_sm,w_pred,0.025);
a_sm_95_u=w_qtile(a_sm,w_pred,0.975);
a_sm_75_l=w_qtile(a_sm,w_pred,0.125);
a_sm_75_u=w_qtile(a_sm,w_pred,0.875);

lna_sm=log(a_sm);
[lna_sm_mu,lna_sm_std]=w_mom(lna_sm,w_pred);
lna_sm_median=w_qtile(lna_sm,w_pred,0.5);
lna_sm_95_l=w_qtile(lna_sm,w_pred,0.025);
lna_sm_95_u=w_qtile(lna_sm,w_pred,0.975);
lna_sm_75_l=w_qtile(lna_sm,w_pred,0.125);
lna_sm_75_u=w_qtile(lna_sm,w_pred,0.875);

%%%%%%%%%%%%%%%%%

%f9=figure('units','centimeters','Position', [5 -10 14 24]);
%tl = tiledlayout(2,1);

f9=figure('units','centimeters','Position', [5 5 18.5 7]);
tl = tiledlayout(1,2);

nexttile;
box on
hold on
title('Predicted State $\ln(a(n))$')
scatter(n(2:51),lna_meas(1:50),10,'g','filled')
plot(n,lna_true,'k','LineWidth',1)
plot(n,lna_sm_mu,'c')
plot(n,lna_sm_95_l,'r')
plot(n,lna_sm_75_l,'b')
plot(n,lna_sm_95_u,'r')
plot(n,lna_sm_75_u,'b')
axis([0 n_e 0 max([lna_meas(end),lna_true(end)])*1.05])
ylabel('$ln(a(n))$')
xlabel('$n$')

nexttile;
box on
hold on
title('Predicted State $a(n)$')
scatter(n(2:51),a_meas(1:50),10,'g','filled')
plot(n,a_true,'k','LineWidth',1)
plot(n,a_sm_mu,'c')
plot(n,a_sm_95_l,'r')
plot(n,a_sm_75_l,'b')
plot(n,a_sm_95_u,'r')
plot(n,a_sm_75_u,'b')
axis([0 n_e 0 max([a_meas(end),a_true(end)])])
ylabel('$a(n)$')
xlabel('$n$')

lg=legend('measured state','true state','mean state estimate','$95\%$ interval','$75\%$ interval','interpreter','latex','Orientation','Horizontal');
lg.Layout.Tile = 'south';

%% nested functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sf=state_fun(x,N,bound)

nN=length(N);
ns=size(x,1);
sf=zeros(ns,nN);
for i=1:nN
    a=((1-x(:,4)/2).*exp(x(:,3)).*x(:,2).^x(:,4).*pi.^(x(:,4)/2).*N(i)+x(:,1).^(1-x(:,4)/2)).^(1./(1-x(:,4)/2));
    ind=imag(a)~=0 | a>bound;
    a(ind)=bound;
    sf(:,i)=a;
end

end

function lh=likelihood_eval(a,meas,mu_meas,sig_meas)

lh=prod(normpdf(log(a)-log(meas),mu_meas,sig_meas),2);

end