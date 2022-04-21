%{
 Off-line solution of the filtering problem for the fatigue crack growth model
 with the IMH-GM-based SMC filter
%}

clear; clc; close all;
rng(10)
set(0,'DefaultTextInterpreter','latex')
load('cgm_stat_rng1.mat')

%% filter
rng('shuffle')

N =1000;       % number of samples per level
w=1/N*ones(N,1);

state_fun_x=@(x,N)state_fun(x,N,bound);
lh_a=@(a,meas)likelihood_eval(a,meas,mu_meas,sig_meas);

lh_x_cur=@(x)lh_a(state_fun_x(x,n(2:end)),a_meas);
log_likelihood_cur = @(theta) log(lh_x_cur(theta)+realmin);

nGM=8;
burn=0;
tarCoV=sqrt(1/0.5-1);
p=1;

[samplesU, samplesX, q, k_fin, logcE] = SMC_GM(N, p, log_likelihood_cur, prior, nGM, burn, tarCoV);

part_x=samplesX{end};
part_u=samplesU{end};

logcE

%%

x_p_mu=mean(part_x)
x_p_std=std(part_x)
x_p_skewness=skewness(part_x)
x_p_kurtosis=kurtosis(part_x)
x_p_corr=corr(part_x)

%%%%%%%%%%%%%%%%%

figure()
tiledlayout(2,2);
sgtitle('Final estimates for time invariant parameters (physical space)','interpreter','latex')

nexttile;
hold on
xlabel('$a_0$')
ylabel('$f(a_0)$')
histogram(part_x(:,1),'Normalization','PDF');
xi_x1=0:0.05:2.5;
f_eval=w_kde(part_x(:,1),w,xi_x1);
plot(xi_x1,f_eval,'r','LineWidth',1)
plot(xi_x1,a0_dist.pdf(xi_x1),'k','LineWidth',1);
scatter(a0_mu,0,'k','filled')
scatter(x_p_mu(end,1),0,'r','filled')
scatter(a0,0,15,'b','filled')

nexttile;
hold on
xlabel('$dS$')
ylabel('$f(dS)$')
histogram(part_x(:,2),'Normalization','PDF');
xi_x2=20:0.1:100;
f_eval=w_kde(part_x(:,2),w,xi_x2);
plot(xi_x2,f_eval,'r','LineWidth',1)
plot(xi_x2,dS_dist.pdf(xi_x2),'k','LineWidth',1);
scatter(dS_mu,0,'k','filled')
scatter(x_p_mu(end,2),0,'r','filled')
scatter(dS,0,15,'b','filled')

nexttile;
hold on
xlabel('$ln(C)$')
ylabel('$f(ln(C))$')
histogram(part_x(:,3),'Normalization','PDF');
xi_x3=-34.5:0.01:-31.5;
f_eval=w_kde(part_x(:,3),w,xi_x3);
plot(xi_x3,f_eval,'r','LineWidth',1)
plot(xi_x3,lnC_dist.pdf(xi_x3),'k','LineWidth',1);
scatter(lnC_mu,0,'k','filled')
scatter(x_p_mu(end,3),0,'r','filled')
scatter(lnC,0,15,'b','filled')

nexttile;
hold on
xlabel('$m$')
ylabel('$f(m)$')
histogram(part_x(:,4),'Normalization','PDF');
xi_x4=2.6:0.005:4.4;
f_eval=w_kde(part_x(:,4),w,xi_x4);
plot(xi_x4,f_eval,'r','LineWidth',1)
plot(xi_x4,m_dist.pdf(xi_x4),'k','LineWidth',1);
scatter(m_mu,0,'k','filled')
scatter(x_p_mu(end,4),0,'r','filled')
scatter(m,0,15,'b','filled')

lg=legend('normalized histogram of posterior','smoothed posterior','prior','prior mean','posterior mean','true value','interpreter','latex','Orientation','Horizontal','NumColumns',3);
lg.Layout.Tile = 'south';

%% smoothing

part_sm=part_x;
w_sm=w;

% state computation
a_sm=state_fun_x(part_sm,n);

% computation of estimates
[a_sm_mu,~,a_sm_std,~]=w_stat(a_sm,w_sm);
a_sm_median=w_qtile(a_sm,w_sm,0.5);
a_sm_95_l=w_qtile(a_sm,w_sm,0.025);
a_sm_95_u=w_qtile(a_sm,w_sm,0.975);
a_sm_75_l=w_qtile(a_sm,w_sm,0.125);
a_sm_75_u=w_qtile(a_sm,w_sm,0.875);

lna_sm=log(a_sm);
[lna_sm_mu,~,lna_sm_std,~]=w_stat(lna_sm,w_sm);
lna_sm_median=w_qtile(lna_sm,w_sm,0.5);
lna_sm_95_l=w_qtile(lna_sm,w_sm,0.025);
lna_sm_95_u=w_qtile(lna_sm,w_sm,0.975);
lna_sm_75_l=w_qtile(lna_sm,w_sm,0.125);
lna_sm_75_u=w_qtile(lna_sm,w_sm,0.875);

%%%%%%%%%%%%%%%%%
figure()

subplot(1,2,1)
hold on
title('Smoothing of $a$')
scatter(n(2:end),a_meas,15,'g','filled')
plot(n,a_sm_mu,'c')
plot(n,a_sm_median,'color',[0.3,0.3,0.3],'LineStyle',':')
plot(n,a_true,'k','LineWidth',1)
plot(n,a_sm_95_l,'r')
plot(n,a_sm_75_l,'b')
plot(n,a_sm_95_u,'r')
plot(n,a_sm_75_u,'b')
legend('measured state','mean estimated state','estimated median state','true state $a$','$95\%$ interval','$75\%$ interval','interpreter','latex','location','northwest')
axis([0 n_e 0 max([a_meas(end),a_true(end)])])
ylabel('$a$')
xlabel('$N$')

%figure()
subplot(1,2,2)
hold on
title('Smoothing of $ln(a)$')
scatter(n(2:end),lna_meas,15,'g','filled')
plot(n,lna_sm_mu,'c')
plot(n,lna_sm_median,'color',[0.3,0.3,0.3],'LineStyle',':')
plot(n,lna_true,'k','LineWidth',1)
plot(n,lna_sm_95_l,'r')
plot(n,lna_sm_75_l,'b')
plot(n,lna_sm_95_u,'r')
plot(n,lna_sm_75_u,'b')
legend('measured state','mean estimated state','estimated median state','true state $a$','$95\%$ interval','$75\%$ interval','interpreter','latex','location','northwest')
axis([0 n_e 0 max([lna_meas(end),lna_true(end)])*1.05])
ylabel('$ln(a)$')
xlabel('$N$')

%%

f4=figure('units','centimeters','Position', [5 5 15 12]);
hold on
[S,AX,BigAx,H,HAx] = plotmatrix(part_x);
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