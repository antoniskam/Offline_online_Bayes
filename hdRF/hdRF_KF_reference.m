clear; clc; close all;
rng(40)
set(0,'DefaultTextInterpreter','latex')

%% Parameters of the problem

% plot at time
t_plot=50;
ind_plot=t_plot+1;

% measurement noise
mu_meas=0; std_meas=0.1;

% incorporated measurements
ind_s=1:10;
%ind_s=[1,4,7,10];
%ind_s=[4,7];
ns_sol=length(ind_s);

% number of elements for solution
% n_ele=25;
% mp=5; np=5;

n_ele=50;
mp=8; np=7;

% n_ele=75;
% mp=10; np=8;

% n_ele=100;
% mp=10; np=10;

% correlation length
l = 2;
% number of KL-modes
mkl = 400;

% length of domain
L = 4;
dx=L/n_ele;
% sensor location
xs = 0.2:0.4:3.8;
ns=length(xs);

nt=50;
t=1:nt;
t0=0:nt;

cl=0.05; cu=0.95;

% Solve Integral eigenvalue problem for the exponential kernel
[lambda,phi] = EigFcnKL(mkl,L,l);
err_var = 1-sum(lambda)/L;

CoV_A=0.3; sig_lnA=sqrt(log(CoV_A^2+1));
mu_A = 0.8; mu_lnA=log(mu_A)-sig_lnA^2/2;
mu_B=0.8; CoV_B=0.15; std_B=mu_B*CoV_B;
mu_lnD = mu_lnA+log(t)*mu_B;
sig_lnD = sqrt(sig_lnA^2+log(t).^2*std_B^2);
mu_D = exp(mu_lnD+sig_lnD.^2/2);
mu_lnD_final = mu_lnA+log(50)*mu_B;
sig_lnD_final = sqrt(sig_lnA^2+log(50)^2*std_B^2);
mu_D_final = exp(mu_lnD_final+sig_lnD_final^2/2);

% compute kl values
kl_s=zeros(ns,mkl);
for i = 1:mkl 
    kl_s(:,i) = sqrt(lambda(i))*phi{i}(xs); 
end

% compute kl values
nf=2000;
xf = linspace(0,L,nf);
kl_f=zeros(nf,mkl);
for i = 1:mkl 
    kl_f(:,i) = sqrt(lambda(i))*phi{i}(xf); 
end

% compute RF realization for A
xi_A = randn(mkl,1);
lnA_s = mu_lnA+sig_lnA*kl_s*xi_A;
lnA_f = mu_lnA+sig_lnA*kl_f*xi_A;
A_s = exp(lnA_s);
A_f = exp(lnA_f);

sig_lnA_s=sqrt(sum((sig_lnA*kl_s).^2,2));
sig_lnA_f=sqrt(sum((sig_lnA*kl_f).^2,2));

% compute RF realization for B
xi_B = randn(mkl,1);
xi_true = [xi_A; xi_B];
B_s = mu_B+std_B*kl_s*xi_B;
B_f = mu_B+std_B*kl_f*xi_B;

sig_B_s=sqrt(sum((std_B*kl_s).^2,2));
sig_B_f=sqrt(sum((std_B*kl_f).^2,2));

% compute RF realization for D
D_s=A_s.*t.^B_s;
D_f=A_f.*t.^B_f;
lnD_s=log(D_s);
lnD_f=log(D_f);

meas_noise=normrnd(mu_meas,std_meas,ns,nt);
D_s_meas=D_s.*exp(meas_noise);
lnD_s_meas=log(D_s_meas);

sig_lnD_s=sqrt(sig_lnA_s.^2+log(t).^2.*sig_B_s.^2);
sig_lnD_f=sqrt(sig_lnA_f.^2+log(t).^2.*sig_B_f.^2);

%%

figure()
subplot(1,3,1)
hold on
plot(xf,lnA_f,'color','#0072BD')
plot([0,L],[norminv(cu,mu_lnA,sig_lnA),norminv(cu,mu_lnA,sig_lnA)],'k')
plot([0,L],[norminv(cl,mu_lnA,sig_lnA),norminv(cl,mu_lnA,sig_lnA)],'k')
plot([0,L],[mu_lnA,mu_lnA],'k','linestyle','--')
axis([0 L mu_lnA-3*sig_lnA mu_lnA+3*sig_lnA])
title('$\ln(A)$')

subplot(1,3,2)
hold on
plot(xf,B_f,'color','#0072BD')
plot([0,L],[norminv(cu,mu_B,std_B),norminv(cu,mu_B,std_B)],'k')
plot([0,L],[norminv(cl,mu_B,std_B),norminv(cl,mu_B,std_B)],'k')
plot([0,L],[mu_B,mu_B],'k','linestyle','--')
axis([0 L mu_B-3*std_B mu_B+3*std_B])
title('$B$')

subplot(1,3,3)
hold on
plot(xf,lnD_f(:,end),'color','#0072BD')
plot([0,L],[norminv(cu,mu_lnD_final,sig_lnD_final),norminv(cu,mu_lnD_final,sig_lnD_final)],'k')
plot([0,L],[norminv(cl,mu_lnD_final,sig_lnD_final),norminv(cl,mu_lnD_final,sig_lnD_final)],'k')
plot([0,L],[mu_lnD_final,mu_lnD_final],'k','linestyle','--')
axis([0 L mu_lnD(end)-3*sig_lnD(end) mu_lnD(end)+3*sig_lnD(end)])
title('$\ln(D(t=50))$')

%% 

figure()
sgtitle('Prior $\ln(D(t))$')
for i=1:ns
    subplot(2,5,i)
    hold on
    plot(t,lnD_s(i,:),'color','#0072BD')
    plot(t,norminv(cu,mu_lnD,sig_lnD),'k')
    plot(t,norminv(cl,mu_lnD,sig_lnD),'k')
    plot(t,mu_lnD,'k','linestyle','--')
    scatter(t,lnD_s_meas(i,:),5,'filled','k')
    xlim([0,50])
    title(['Location ',num2str(i)])
    xlabel('$t$')
    ylabel('$\ln(D(t))$')
end

%% filter

xm = dx/2:dx:L;
x0 = xm+1e-8-dx/2;
x1 = xm+1e-8+dx/2;
Rho=zeros(n_ele);
for i=1:n_ele
    Rho(i,:)=exp(-abs(xm(i)-xm)/l);
end
Rho_comp=[Rho,zeros(n_ele);zeros(n_ele),Rho];
mu_par=[mu_lnA*ones(n_ele,1);mu_B*ones(n_ele,1)];
sig_par=[sig_lnA*ones(n_ele,1);std_B*ones(n_ele,1)];
cov_par=(sig_par*sig_par').*Rho_comp;
chol_L=chol(cov_par)';

ind_sm=zeros(1,ns);
for i=1:ns
    ind_sm(i)=find(x0<xs(i) & xs(i)<=x1);
end
ind_sm_sol=ind_sm(ind_s);

mu_u_t=zeros(nt+1,n_ele*2);
sig_u_t=ones(nt+1,n_ele*2);
cov_u_t=cell(nt+1,1);
corr_u_t=cell(nt+1,1);
cov_u_t{1}=eye(n_ele*2);
corr_u_t{1}=eye(n_ele*2);
corr_u_t_vec=zeros(nt+1,n_ele);

mu_x_t=zeros(nt+1,n_ele*3);
sig_x_t=zeros(nt+1,n_ele*3);
cov_x_t=cell(nt+1,1);
corr_x_t=cell(nt+1,1);
corr_x_t_vec=zeros(nt+1,n_ele);

logcE_t=zeros(nt+1,1);

R=eye(ns_sol)*std_meas^2;

x=mu_par;
P=cov_par;

mu_x_t(1,:)=[mu_par; -Inf*ones(n_ele,1)];
sig_x_t(1,:)=[sig_par; Inf*ones(n_ele,1)];
cov_x_t{1}=cov_par;
corr_x_t{1}=Rho_comp;

A_ind=[eye(n_ele),log(t(t_plot))*eye(n_ele)];
mu_x_t_lnD_ind=zeros(nt+1,n_ele);
sig_x_t_lnD_ind=zeros(nt+1,n_ele);
mu_x_t_lnD_ind(1,:)=A_ind*x;
sig_x_t_lnD_ind(1,:)=sqrt(diag(A_ind*P*A_ind'));

for i=2:nt+1

    % prediction state
    a=[eye(2*n_ele);eye(n_ele),eye(n_ele)*log(t(i-1))];
    x=a*x;
    P=a*P*a';
    
    % measurement equation
    H=[zeros(n_ele,2*n_ele), eye(n_ele)];
    H=H(ind_sm_sol,:);
    z=lnD_s_meas(ind_s,i-1)-H*x;
    S=H*P*H'+R;
    K=P*H'/S;
    x=x+K*z;
    P=P-K*H*P;
    
    % evidence
    z_post=lnD_s_meas(ind_s,i-1)-H*x;
    logcE_t(i)=-0.5*(z'/S*z+log(det(S))+ns_sol*log(2*pi));
    
    mu_x_t(i,:)=x;
    sig_x_t(i,:)=sqrt(diag(P));
    cov_x_t{i}=P;
    corr_x_t{i}=P./(sig_x_t(i,:)'*sig_x_t(i,:));
    corr_x_t_vec(i,:)=diag(corr_x_t{i}(1:n_ele,n_ele+1:2*n_ele));
    
    x=x(1:2*n_ele);
    P=P(1:2*n_ele,1:2*n_ele);
    
    mu_u_t(i,:)=chol_L\(x-mu_par);
    cov_u_t{i}=chol_L\P/chol_L';
    sig_u_t(i,:)=sqrt(diag(cov_u_t{i}));
    corr_u_t{i}=cov_u_t{i}./(sig_u_t(i,:)'*sig_u_t(i,:));
    corr_u_t_vec(i,:)=diag(corr_u_t{i}(1:n_ele,n_ele+1:2*n_ele));
    
    mu_x_t_lnD_ind(i,:)=A_ind*x;
    sig_x_t_lnD_ind(i,:)=sqrt(diag(A_ind*P*A_ind'));
    
end

logcE_t=cumsum(logcE_t);

figure()
plot(t0,logcE_t)
xlim([0,nt])
xlabel('$t$')
ylabel('$\ln(c_E)$')
title('Model Evidence')

figure()
sgtitle('Filtered $\ln(A)$')
for i=1:ns
    subplot(2,5,i)
    hold on
    plot([0,nt],[lnA_s(i),lnA_s(i)]);
    plot(t0,mu_x_t(:,ind_sm(i)),'k','Linestyle','--')
    plot(t0,norminv(cl,mu_x_t(:,ind_sm(i)),sig_x_t(:,ind_sm(i))),'k')
    plot(t0,norminv(cu,mu_x_t(:,ind_sm(i)),sig_x_t(:,ind_sm(i))),'k')
    xlim([0,nt])
    title(['Location ',num2str(i)])
    xlabel('$t$')
    ylabel('$\ln(A)$')
end

figure()
sgtitle('Filtered $B$')
for i=1:ns
    subplot(2,5,i)
    hold on
    plot([0,nt],[B_s(i),B_s(i)]);
    plot(t0,mu_x_t(:,ind_sm(i)+n_ele),'k','Linestyle','--')
    plot(t0,norminv(cl,mu_x_t(:,ind_sm(i)+n_ele),sig_x_t(:,ind_sm(i)+n_ele)),'k')
    plot(t0,norminv(cu,mu_x_t(:,ind_sm(i)+n_ele),sig_x_t(:,ind_sm(i)+n_ele)),'k')
    xlim([0,nt])
    title(['Location ',num2str(i)])
    xlabel('$t$')
    ylabel('$B$')
end

figure()
sgtitle('Filtered $\ln(D(t))$')
for i=1:ns
    subplot(2,5,i)
    hold on
    plot(t,lnD_s(i,:));
    plot(t,mu_x_t(2:end,ind_sm(i)+2*n_ele),'k','Linestyle','--')
    plot(t,norminv(cl,mu_x_t(2:end,ind_sm(i)+2*n_ele),sig_x_t(2:end,ind_sm(i)+2*n_ele)),'k')
    plot(t,norminv(cu,mu_x_t(2:end,ind_sm(i)+2*n_ele),sig_x_t(2:end,ind_sm(i)+2*n_ele)),'k')
    scatter(t,lnD_s_meas(i,:),5,'filled','k')
    xlim([0,nt])
    title(['Location ',num2str(i)])
    xlabel('$t$')
    ylabel('$\ln(D(t))$')
end


%% posterior field plots

% % plot at time
% t_plot=50;
% ind_plot=t_plot+1;


figure()
subplot(1,3,1)
hold on
plot(xf,lnA_f)
plot(xm,mu_x_t(ind_plot,1:n_ele),'k','Linestyle','--')
plot(xm,norminv(cl,mu_x_t(ind_plot,1:n_ele),sig_x_t(ind_plot,1:n_ele)),'k')
plot(xm,norminv(cu,mu_x_t(ind_plot,1:n_ele),sig_x_t(ind_plot,1:n_ele)),'k')
plot([0,L],[norminv(cl,mu_lnA,sig_lnA),norminv(cl,mu_lnA,sig_lnA)],'k','Linestyle','-.')
plot([0,L],[norminv(cu,mu_lnA,sig_lnA),norminv(cu,mu_lnA,sig_lnA)],'k','Linestyle','-.')
axis([0 L mu_lnA-3*sig_lnA mu_lnA+3*sig_lnA])
title(['Field $\ln(A)$ at $t=',num2str(t_plot),'$'])

subplot(1,3,2)
hold on
plot(xf,B_f)
plot(xm,mu_x_t(ind_plot,n_ele+1:2*n_ele),'k','Linestyle','--')
plot(xm,norminv(cl,mu_x_t(ind_plot,n_ele+1:2*n_ele),sig_x_t(ind_plot,n_ele+1:2*n_ele)),'k')
plot(xm,norminv(cu,mu_x_t(ind_plot,n_ele+1:2*n_ele),sig_x_t(ind_plot,n_ele+1:2*n_ele)),'k')
plot([0,L],[norminv(cl,mu_B,std_B),norminv(cl,mu_B,std_B)],'k','Linestyle','-.')
plot([0,L],[norminv(cu,mu_B,std_B),norminv(cu,mu_B,std_B)],'k','Linestyle','-.')
axis([0 L mu_B-3*std_B mu_B+3*std_B])
title(['Field $B$ at $t=',num2str(t_plot),'$'])

subplot(1,3,3)
hold on
plot(xf,lnD_f(:,t_plot))
plot(xm,mu_x_t(ind_plot,2*n_ele+1:end),'k','Linestyle','--')
plot(xm,norminv(cl,mu_x_t(ind_plot,2*n_ele+1:end),sig_x_t(ind_plot,2*n_ele+1:end)),'k')
plot(xm,norminv(cu,mu_x_t(ind_plot,2*n_ele+1:end),sig_x_t(ind_plot,2*n_ele+1:end)),'k')
plot([0,L],[norminv(cl,mu_lnD(t_plot),sig_lnD(t_plot)),norminv(cl,mu_lnD(t_plot),sig_lnD(t_plot))],'k','Linestyle','-.')
plot([0,L],[norminv(cu,mu_lnD(t_plot),sig_lnD(t_plot)),norminv(cu,mu_lnD(t_plot),sig_lnD(t_plot))],'k','Linestyle','-.')
axis([0 L mu_lnD(t_plot)-3*sig_lnD(t_plot) mu_lnD(t_plot)+3*sig_lnD(t_plot)])
title(['Field $\ln(D(t=',num2str(t_plot),'))$ at $t=',num2str(t_plot),'$'])

