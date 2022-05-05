clear; clc; close all;
rng(1)
set(0,'DefaultTextInterpreter','latex')

% Import the KF reference solution Matlab data: 'rng_numelements_numsensors_s'
name='rng40_25_10_s.mat';
load(name)

%% filter

N=2000;
nGM=8;
N_T=0.1*N;

corr_prior=[Rho,zeros(n_ele);zeros(n_ele),Rho];
marg=[repmat(ERADist('normal','MOM',[mu_lnA,sig_lnA]),n_ele,1);repmat(ERADist('normal','MOM',[mu_B,std_B]),n_ele,1)];
prior=ERANataf(marg,corr_prior);

part_x=prior.random(N);
part_u=prior.X2U(part_x);
w=1/N*ones(N,1);
logw=log(w);

% pre-allocation of parameter estimate statistics in physical space
x_p_mu=zeros(nt+1,n_ele*2);
x_p_mu(1,:)=mean(part_x);
x_p_std=zeros(nt+1,n_ele*2);
x_p_std(1,:)=std(part_x);
x_p_corr=cell(nt+1,1);
x_p_corr{1}=corr(part_x);
x_p_corr_vec=zeros(nt+1,n_ele);
x_p_corr_vec(1,:)=diag(x_p_corr{1}(1:n_ele,n_ele+1:2*n_ele));
x_p_median=zeros(nt+1,n_ele*2);
x_p_median(1,:)=median(part_x);
x_p_90_l=zeros(nt+1,n_ele*2);
x_p_90_l(1,:)=quantile(part_x,0.05);
x_p_90_u=zeros(nt+1,n_ele*2);
x_p_90_u(1,:)=quantile(part_x,0.95);

% pre-allocation of parameter estimate statistics in standard normal space
u_p_mu=zeros(nt+1,n_ele*2);
u_p_mu(1,:)=mean(part_u);
u_p_std=zeros(nt+1,n_ele*2);
u_p_std(1,:)=std(part_u);
u_p_corr=cell(nt+1,1);
u_p_corr{1}=corr(part_u);
u_p_corr_vec=zeros(nt+1,n_ele);
u_p_corr_vec(1,:)=diag(u_p_corr{1}(1:n_ele,n_ele+1:2*n_ele));
u_p_median=zeros(nt+1,n_ele*2);
u_p_median(1,:)=median(part_u);
u_p_90_l=zeros(nt+1,n_ele*2);
u_p_90_l(1,:)=quantile(part_u,0.05);
u_p_90_u=zeros(nt+1,n_ele*2);
u_p_90_u(1,:)=quantile(part_u,0.95);

x_p_corr_vec_comp=zeros(nt+1,(4*n_ele^2-2*n_ele)/2);
u_p_corr_vec_comp=zeros(nt+1,(4*n_ele^2-2*n_ele)/2);
x_p_corr_vec_comp_cur=[];
u_p_corr_vec_comp_cur=[];
for j=1:2*n_ele-1
    x_p_corr_vec_comp_cur=[x_p_corr_vec_comp_cur,x_p_corr{1}(j,j+1:end)]; %#ok<AGROW>
    u_p_corr_vec_comp_cur=[u_p_corr_vec_comp_cur,u_p_corr{1}(j,j+1:end)]; %#ok<AGROW>
end
x_p_corr_vec_comp(1,:)=x_p_corr_vec_comp_cur;
u_p_corr_vec_comp(1,:)=u_p_corr_vec_comp_cur;

% pre-allocation of state estimate statistics (physical & log space)
lnD_p_mu=zeros(nt,n_ele);
lnD_p_std=zeros(nt,n_ele);
lnD_p_corr=cell(nt,1);
lnD_p_median=zeros(nt,n_ele);
lnD_p_90_l=zeros(nt,n_ele);
lnD_p_90_u=zeros(nt,n_ele);

N_eff_seq=zeros(nt+1,1);  % effective sample size (based on weights)
N_eff_seq(1)=N;
N_part=N_eff_seq;       % number of different particles

logcE_p=0;
logcE_pt=zeros(nt+1,1);
res_i=[];

count=ones(nt,1);
temp=zeros(nt,10);
N_eff_t=temp;

for i=2:nt+1
    
    lnA_part = part_x(:,1:n_ele);
    B_part = part_x(:,n_ele+1:end);
    lnD_cur = lnA_part+B_part*log(t(i-1));
    
    % likelihood evaluation
    log_lh_cur=loglikelihood(part_x, ind_sm_sol, std_meas, t(i-1), lnD_s_meas(ind_s,i-1));
    q=0;
    
    while q<1
        
        logw_cur=logw+(1-q)*log_lh_cur;
        logcE_p_t_cur=logsumexp(logw_cur);
        logcE_p_cur=logcE_p+logcE_p_t_cur;
        logw_cur=logw_cur-logcE_p_t_cur;
        w_cur=exp(logw_cur);
        
        % effective sample size
        N_eff_cur=1/sum(w_cur.^2);
        
        if N_eff_cur < N_T
            res_i=[res_i,i]; %#ok<AGROW>
            
            fun = @(dq) exp(2*logsumexp(logw+abs(dq)*log_lh_cur)-logsumexp(2*(logw+abs(dq)*log_lh_cur))) - N_T; 
            [dq,~,flag] = fzero(fun, 0);
            dq=abs(dq);
            temp(i-1,count(i-1))=q+dq;
            q_new = min(1, q+dq);
            dq=q_new-q;
            q=q_new;
            
            logw=logw+dq*log_lh_cur;
            logcE_p_t=logsumexp(logw);
            logcE_p=logcE_p+logcE_p_t;
            logw=logw-logcE_p_t;
            w=exp(logw);
            N_eff_t(i-1,count(i-1))=exp(2*logsumexp(logw)-logsumexp(2*logw));
            
            [mu, si, p] = EMGM(part_u',w,nGM);
            part_u = gmm_rand(N, mu, si, p);
            part_x = prior.U2X(part_u);
            if q<1
                log_lh_cur = loglikelihood(part_x, ind_sm_sol, std_meas, t(i-1), lnD_s_meas(ind_s,i-1));
            end
            
            w=1/N*ones(N,1); % reweighting
            logw=log(w);
            count(i-1)=count(i-1)+1;
        else
            q=1;
            logcE_p=logcE_p_cur;
            logw=logw_cur;
            w=w_cur;
            N_eff_t(i-1,count(i-1))=N_eff_cur;
            temp(i-1,count(i-1))=1;
        end
    end
    logcE_pt(i)=logcE_p;
    N_eff_seq(i)=N_eff_cur;
    
    % parameter estimate statistics in physical space
    [mu_x_p,cov_x_p,std_x_p,corr_x_p]=w_stat(part_x,w);
    x_p_mu(i,:)= mu_x_p;
    x_p_std(i,:)= std_x_p;
    x_p_corr{i}= corr_x_p;
    x_p_corr_vec(i,:)=diag(x_p_corr{i}(1:n_ele,n_ele+1:2*n_ele));
    x_p_median(i,:)=w_qtile(part_x,w,0.5);
    x_p_90_l(i,:)=w_qtile(part_x,w,0.05);
    x_p_90_u(i,:)=w_qtile(part_x,w,0.95);
    
    % parameter estimate statistics in standard normal space
    [mu_u_p,cov_u_p,std_u_p,corr_u_p]=w_stat(part_u,w);
    u_p_mu(i,:)= mu_u_p;
    u_p_std(i,:)= std_u_p;
    u_p_corr{i}= corr_u_p;
    u_p_corr_vec(i,:)=diag(u_p_corr{i}(1:n_ele,n_ele+1:2*n_ele));
    u_p_median(i,:)=w_qtile(part_u,w,0.5);
    u_p_90_l(i,:)=w_qtile(part_u,w,0.05);
    u_p_90_u(i,:)=w_qtile(part_u,w,0.95);
    
    x_p_corr_vec_comp_cur=[];
    u_p_corr_vec_comp_cur=[];
    for j=1:2*n_ele-1
        x_p_corr_vec_comp_cur=[x_p_corr_vec_comp_cur,x_p_corr{i}(j,j+1:end)]; %#ok<AGROW>
        u_p_corr_vec_comp_cur=[u_p_corr_vec_comp_cur,u_p_corr{i}(j,j+1:end)]; %#ok<AGROW>
    end
    x_p_corr_vec_comp(i,:)=x_p_corr_vec_comp_cur;
    u_p_corr_vec_comp(i,:)=u_p_corr_vec_comp_cur;
    
    % state estimate statistics (physical & log space)
    [mu_lnD_p,~,std_lnD_p,corr_lnD_p]=w_stat(lnD_cur,w);
    lnD_p_mu(i-1,:)= mu_lnD_p;
    lnD_p_std(i-1,:)= std_lnD_p;
    lnD_p_corr{i-1}=corr_lnD_p;
    lnD_p_median(i-1,:)=w_qtile(lnD_cur,w,0.5);
    lnD_p_90_l(i-1,:)=w_qtile(lnD_cur,w,0.05);
    lnD_p_90_u(i-1,:)=w_qtile(lnD_cur,w,0.95);
end

model_eval=sum(count)*N;

%% figures

figure()
hold on
plot(t0,N_eff_seq,'k')
plot([0,nt],[N_T,N_T],'r')
xlim([0,nt])
ylim([0,N*1.1])
xlabel('$t$')
ylabel('$N_{Eff}$')
title('Effective Sample Size')

figure()
hold on
plot(t0,logcE_t,'k')
plot(t0,logcE_pt,'r')
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
    plot(t0,x_p_mu(:,ind_sm(i)),'r','Linestyle','--')
    plot(t0,x_p_90_l(:,ind_sm(i)),'r')
    plot(t0,x_p_90_u(:,ind_sm(i)),'r')
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
    plot(t0,x_p_mu(:,ind_sm(i)+n_ele),'r','Linestyle','--')
    plot(t0,x_p_90_l(:,ind_sm(i)+n_ele),'r')
    plot(t0,x_p_90_u(:,ind_sm(i)+n_ele),'r')
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
    plot(t,lnD_p_mu(:,ind_sm(i)),'r','Linestyle','--')
    plot(t,lnD_p_90_l(:,ind_sm(i)),'r')
    plot(t,lnD_p_90_u(:,ind_sm(i)),'r')
    xlim([0,nt])
    title(['Location ',num2str(i)])
    xlabel('$t$')
    ylabel('$\ln(D(t))$')
end


%% posterior field plots

% % plot at time
% t_plot=50;
% ind_plot=t_plot+1;

X_sol=repmat(xm,n_ele,1);

figure()
subplot(1,3,1)
hold on
plot(xf,lnA_f)
plot(xm,mu_x_t(ind_plot,1:n_ele),'k','Linestyle','--')
plot(xm,norminv(cl,mu_x_t(ind_plot,1:n_ele),sig_x_t(ind_plot,1:n_ele)),'k')
plot(xm,norminv(cu,mu_x_t(ind_plot,1:n_ele),sig_x_t(ind_plot,1:n_ele)),'k')
plot(xm,x_p_mu(ind_plot,1:n_ele),'r','Linestyle','--')
plot(xm,x_p_90_l(ind_plot,1:n_ele),'r')
plot(xm,x_p_90_u(ind_plot,1:n_ele),'r')
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
plot(xm,x_p_mu(ind_plot,n_ele+1:end),'r','Linestyle','--')
plot(xm,x_p_90_l(ind_plot,n_ele+1:end),'r')
plot(xm,x_p_90_u(ind_plot,n_ele+1:end),'r')
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
plot(xm,lnD_p_mu(ind_plot-1,:),'r','Linestyle','--')
plot(xm,lnD_p_90_l(ind_plot-1,:),'r')
plot(xm,lnD_p_90_u(ind_plot-1,:),'r')
plot([0,L],[norminv(cl,mu_lnD(t_plot),sig_lnD(t_plot)),norminv(cl,mu_lnD(t_plot),sig_lnD(t_plot))],'k','Linestyle','-.')
plot([0,L],[norminv(cu,mu_lnD(t_plot),sig_lnD(t_plot)),norminv(cu,mu_lnD(t_plot),sig_lnD(t_plot))],'k','Linestyle','-.')
axis([0 L mu_lnD(t_plot)-3*sig_lnD(t_plot) mu_lnD(t_plot)+3*sig_lnD(t_plot)])
title(['Field $\ln(D(t=',num2str(t_plot),'))$ at $t=',num2str(t_plot),'$'])


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lh=loglikelihood(xi, ind, std_meas, t, meas)
n=size(xi,2);
lnA = xi(:,ind);
B = xi(:,ind+n/2);
lnD = lnA+B*log(t);
lh=-0.5*((lnD-meas')/std_meas).^2+log(1/(sqrt(2*pi)*std_meas));
lh = sum(lh,2);

end

function s = logsumexp(x, dim)
% Compute log(sum(exp(x),dim)) while avoiding numerical underflow.
% By default dim = 1 (columns).
% Written by Michael Chen (sth4nth@gmail.com).

if nargin == 1
    % Determine which dimension sum will use
    dim = find(size(x)~=1,1);
    if isempty(dim)
        dim = 1;
    end
end

% subtract the largest in each column
y = max(x,[],dim);
x = bsxfun(@minus,x,y);
s = y + log(sum(exp(x),dim));
i = find(~isfinite(y));
if ~isempty(i)
    s(i) = y(i);
end

end
