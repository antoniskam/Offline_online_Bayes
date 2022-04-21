function X=gmm_rand(n, mu, si, pi)

n_GM=length(pi);
n_dim=size(mu,1);
ind = randsample(n_GM, n, true,pi);
n_ind=hist(ind,1:n_GM); %#ok<HIST>
n_ind_cum=[0,cumsum(n_ind)];
X=zeros(n,n_dim);
for i=1:n_GM
    X(n_ind_cum(i)+1:n_ind_cum(i+1),:)=mvnrnd(mu(:,i),si(:,:,i),n_ind(i));
end
ind=randperm(n,n);
X=X(ind,:);

end