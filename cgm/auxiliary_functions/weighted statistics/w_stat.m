function [mu_x,cov_x,std_x,corr_x]=w_stat(x,w)

w=w./sum(w);

if iscolumn(w)
    w=w';
end

mu_x=w*x;
cov_x=((x-mu_x).*w')'*(x-mu_x);
std_x=sqrt(diag(cov_x));
corr_x=cov_x./(std_x*std_x');

end