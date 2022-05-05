function [mu_x,std_x]=w_mom(x,w)

w=w./sum(w);

if iscolumn(w)
    w=w';
end

mu_x=w*x;
n=length(mu_x);
std_x=zeros(1,n);
for i=1:n
    std_x(i)=sqrt(w*(x(:,i)-mu_x(i)).^2);
end

end