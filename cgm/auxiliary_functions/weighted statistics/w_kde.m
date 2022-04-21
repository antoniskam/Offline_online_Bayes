function f_eval=w_kde(x,w,x_eval)
%{
Kernel density estimation of the pdf with weighted samples.

Input:
x           samples (vector)
w           sample weights (vector)
x_eval      density evaluation points (vector)

Output:
f_eval      probability density evaluations
%}

% adjustment of input
if isrow(x)
    x=x';
end
if isrow(w)
    w=w';
end

[x,i]=sort(x);
w=w(i);
w=w/sum(w);

% interquartile range R
wsum=cumsum(w);
i25=find(wsum>0.25,1);
i75=find(wsum>0.75,1);
x25=x(i25-1)+(0.25-wsum(i25-1))/(wsum(i25)-wsum(i25-1))*(x(i25)-x(i25-1));
x75=x(i75-1)+(0.75-wsum(i75-1))/(wsum(i75)-wsum(i75-1))*(x(i75)-x(i75-1));
R=x75-x25;

% moments
mu_x=x'*w;
sig_x=sqrt(w'*(x-mu_x).^2);

% sample size
% N=1/sum(w.^2); % effective sample size
N=length(x);

% optimal bandwidth
h=1.06*min([sig_x,R])*N^-0.2;

% density evaluation
nx=length(x_eval);
f_eval=zeros(size(x_eval));
for i=1:nx
    f_eval(i)=w'*normpdf((x_eval(i)-x)/h)/h;
end

end