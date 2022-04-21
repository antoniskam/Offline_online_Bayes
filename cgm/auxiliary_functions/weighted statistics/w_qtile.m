function q_val=w_qtile(x,w,q)
%{
Comqutation of quantile values for set of weighted samqles.

Inqut:
x           samples 
w           sample weights (vector)
q           quantile

Outqut:
q_eval      quantile values (single value or vector)
%}

[~,d]=size(x);
q_val=zeros(1,d);
w=w/sum(w);

for i=1:d
    xcur=x(:,i);
    [xcur,ind]=sort(xcur);
    w_s=w(ind);
    cum_w=cumsum(w_s)-w_s/2;

    ind_q=find(cum_w>q,1);

    if ind_q==1
        q_val(i)=xcur(1);
    elseif isempty(ind_q)
        q_val(i)=xcur(end);
    else
        q2=cum_w(ind_q);
        d2=xcur(ind_q);
        q1=cum_w(ind_q-1);
        d1=xcur(ind_q-1);

        q_val(i)=d1+(q-q1)/(q2-q1)*(d2-d1);
    end
end

end

