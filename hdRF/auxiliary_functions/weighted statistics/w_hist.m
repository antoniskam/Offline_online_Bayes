function w_hist(x,w,nbins)
%{
Creation of a PDF normalized histogram of weighted samples.

Input:
x           samples (vector)
w           sample weights (vector)
nbins       number of bins (single value)

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

bin_width=(x(end)-x(1))/nbins;
edges=linspace(x(1),x(end),nbins+1);
counts=zeros(nbins,1);
for i=1:nbins
    ind=edges(i)<x & x<=edges(i+1);
    counts(i)=sum(w(ind));
end
counts=counts/bin_width;

% plot histogram
histogram('BinEdges',edges,'BinCounts',counts)

end