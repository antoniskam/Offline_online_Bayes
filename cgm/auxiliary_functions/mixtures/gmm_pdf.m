function pdf=gmm_pdf(X, mu, si, pi)

n_GM=length(pi);
pdf=pi(1)*mvnpdf(X,mu(:,1)',si(:,:,1));
for i=2:n_GM
    pdf=pdf+pi(i)*mvnpdf(X,mu(:,i)',si(:,:,i));
end

end