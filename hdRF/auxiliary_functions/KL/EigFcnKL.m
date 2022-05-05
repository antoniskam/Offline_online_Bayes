function [lambda , phi] = EigFcnKL(m , L , l)

% Returns the eigenvalues and eigenfunctions of the KL expansion with
% exponential kernel
% m: number of random variables
% L: Length of the domain
% l: correlation length
% lambda: vector of eigenvalues
% phi: cell array of eigenfunctions

lambda = zeros(m,1);
phi = cell(m,1);

for i = 1:m
    if (mod(i,2) == 0) % even
        lowerBound = (i-0.999)*pi/L;
        upperBound = 0.999*i*pi/L;
        w = fzero(@(w)1/l*tan(w*L/2)+w,[lowerBound,upperBound]);
        
        lambda(i) = 2*l / (1 + w^2*l^2);
        
        alpha = 1 / sqrt(L/2 - sin(w*L)/(2*w));
        phi{i} = @(z)alpha * sin(w*(z - L/2));
        
    else % odd
        lowerBound = (i-0.999)*pi/L;
        upperBound = 0.999*i*pi/L;
        w = fzero(@(w)1/l-w*tan(w*L/2),[lowerBound,upperBound]);
        
        lambda(i) = 2*l / (1 + w^2*l^2);
        
        alpha = 1 / sqrt(L/2 + sin(w*L)/(2*w));
        phi{i} = @(z)alpha * cos(w*(z - L/2));
        
    end
end

end