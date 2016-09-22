%% Input:
%      data         A row vector of observations.
%      alpha        A scalar. The concentration parameter for the Dirichlet
%                   process.
%      a            A scalar. The shape for the Gamma distribution.
%      b            A scalar. The rate for the Gamma distribution.
%      T            A scalar. The truncation value of the infinite vector.
%      maxIter      A scalar. The maximum number of iterations.
%% Output:
%      mu           A struct. The first field contains the
%                   means of the means of each atoms and the second field 
%                   contains the precisions of the mean of each atom.
%      lambda       A row struct. The first field contains the
%                   shapes of lambda, and the second field contain the 
%                   rates of lambda.
%      pi           A struct. The two fields contain the two parameters of
%                   pi.
%      z            A matrix of size length(data) * T. The variational
%                   parameters for the indicators of each observation.

function [z, mu, lambda, pi] = vb(data, alpha, a, b, T, maxIter)
if nargin < 2
    alpha = 1;
end
if nargin < 3
    a = 1;
end
if nargin < 4
    b = 1;
end
if nargin < 5
    T = 20;
end
if nargin < 6
    maxIter = 100;
end

%% Init: 

% each row is corresponding to a observation and each column is
% corresponding to a cluster
z = ones(length(data), T);
z = z / T;

% the parameters for the cluster means (normal distribution)
mu = repmat(struct('mean', mean(data), 'precision', 1), T, 1);

% the parameters for the cluster precisions (Gamma distribution)
lambda = repmat(struct('a', 1, 'b', 1), T, 1);

% the variational parameters for pi (Beta distribution)
pi = repmat(struct('a1', 1, 'a2', 1), T, 1);

% record the mean for every iteration
mean_vec = zeros(1, maxIter) + mean(data);

%%
for iter = 1:maxIter
    
    tic
    
    % two useful vectors
    z_agg1 = sum(z, 1);
    z_agg2 = cumsum([z_agg1, 0], 'reverse');
    
    % parameters of pi
    for j = 1:T
        pi(j).a1 = 1 + z_agg1(j);
        pi(j).a2 = alpha + z_agg2(j+1);
    end
    
    % parameters of mu
    for j = 1:T
        mu(j).mean = data * z(:,j) / (1 + z_agg1(j));
        mu(j).precision = (1 + z_agg1(j)) * lambda(j).a / lambda(j).b;
    end
    mean_vec(iter) = mu(1).mean;
    
    % parameters of lambda
    for j = 1:T
        lambda(j).a = a + 1/2 + 1/2 * z_agg1(j);
        E_squre = mu(j).mean^2 + 1 / mu(j).precision;
        lambda(j).b = b + 1/2 * E_squre + 1/2 * (data .^ 2) * z(:,j)...
            + 1/2 * z_agg1(j) * E_squre - data * z(:,j) * mu(j).mean;
    end
    
    % parameters of indicators
    for i = 1:length(data)
        for j = 1:T
            first = psi(pi(j).a1) - psi(pi(j).a1 + pi(j).a2);
            
            second = 0;
            if j > 1
                for t = 1:j-1
                    second = second + psi(pi(j).a2) - psi(pi(j).a1 + pi(j).a2);
                end
            end
            
            third = psi(lambda(j).a) - log(lambda(j).b);
            
            fourth = 1/2 * lambda(j).a / lambda(j).b * ...
                ( data(i)^2 + mu(j).mean^2 - 1 / mu(j).precision - 2 * data(i) * mu(j).mean);
            
            z(i,j) = exp(first + second + third - fourth);
        end
    z = z ./ repmat(sum(z, 2), 1, T);    
    end
    
    fprintf(['iter ', num2str(iter), ' done \n'])
    
    toc
end

plot(mean_vec)
    
end







