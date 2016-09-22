%% Input:
%     K         The number of clusters
%     n         The number of observations in each cluster
%% Output:
%     data      A row vector of length K * n. 
function [data, z] = data_generate(K, n)
% generate probabilities for each cluster
p = rand(1, K);
p = p / sum(p);
p = [0, cumsum(p)];

% generate mean
mu = 1 : 5 : 5*K;
mu = mu - mean(mu);

% generate variance
theta = gamrnd(ones(1,K), 1/2);
% get standard deviation
theta = sqrt(theta);

% sample a cluster for each observation
[~, ~, z] = histcounts(rand(1, K*n), p);

% generate data
data = zeros(1, n*K);
for i = 1:length(data)
    data(i) = randn(1) * theta(z(i)) + mu(z(i));
end

end


