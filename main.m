clear
clc

% Generate and exhibite dataset
[data, z_true] = data_generate(3, 300);
tabulate(z_true)

% Inference
[z, mu, lambda, pi] = vb(data);

% Post processing
[~, z] = max(z, [], 2);

% The number of observations in each cluster
tabulate(z)
figure(2)
hold on
for i = 1:max(z)
    plot(data(z == i), ones(1, sum(z == i)) * i, 'o')
end
hold off