function out=alg_SOCP(net_original, net_noise) 
%% INITIALIZATION

K = net_noise.K;
m = net_noise.anchors;
N = K - m;
n = size(net_noise.Matrices.X_real, 1);
F = net_noise.GI.F;
bias = net_noise.GI.bias;
amatrix = net_noise.Matrices.noised_anchors;
edges = net_original.net.GI.edges;
num_edges = length(edges);
noised_distnaces = net_noise.Matrices.noised_distances;

%% ALGORITHM
tic
cvx_begin
variable X(n, K)
variables v s(m) t(num_edges) q(num_edges)
minimize v
subject to
norm([t; s]) <= v;
for i = 1:num_edges
    abs(q(i) - noised_distnaces(edges(i, 1), edges(i, 2))) <= t(i);
    norm(X(:, edges(i, 1)) - X(:, edges(i, 2))) <= q(i);
end
for i=1:m
    norm(X(:, N + i) - amatrix(:, i)) <= s(i);
end
cvx_end
time_end = toc;

%% OUTPUT

out = struct;
out.location = X;
out.time = time_end;
out.fun_val = F(X)
out.norm_bias = bias(X)
disp('SOCP is Done!')