function out=alg_SGO(net_original, net_noise, x0, max_iter)
%% INITIALIZATION

%IMnet = net_noise.Matrices.IM_net;
K = net_noise.K;
m = net_noise.anchors;
N = K - m;
n = size(net_noise.Matrices.X_real, 1);
amatrix = net_noise.Matrices.noised_anchors;
F = net_noise.GI.F;
bias = net_noise.GI.bias;
ball_radius = net_noise.GI.ball_radius;  % radius of ball uncertainty (0 or greater)
edges = net_original.net.GI.edges;
Ne = length(edges(:,1));

bias_v = zeros(max_iter, 1); fun_v = bias_v; t = bias_v; tp = bias_v;
x0 = reshape(x0, n, K);
x0(:, N+1:end) = amatrix;
xk = x0;
iter = 0;

%% ALGORITHM
fprintf('SGO | Iter=%2d | Bias=%10.10f | F=%10.10f\n', 0, bias(x0), F(x0))  

while iter < max_iter
    starta=tic;

    for i = 1:N
        neigh_i = net_noise.node{i}.neighbors;
        xk_neigh_i = xk(:, neigh_i);
        dist_i = full(net_noise.node{i}.distance);
        diff_ij = xk(:,i) - xk_neigh_i;
        norm_i = sqrt(sum(diff_ij.^2));
        v_j = sum(xk_neigh_i + dist_i.*(diff_ij./norm_i), 2);
        xk(:,i) = (1/length(neigh_i))*v_j;
    end
    for i = N+1:K
        neigh_i = net_noise.node{i}.neighbors;
        xk_neigh_i = xk(:,neigh_i);
        dist_i = full(net_noise.node{i}.distance);
        diff_ij = xk(:,i) - xk_neigh_i;
        norm_i = sqrt(sum(diff_ij.^2));
        v_j = sum(xk_neigh_i + dist_i.*(diff_ij./norm_i), 2);
        amatrix_i = amatrix(:,i-N);
        x_i_hat = (1/(length(neigh_i)+1))*(v_j + amatrix_i);
        diff_i = x_i_hat - amatrix_i;
        xk(:,i) = amatrix_i + ball_radius*(diff_i/max(norm(diff_i), ball_radius));
    end

    iter = iter + 1;

    t(iter) = toc(starta); % upper bound on the non-parallelizable running time
    tp(iter) = t(iter); % running time estimation
    bias_v(iter) = bias(xk);
    fun_v(iter) = F(xk);
    if ~mod(iter, 1000) || iter == 1
        fprintf('SGO | Iter=%2d | Bias=%10.10f | F=%10.10f\n', iter, bias_v(iter), fun_v(iter))
    end
end

%% OUTPUT

out = struct;
out.location = xk;
out.time = sum(t);
out.parllel_time = sum(tp);
out.fun_val = [F(x0); fun_v];
out.norm_bias = [bias(x0); bias_v];
disp('SGO is Done!')