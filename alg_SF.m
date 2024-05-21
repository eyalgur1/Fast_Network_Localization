function out = alg_SF(net_original, net_noise, x0, max_iter)
%% INITIALIZATION

IMnet = net_noise.Matrices.IM_net;
IMDnet = net_noise.Matrices.IMD_net;
K = net_noise.K;
m = net_noise.anchors;
N = K - m;
n = size(net_noise.Matrices.X_real, 1);
F = net_noise.GI.F;
bias = net_noise.GI.bias;
amatrix = net_noise.Matrices.noised_anchors;
a = reshape(amatrix, n*m, 1);
L = net_original.net.GI.Lip_dist - 1;
nonanchor_indices = 1:N;
Laplacian = sparse(IMDnet'*IMDnet);  % Laplacian matrix of the entire network with distances
QTQtilde = sparse(IMnet'*IMnet);  % Laplacian matrix of the entire network

E = {};
for i = 1:K % sets the sets of neighbors
    E{end+1} = find(QTQtilde(i, :) == -1);  %finding all the neigbors in one go
end

bias_v = zeros(max_iter, 1); fun_v = bias_v; t = bias_v; tp = bias_v;
x0 = reshape(x0, n, K); 
x0(:, N+1:end) = amatrix;
Fx0 = F(x0); biasx0 = bias(x0);
x0 = x0(:, 1:N);
xk = [reshape(x0, n*N, 1); a];
xk_prev = xk;
iter = 0;

%% Algorithm
fprintf('SF | Iter=%2d | Bias=%10.10f | F=%10.10f\n', 0, biasx0, Fx0)  

while iter < max_iter
    starta = tic;

    w = xk + ((iter-2)/(iter+1))*(xk - xk_prev);  % Nesterov's optimal gradient method
    grad = zeros(K*n,1);
    f_v = 0;
    enda = toc(starta);  % for parallel time

    endi = zeros(1, length(nonanchor_indices));  % time for each sensor i at each iteration
    for i = nonanchor_indices  % considering non-anchors only
        starti = tic;  % for parallel time
        E_i = [E{i}];  % neighbor set of non-anchor i
        i_index = (i-1)*n+1:i*n;
        for j = E_i
            j_index = (j-1)*n + 1:j*n;
            absLap = abs(Laplacian(i, j));
            sqrtabsLap = sqrt(absLap);
            normdiff = norm(w(i_index) - w(j_index));
            if normdiff > sqrtabsLap
                grad(i_index) = grad(i_index) + (w(i_index) - w(j_index))*(1 - sqrtabsLap/normdiff); %and easy calculation of the gradient
                f_v = f_v + 0.5*(normdiff - absLap);
            end
        end
        endi(i) = toc(starti);  % for parallel time
    end

    % update the algorithm
    startb = tic;  % for parallel time
    xk_prev = xk;
    xk = w - (1/L)*grad;
    iter = iter+1;
    endb = toc(startb);  % for parallel time

    t(iter) = toc(starta);  % for non-parallel time
    tp(iter) = max(endi) + (enda+endb)/K;  % for parallel time
    xk_matrix = reshape(xk, n, K); 
    bias_v(iter) = bias(xk_matrix);
    fun_v(iter) = F(xk_matrix);

    if ~mod(iter, 100) || iter == 1
        fprintf('SF | Iter=%2d | Bias=%10.10f | F=%10.10f\n', iter, bias_v(iter), fun_v(iter))
    end
end


%% OUTPUT
out = struct;
out.location = xk;
out.time = sum(t);
out.parllel_time = sum(tp);
out.fun_val = [Fx0; fun_v];
out.norm_bias = [biasx0; bias_v];
disp('SF is done!')