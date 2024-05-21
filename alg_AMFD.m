function out = alg_AMFD(net_noise, x0, max_iter)
%% INITIALIZATION

K = net_noise.K;
m = net_noise.anchors;
N = K - m;
n = size(net_noise.Matrices.X_real, 1);
Atilde = net_noise.Matrices.A_tilde;
Qtilde = net_noise.Matrices.Q_tilde;
Btilde = net_noise.Matrices.B_tilde;
ADtilde = net_noise.Matrices.AD_tilde;
QDtilde = net_noise.Matrices.QD_tilde;
BDtilde = net_noise.Matrices.BD_tilde;
Ptilde = net_noise.Matrices.P_tilde;
amatrix = net_noise.Matrices.noised_anchors;
F = net_noise.GI.F;
bias = net_noise.GI.bias;
v_x_coef = [QDtilde; ADtilde]';
v_const = [zeros(2, size(QDtilde,1)), -amatrix*BDtilde'];
Ne = size(Qtilde, 1) + size(Atilde, 1);
Neq = size(Qtilde, 1);

bias_v = zeros(max_iter, 1); fun_v = bias_v; t = bias_v; tp = bias_v;
x0 = reshape(x0, n, K);
x0(:, N+1:end) = amatrix;  % set amatrix as the initial point for the anchors
xk = x0(:, 1:N);
uk = zeros(n, Ne);
iter = 0;

% generating N clusters (each sensor is a cluster)
clusters_strct = MakeClusters(net_noise);
clusters = clusters_strct.clusters;

% setting coefficients according to the clusters
for i=1:N
    S_i=clusters{i}; % cluster S_i
    clusters{i}.Ptilde = Ptilde(S_i.nodes, S_i.nodes);
    clusters{i}.invPtilde = inv(clusters{i}.Ptilde);
    invPtilde = clusters{i}.invPtilde;
    clusters{i}.x_const = amatrix(:, S_i.anchor_neigbhors)*(invPtilde*Atilde(S_i.anchor_edges, S_i.nodes)'*Btilde(S_i.anchor_edges, S_i.anchor_neigbhors))'; % linear coeeficient of x_i: invP*A'*B*a
    clusters{i}.x_u_coef = (invPtilde*[QDtilde(S_i.nonanchor_edges, S_i.nodes)', ADtilde(S_i.anchor_edges, S_i.nodes)'])'; % quadratic coefficient of x_i.u_i: invP*[QD' AD']
    clusters{i}.x_x_coef = -(invPtilde*Ptilde(S_i.nonanchor_neigbhors, S_i.nodes)')'; % quadratic coefficient of x_i.x_i: -invP
end

%% ALGORITHM
fprintf('AMFD | Iter=%2d | Bias=%10.10f | F=%10.10f\n', 0, bias(x0), F(x0))  

while iter < max_iter
    starta = tic;
    iter=iter+1;

    % x update and parallel run time calculations
    endi = zeros(1, N); % setting run time of each cluster
    for i = 1:N % update the x block
        starti = tic; % run time for each cluster
        S_i = clusters{i};
        rel_edges = [S_i.nonanchor_edges; S_i.anchor_edges + Neq]; % cluster edges
        if ~isempty(S_i.nonanchor_neigbhors) % if the cluster has non-anchor neighbors
            xk(:, S_i.nodes) = xk(:,S_i.nonanchor_neigbhors)*S_i.x_x_coef + uk(:,rel_edges)*S_i.x_u_coef + S_i.x_const; % update x of each cluster
        else % if the cluster has only anchor neighbors
            xk(:, S_i.nodes) = uk(:, rel_edges)*S_i.x_u_coef + S_i.x_const; % update x of each cluster
        end
        endi(i) = toc(starti); % run time of the cluster (sensor)
    end

    [uk, u_time, ~] = u_update(xk, uk, v_x_coef, v_const, Ne, n, K, {}); 
    u_max_time = max(u_time(1:N)); 

    % update outputs
    t(iter) = toc(starta);
    bias_v(iter) = bias([xk amatrix]);
    fun_v(iter) = F([xk amatrix]);
    %anc_AMFD = amatrix(:, 1)
    tp(iter) = sum(endi) + u_max_time;  % parallel run time

    if ~mod(iter,100) || iter==1
        fprintf('AMFD | Iter=%2d | Bias=%10.10f | F=%10.10f\n', iter, bias_v(iter), F([xk amatrix]))
    end
end

%% OUTPUT

out = struct;
out.location = [xk amatrix];
out.time = sum(t);
out.parllel_time = sum(tp);
out.fun_val = [F(x0); fun_v];
out.norm_bias = [bias(x0); bias_v];

disp('AMFD is Done!')
end



%% MakeClustrs function
function cluster_stc = MakeClusters(net_noise)

K = net_noise.K;
m = net_noise.anchors;
N = K - m;
Atilde = net_noise.Matrices.A_tilde;
Qtilde = net_noise.Matrices.Q_tilde;
Btilde = net_noise.Matrices.B_tilde;

cluster_stc = struct;
clusters = cell(N, 1);
perm_nodes = randperm(N);
for i=1:N-1
    clusters{i}.nodes = setdiff(perm_nodes((i - 1) + 1:i), N + 1:K);
end
clusters{N}.nodes = setdiff(perm_nodes((N - 1) + 1:end), N + 1:K);

for i = 1:N % setting data for each cluster
    nodes = clusters{i}.nodes; % relative indices of nodes in nonanchor_indices
    edges_with_nodes = find(sum(Atilde(:, nodes), 2)==1);
    clusters{i}.anchor_edges = edges_with_nodes; % finding edges with anchors not in cluster
    edges_with_nodes = find(sum(Qtilde(:, nodes)~=0, 2)); % finding edges with nonanchors not in cluster
    clusters{i}.nonanchor_edges = edges_with_nodes;
    if length(clusters{i}.anchor_edges) > 1
        connected_anchors = find(sum(Btilde(clusters{i}.anchor_edges, :)) > 0);
    else
        connected_anchors = find(Btilde(clusters{i}.anchor_edges, :) > 0);
    end
    clusters{i}.anchor_neigbhors = connected_anchors;
    relavent_nonanchors = setdiff(1:N, nodes);
    connected_nonanchors = (sum(abs(Qtilde(clusters{i}.nonanchor_edges, relavent_nonanchors))) > 0);
    clusters{i}.nonanchor_neigbhors = relavent_nonanchors(connected_nonanchors);
end
cluster_stc.clusters = clusters;
end