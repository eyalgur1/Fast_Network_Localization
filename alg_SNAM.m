function out = alg_SNAM(net_noise, net_original, x0, r, s, max_iter, cen)
%% INITIALIZATION

IMDnet = net_noise.Matrices.IMD_net;
K = net_noise.K;
m = net_noise.anchors;
N = K - m;
n = size(net_noise.Matrices.X_real, 1);
amatrix = net_noise.Matrices.noised_anchors;
F = net_noise.GI.F;
bias = net_noise.GI.bias;
ball_radius = net_noise.GI.ball_radius;  % radius of ball uncertainty (0 or greater)
S = net_original.net.Matrices.S;
P = net_original.net.Matrices.P;
edges = net_original.net.GI.edges;
edge_1 = edges(:,1); edge_2 = edges(:,2);
Ne = length(edge_1);
lin_const = [zeros(n*N,1);reshape(amatrix,n*m,1)];  % A*a
lin_coeff = kron(diag(0.5*sum(abs(IMDnet),2)), eye(2))*S;  % D*S
II = eye(n); e1 = II(:,1);
if cen
    L = net_original.net.GI.Lip_cen;  % centralized Lip constant
else 
    L = net_original.net.GI.Lip_dist;  % distributed Lip constant (max(diag(P)) is the maximal node degree)
    % L = 4*K+1; % very distributed 
end

bias_v = zeros(max_iter, 1); fun_v = bias_v; t = bias_v; tp = bias_v;
x0 = reshape(x0,n,K);
x0(:, N+1:end) = amatrix;  % set amatrix as the initial point for the anchors
xj = x0;
iter = 0;
out_iter=0;

%% ALGORITHM
fprintf('SNAM | Iter=%2d | Bias=%10.10f | F=%10.10f\n', 0, bias(x0), F(x0))  

while iter < max_iter  

    starta = tic;
    out_iter = out_iter + 1;

    % update u vectors efficiently
    u = xj(:, edge_1) - xj(:, edge_2);  % all current u vectors
    norm_x = sqrt(sum(u.^2));
    norm_x_pos = find(norm_x);  % only u vectors with positive norm (differentiable)
    u = u(:,norm_x_pos)./sqrt(sum(u(:,norm_x_pos).^2));  % update differentiable u vectors
    rep = setdiff(1:Ne,norm_x_pos);  % indices of non-differentiable u vectors
    if ~isempty(rep)  % find if there are non-differentiable u vectors, and update them
        u(:,rep) = e1;
    end

    nk = s + 2^(floor(out_iter/r)) - 1;  % set number of inner iterations
    
    yj = reshape(xj,K*n,1); tstep=1;  % set aux variable for FISTA and step-size
    for j = 1:nk
        startb=tic;
        iter = iter + 1;

        tstep_prev=tstep;
        xj_prev=xj;

        cj = yj - (1/L)*(P*yj - (lin_coeff'*reshape(u,n*Ne,1) + lin_const));
        xj = reshape(cj,n,K);  % update non-anchor
        diff_a = xj(:,N+1:end) - amatrix;  % update anchors by projection onto balls
        xj(:,N+1:end) = amatrix + ball_radius*(diff_a./max(sqrt(sum(diff_a.^2)),ball_radius));

        tstep = (1+sqrt(1+4*tstep^2))/2;  % update step-size

        yj = reshape(xj + ((tstep_prev-1)/tstep)*(xj-xj_prev),n*K,1);  % update aux variable

        if j == 1
            t(iter) = toc(starta);  % add the u time calculation only once
        else
            t(iter) = toc(startb);
        end
        tp(iter) = (1/K)*t(iter);  % parallel time estimation

        fun_v(iter) = F(xj);
        bias_v(iter) = bias(xj);

        if ~mod(iter, 1000) || iter == 1
            fprintf('SNAM, cen=%1d, s=%2d | Iter=%2d | Bias=%10.10f | F=%10.10f\n', cen, s, iter, bias_v(iter), fun_v(iter))  
        end

        if iter >= max_iter
            break
        end
    end
    
end

%% OUTPUT
out = struct;
out.location = xj;
out.time = sum(t);
out.parllel_time = sum(tp);
out.fun_val = [F(x0); fun_v];
out.norm_bias = [bias(x0); bias_v];
out.s = s;
out.r = r;
disp('SNAM is Done!')