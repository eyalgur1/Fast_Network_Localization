function out = alg_AMFC(net_original, net_noise, x0, max_iter)
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
invPtilde = net_original.net.Matrices.invP_tilde;
Ne = size(Qtilde, 1) + size(Atilde, 1);
amatrix = net_noise.Matrices.noised_anchors;
F = net_noise.GI.F;
bias = net_noise.GI.bias;
x_coef = amatrix*(Btilde'*Atilde);
v_x_coef = [QDtilde; ADtilde]';  % [QD;AD] for the objective function
x_const = x_coef*invPtilde;  % (invP*A'*B)*a' for the objective function
x_u_coef = v_x_coef'*invPtilde;  % invP*[QD', AD'] for the objective function
v_const = [zeros(2, size(QDtilde,1)), -amatrix*BDtilde'];  % [zeros(size(QD,1),1);-BD*a] for the objective function

bias_v = zeros(max_iter, 1); fun_v = bias_v; t = bias_v; tp = bias_v;
x0 = reshape(x0, n, K);
x0(:, N+1:end) = amatrix;  % set amatrix as the initial point for the anchors
uk = zeros(n, Ne);
iter = 0;

%% ALGORITHM
fprintf('AMFC | Iter=%2d | Bias=%10.10f | F=%10.10f\n', 0, bias(x0), F(x0))  

while iter < max_iter
    
    starta=tic;
    iter = iter + 1;
    
    starti = tic;
    xk = x_const + uk*x_u_coef; % invP*([QD', AD']*u + A'*B*a) - update the x block   
    endi = toc(starti);
 
    [uk, ~, endb] = u_update(xk, uk, v_x_coef, v_const, Ne, n, K, []); % update the u block
    
    t(iter) = toc(starta);  % upper bound on the non-parallelizable running time
    tp(iter) = endi + endb;  % running time estimation
    bias_v(iter) = bias([xk amatrix]);
    fun_v(iter) = F([xk amatrix]);
    %anc_AMFC = amatrix(:, 1)

    if ~mod(iter, 1000) || iter == 1
        fprintf('AMFC | Iter=%2d | Bias=%10.10f | F=%10.10f\n', iter, bias_v(iter), F([xk amatrix]))  
    end
end

%% OUTPUT

out = struct;
out.location = [xk amatrix];
out.time = sum(t);
out.parllel_time = sum(tp);
out.fun_val = [F(x0); fun_v];
out.norm_bias = [bias(x0); bias_v];
disp('AMFC is Done!')