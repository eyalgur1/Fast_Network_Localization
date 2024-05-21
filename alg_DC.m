function out = alg_DC(net_original, net_noise, x0, max_iter, regul)
%% INITIALIZATION

% Hyper-parameters
K = net_noise.K;
m = net_noise.anchors;
N = K - m;
n = size(net_noise.Matrices.X_real, 1);
mat = [ones(n,1) -ones(n,1)];


% Independent of realization
sum_Bl = sparse(diag([zeros(1,N) 2*ones(1,m)]));
H = 2*net_original.net.Matrices.Laplacian_net + sum_Bl + sparse(2*regul*eye(K));  % should be set differently if the noises are known to the user


% Depends on realization
noised_distances = net_noise.Matrices.noised_distances;
[row,col] = find(tril(net_original.net.Matrices.true_distances));  % find indices of neighboring sensors i~j witout j~i
noised_distances_vec = zeros(length(row),1);
for i =1:length(row)
    noised_distances_vec(i) = noised_distances(col(i), row(i));
end
amatrix = net_noise.Matrices.noised_anchors;
barX = sparse([zeros(n,N) amatrix]);  % each anchor l column should be multiplied by inv(covariance_l) if we posses this knowledge
update_X = @(Y)(Y + barX*sum_Bl)/H;  % the update of the location matrix X
F = net_noise.GI.F;  % objective function
bias = net_noise.GI.bias;  % bias function


iter = 0;
X = reshape(x0, n, K);  % starting point
X(:, N+1:end) = amatrix;
bias_v = zeros(max_iter, 1); fun_v = bias_v; t = bias_v;
fprintf('DC | Iter=%2d | Bias=%10.10f | F=%10.10f\n', 0, bias(X), F(X))  

%% Run Algorithm

starta = tic;
while iter < max_iter
    iter = iter + 1;

    Y = update_Y(X, noised_distances_vec, row, col, regul, n, K, mat);
    X = update_X(Y);

    t(iter) = toc(starta);
    fun_v(iter) = F(X);
    bias_v(iter) = bias(X);

    if ~mod(iter, 1000) || iter == 1
            fprintf('DC | Iter=%2d | Bias=%10.10f | F=%10.10f\n', iter, bias_v(iter), fun_v(iter))  
    end
end


%% Output

out = struct;
out.location = X;
out.time = sum(t);
out.fun_val = [F(reshape(x0, n, K)); fun_v];
out.norm_bias = [bias(reshape(x0, n, K)); bias_v];
disp('DC is Done!')
end


%% Auxiliary Function (updating matrix Y)
function Y = update_Y(X, noised_distances_vec, row, col, r, n, K, mat)
Y = zeros(n,K);
v = X(:,col) - X(:,row);
norm_v = sqrt(sum(v.^2));
if sum(norm_v==0)>0
    norm_v(norm_v==0) = 1;
end
v = v./norm_v;
for i = 1:length(row)
    inde  = [col(i),row(i)];
    Y(:,inde) = Y(:,inde) + noised_distances_vec(i)*mat.*v(:,i);
end
Y = 2*(Y+r*X);
end