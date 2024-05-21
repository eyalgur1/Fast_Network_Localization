function [net_noise, x0]=create_realization(net_original,sigma_d,sigma_a,ball_radius,x0_bound)
%% INFORMATION
% DESCRIPTION: creates a truncated benchmark network structure based on the networks
%              in https://web.stanford.edu/~yyye/Col.html.

% INPUTS:
% * sigma_d - if greater than 0 then noises the distances according to sigma_d,
%           otherwise adds no noise.
% * sigma_a - if greater than 0 then noises the anchors according to sigma_a,
%           otherwise adds no noise.

%% CREATING THE NETWORK

m = net_original.net.anchors;
PP = net_original.net.Matrices.X_real;
dd = net_original.net.Matrices.true_distances;

n=size(PP,1);
K=size(PP,2);
N = K - m;
SS=PP'; % Location matrix, each row is a sensor

distances=zeros(K);
for i=1:n % calculated the squared real distances beween all sensors
    distances=distances+(SS(:,i)*ones(1,K)-ones(K,1)*SS(:,i)').^2;
end
distances=sqrt(distances).*(dd>0); % returns the true distances between all pairs of neighbors

truncated_given_distances=dd; % dd is a given distance measurement between neighbors (noised or not)

% If choosing to noise, then new distance measurements are created between
% all pairs. However, only the distances between the given set of neighbors
% (as defined by the input matrix dd) are used in the code. If choosing not
% to noise, then the given noises in the input matrix dd are used.
if sigma_d>0
    distance_noise=sigma_d*full(sprandsym(dd)); % creates a symmetric matrix with the structure of dd with non-zero elemnts from Gauss(0,sigma)
    noised_distances=abs(distances+distance_noise); % sesnor adjacency (symmetric) matrix of with noised distances, with no negative entries
else
    noised_distances=abs(truncated_given_distances); % sesnor adjacency (symmetric) matrix of with noised distances, with no negative entries (in this case dd is already noised)
end

true_distances=distances.*(dd>0);

if sigma_a>0
    noises = sigma_a*randn(2,m);
    noise_a = noises;
    ind_rep = find(sqrt(sum(noises.^2)) > ball_radius);
    noise_a(:,ind_rep) = ball_radius*noises(:,ind_rep)./sqrt(sum(noises(:,ind_rep).^2));
    noised_anchors = PP(:, end-m+1: end) + noise_a;
end

% if ~isempty(varargin)
%     R=varargin{1}; % radius as given by the user
% else
%     R=max(max(noised_distances)); % calculate real radius
% end

num_edges=nnz(truncated_given_distances)/2; % number of edges

net_noise=struct;
net_noise.K=K;
net_noise.anchors=m;
if n==2
    net_noise.positions=(SS(:,1)+SS(:,2)*1i)';
else
    net_noise.positions=SS';
end
net_noise.is_anchor=logical([zeros(1,K-m),ones(1,m)]);

IM=zeros(num_edges,K);
IMD=IM;
Nodes = cell(1,K);
cur_l=0;
for i=1:K
    %create edges for node i
    i_neighbors=find(truncated_given_distances(i,i+1:K))+i;
    n_i_neighbors=length(i_neighbors);
    %update the incidence matrix and the incidence matrix with distances
    IM(cur_l+1:cur_l+n_i_neighbors,i)=1;
    IM(cur_l+1:cur_l+n_i_neighbors,i_neighbors)=-eye(n_i_neighbors);
    IMD(cur_l+1:cur_l+n_i_neighbors,i)=noised_distances(i,i_neighbors)';
    IMD(cur_l+1:cur_l+n_i_neighbors,i_neighbors)=-diag(noised_distances(i,i_neighbors)');
    cur_l=cur_l+n_i_neighbors;
    
    %create the node
    Nodes{i}.neighbors=[find(truncated_given_distances(i,1:i-1)),i_neighbors];
    Nodes{i}.N_neigh=length(Nodes{i}.neighbors);
    Nodes{i}.distance=abs(noised_distances(i,Nodes{i}.neighbors));
    Nodes{i}.real_distance=abs(distances(i,Nodes{i}.neighbors));
    Nodes{i}.position=net_noise.positions(i);
    if i>K-m
        Nodes{i}.is_anchor=true;
    else
        Nodes{i}.is_anchor=false;
    end
end
for j=1:K
    for l=1:Nodes{j}.N_neigh
        %updating the position of current node in neighbor's neighbor list
        neighbor_index=Nodes{j}.neighbors(l);
        Nodes{j}.neigh_pos(l)=find(Nodes{neighbor_index}.neighbors==j);
    end
end

%% OUTPUT

edges_with_anchor=(sum(abs(IM(:,K-m+1:K)),2)==1);
A_tilde=abs(IM(edges_with_anchor,1:K-m));
B_tilde=abs(IM(edges_with_anchor,K-m+1:K));
AD_tilde=abs(IMD(edges_with_anchor,1:K-m));
BD_tilde=abs(IMD(edges_with_anchor,K-m+1:K));
a=reshape(SS(K-m+1:K,:)',m*n,1); %locations of anchors
X_real=SS';
x_real=reshape(X_real,K*n,1);
edges_of_nona_nona=(sum(abs(IM(:,1:K-m)),2)==2);
Q_tilde=IM(edges_of_nona_nona,1:K-m);
QD_tilde=IMD(edges_of_nona_nona,1:K-m);
P_tilde=Q_tilde'*Q_tilde+A_tilde'*A_tilde;
Laplacian=IM'*IM;

GI = struct;
Matrices = struct;

net_noise.node = Nodes;
GI.num_anchors = m;
GI.num_non_anchors = K-m;
GI.num_rel_edges = sum(edges_with_anchor)+sum(edges_of_nona_nona);
GI.anchor_indices = K-m+1:K;
GI.ball_radius = ball_radius;
GI.avg_neigh = mean(diag(Laplacian));

% If choosing not to moise, then noise_D is the std of the noise.
% Otherwise, the std is given by the input matrix dd after truncation. In
% this case the std must be calculated.
if sigma_d>0
    GI.Noise_d = sigma_d;
else
    true_truncated_distances=logical(truncated_given_distances).*distances;
    GI.Noise=std(unique(true_truncated_distances-truncated_given_distances));
end

if sigma_a>0
    GI.Noise_a = sigma_a;
end

Matrices.Q_tilde = sparse(Q_tilde);
Matrices.QD_tilde = sparse(QD_tilde);
Matrices.A_tilde = sparse(A_tilde); %indicator matrix - rows are edges between non anchors and anchors, columns are non-anchors nodes
Matrices.AD_tilde = sparse(AD_tilde);
Matrices.B_tilde = sparse(B_tilde); %indicator matrix - rows are edges between non anchors and anchors, columns are anchors nodes
Matrices.BD_tilde = sparse(BD_tilde);
Matrices.P_tilde = P_tilde;
Matrices.a_real = a;
Matrices.X_real = X_real;
Matrices.x_real = x_real;
Matrices.IM_net = sparse(IM);
Matrices.IMD_net = sparse(IMD);
Matrices.Laplacian_net = Laplacian;
Matrices.true_distances=sparse(true_distances);
if sigma_d>0
    Matrices.dd_noise=sparse(distance_noise);
end

net_noise.GI = GI;
net_noise.Matrices = Matrices;
net_noise.Matrices.noised_distances=sparse(noised_distances);
net_noise.Matrices.noised_anchors = noised_anchors;
%anc_noise = noised_anchors(:, 1)


[~,edge_1] = find(Matrices.IM_net>0);
[edge_2,~] = find(Matrices.IM_net'<0);
dist=0.5*sum(abs(Matrices.IMD_net),2)';
%F = @(x) 0.5*sum((sqrt(sum((x(:,edge_1) - x(:,edge_2)).^2)) - dist).^2) + 0.5*sum(sum((x(:,N+1:end) - noised_anchors).^2))...
%    + sum(((sqrt(sum((x(:,N+1:end) - noised_anchors).^2)) - GI.ball_radius) > 10^(-10))).*10^(10);
F = @(x) 0.5*sum((sqrt(sum((x(:,edge_1) - x(:,edge_2)).^2)) - dist).^2) + 0.5*sum(sum((x(:,N+1:end) - noised_anchors).^2));
net_noise.GI.F = F;
net_noise.GI.bias = @(x)norm(x - X_real,'fro');

x0 = -x0_bound + 2*x0_bound*rand(n*K, 1);

end