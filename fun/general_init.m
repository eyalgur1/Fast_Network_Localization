function [Qtilde,QDtilde,Atilde,ADtilde,Btilde,BDtilde,Ptilde,a,Xreal,xreal,...
    IMnet,IMDnet,Laplaciannet,K,m,N,n,Ne,Neq,invPtilde,...
    amatrix,F,bias,x_const,x_u_coef,v_x_coef,v_const]=general_init(Network)
%% INFORMATION

% DESCRIPTION: returns matrices and functions which are relavent for all
%              algorithms.

% INPUTS:
% * Network - the true network structure

%% OUTPUT

% Independent of relaizations
Xreal=Network.Matrices.X_real; % true location matrix
xreal=Network.Matrices.x_real; % true location vector
K=Network.K; % total number of sensors
m=Network.GI.num_anchors; % total number of anchors
N=K-m; % total number of non-anchors
n=size(Network.Matrices.X_real,1); % dimension
a=reshape(Xreal(:,N+1:end),n*m,1); % anchors real location
delta = Network.GI.delta; % radius of uncertainty balls
Ne=size(Network.Matrices.Q_tilde,1)+size(Network.Matrices.A_tilde,1); % total number of non-anchor edges
Neq=size(Network.Matrices.Q_tilde,1); % total number of non-anchor to non-anchor edges

% Matrices independent of relaizations for algorithms with anchor certainty
Ptilde=Network.Matrices.P_tilde;
invPtilde=inv(Ptilde);
Atilde=Network.Matrices.A_tilde;
Btilde=Network.Matrices.B_tilde;
Qtilde=Network.Matrices.Q_tilde;
IMnet=Network.Matrices.IM_net;
Laplaciannet=Network.Matrices.Laplacian_net;

% Matrices and data independent of relaizations for algorithms with anchor uncertainty
%S = net.Matrices.S;
%P = net.Matrices.P;  % P = S'*S+diag([zeros(n*N,1);ones(n*m,1)])  % the matrix P required for SNAM
%Lip_cen=net.GI.Lip_cen;  % = 2*norm(0.5*P);
%Lip_dist=net.GI.Lip_dist;  % = 2*max(diag(P))+1; %Lip_very_dist=4*K+1;
%edges = net.GI.edges;
[~,edge_1] = find(IMnet>0);
[edge_2,~] = find(IMnet'<0);


% Matrices dependent of relaizations for all algorithms
amatrix=Network.Matrices.noised_anchors; % anchors noised (if any) location matrix
IMDnet=Network.Matrices.IMD_net;


% Matrices dependent of relaizations for algorithms with anchor certainty
QDtilde=Network.Matrices.QD_tilde;
ADtilde=Network.Matrices.AD_tilde;
BDtilde=Network.Matrices.BD_tilde;
x_coef=(amatrix*Btilde'*Atilde);
v_x_coef=[QDtilde;ADtilde]'; % [QD;AD] for the objective function
x_const=x_coef*invPtilde; % (invP*A'*B)*a' for the objective function
x_u_coef=v_x_coef'*invPtilde; % invP*[QD', AD'] for the objective function
v_const=[zeros(2,size(QDtilde,1)),-amatrix*BDtilde']; % [zeros(size(QD,1),1);-BD*a] for the objective function


% Function value (which depends on the relaization) for anchor ball uncertainty,
% which is applicable for all compared methods
dist=0.5*sum(abs(IMDnet),2)';
F=@(x) 0.5*sum((sqrt(sum((x(:,edge_1) - x(:,edge_2)).^2)) - dist).^2) + 0.5*sum(sum((x(:,N+1:end) - amatrix).^2))...
    + sum(((sqrt(sum((x(:,N+1:end) - amatrix).^2)) - delta) > 10^(-10))).*10^(10);
%F=@(x,u)(sum(sum(((x*Ptilde').*x)))-sum(sum(2*u.*(x*v_x_coef+v_const)))-sum(sum(2*x_coef.*x))); % objective function
%F=@(x) sum((sqrt(sum((x(:,edge_1) - x(:,edge_2)).^2)) - dist).^2);

% Bias (which is independent of the realzation and certainty/uncertianty)
bias=@(x)norm(x - Xreal);
%RMSE=@(x)(sqrt(1/K))*sqrt(sum(sum(([x,amatrix]-Xreal).^2))); % RMSE function

end