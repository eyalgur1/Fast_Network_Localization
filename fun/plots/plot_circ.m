K= 1000;
m = 20;
delta = 0.05/10;
sigma = delta/3;

%% Check for good radius for ball uncertainty

N = K - m;
figure(1)
h = circle(-0.5 + rand(m, 1), -0.5 + rand(m, 1), delta,m);
scatter(-0.5 + rand(N, 1), -0.5 + rand(N,1 ), '.','blue')
hold off


%% Check for good noise for ball uncertainty

rep = 1:K;
X = zeros(2,K);

while ~isempty(rep)
    X(:,rep) = sigma*randn(2,length(rep));
    rep = find(sqrt(sum(X.^2)) > delta);
end

figure(2)
scatter(X(1,:),X(2,:),'.');
hold on
axis equal
x_ax = -delta:delta/100:delta;
y_ax =sqrt(delta^2 - x_ax.^2);
plot(x_ax,y_ax, 'red')
plot(x_ax,-y_ax, 'red')
hold off


%%
% Auxiliary function
function h = circle(x, y, r, m)
hold on
th = 0:pi/50:2*pi;
for i = 1:m
    xunit = r * cos(th) + x(i);
    yunit = r * sin(th) + y(i);
    scatter(x(i), y(i), 'red','.');
    h = plot(xunit, yunit, 'red');
end
axis equal
end