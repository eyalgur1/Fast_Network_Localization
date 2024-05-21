%% %%% Calculate Averages Script %%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc

%% Hyper-Parameters

net_type = 1;
num_real = 25;
s_SNAM = [40];
sigma_d_intervals = 6;
max_iter = 10000;


%% Initialization

if net_type == 1
    methods = {'SNAM_cen', 'SNAM_dist', 'DC', 'SGO', 'SOCP'};
    %sigma_d_lb = 0.02; sigma_d_ub = 0.5;
    sigma_d = logspace(-3, -1.5, 4);
    num_nets = 1; avg_fact = 1/num_real;
    net_str = 'net1000bench';
    out_folder = 'output\2024_01_04_16_25_34\net1000bench\';
    data_folder = 'datasets\net1000bench.mat';
elseif net_type == 2
    methods = {'SNAM_cen', 'SNAM_dist', 'AMFC', 'SGO'};
    sigma_d = logspace(-5, -1.5, sigma_d_intervals);
    num_nets = 30; avg_fact = 1/(num_real + num_nets);
    out_folder = 'output\Rand1000_complete_data\';
    data_folder = 'datasets\net';
elseif net_type == 3
    methods = {'SNAM_cen', 'SNAM_dist', 'DC', 'SGO'};
    sigma_d = logspace(-4, -4, 1);
    num_nets = 1; avg_fact = 1/num_real;
    net_str = 'net10000';
    out_folder = 'output\2024_03_31_16_56_34\net10000rand\';
    data_folder = 'datasets\net10000.mat';
end

if find(contains(methods, 'SNAM_cen')) ~= 0
    methods = method_SNAM(methods, 'SNAM_cen', s_SNAM);
end
if find(contains(methods, 'SNAM_dist')) ~= 0
    methods = method_SNAM(methods, 'SNAM_dist', s_SNAM);
end

avgs = struct;
RMSE_vals = zeros(size(methods, 2) + length(s_SNAM), length(sigma_d));
cov_mat_fun = @(X_out, X_avg)(X_out - X_avg)*(X_out - X_avg)';


% Plot initilization
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultLegendInterpreter' ,'latex');
set(groot, 'defaultLineLineWidth',1)
set(groot, 'defaultAxesFontSize',14)
colors = distinguishable_colors(10);

names_to_legend = {};
for mm = 1:length(methods)
    met = methods{mm};
    if contains(met, 'SNAM_cen')
        names_to_legend{end + 1} = ['$\textrm{SNAM,\ centralized,\ }s=', erase(met, 'SNAM_cen_s'), '$'];
    elseif contains(met, 'SNAM_dist')
        names_to_legend{end + 1} = ['$\textrm{SNAM,\ distributed,\ }s=', erase(met, 'SNAM_dist_s'), '$'];
    else
        names_to_legend{end + 1} = ['$\textrm{', met,'}$'];
    end
end


%% Calculation of Averages

% For each method and sigma, calculate the folloiwng:
%   1. Average function value (over realizations and networks)
%   2. Average norm of bias (over realizations and networks)
%   3. Average total run time (over realizations and networks)
%   4. Average parallel run time (over realizations and networks, if applicable)
%   5. Average RMSE (over realizations and then over networks)

for mm = methods
    method = mm{1};
    load([out_folder, method, '\', method, '.mat'])
    for sigma = sigma_d
        sigma_str = ['sigma_d_', strrep(strrep(num2str(sigma), '.', ''), '-', '')];
        fprintf(['*****', method, ' | sigma=%10.10f*****\n'], sigma)
        avgs.(method).(sigma_str).funV = 0; avgs.(method).(sigma_str).bias = 0; avgs.(method).(sigma_str).time = 0; avgs.(method).(sigma_str).ptime = 0; avgs.(method).(sigma_str).RMSE = 0;
        for nn = 1:num_nets
            if net_type == 2
                net_str = ['net', num2str(nn)];
            end
            avgs.(method).(sigma_str).X.(net_str) = 0;
            for rr = 1:num_real
                avgs.(method).(sigma_str).funV = avgs.(method).(sigma_str).funV + avg_fact*O.(sigma_str).(net_str).realization{rr}.fun_val;
                avgs.(method).(sigma_str).bias = avgs.(method).(sigma_str).bias + avg_fact*O.(sigma_str).(net_str).realization{rr}.norm_bias;
                avgs.(method).(sigma_str).time = avgs.(method).(sigma_str).time + avg_fact*O.(sigma_str).(net_str).realization{rr}.time;
                if isfield(O.(sigma_str).(net_str).realization{1}, 'parllel_time')
                    avgs.(method).(sigma_str).ptime = avgs.(method).(sigma_str).ptime + avg_fact*O.(sigma_str).(net_str).realization{rr}.parllel_time;
                end
                avgs.(method).(sigma_str).X.(net_str) = avgs.(method).(sigma_str).X.(net_str) + (1/num_real)*O.(sigma_str).(net_str).realization{rr}.location;
            end
            avgs.(method).(sigma_str).cov.(net_str) = 0;
            X_avg = avgs.(method).(sigma_str).X.(net_str);
            for rr = 1:num_real
                X_out = O.(sigma_str).(net_str).realization{rr}.location;
                avgs.(method).(sigma_str).cov.(net_str) = avgs.(method).(sigma_str).cov.(net_str) + (1/num_real)*cov_mat_fun(X_out, X_avg);
            end
            avgs.(method).(sigma_str).RMSE = avgs.(method).(sigma_str).RMSE + (1/num_nets)*sqrt(trace(avgs.(method).(sigma_str).cov.(net_str)));
        end
    end
end


% For each sigma, calculate average RMSE (over networks)
if net_type == 1
    load(data_folder); ball_radius = 0.005;
elseif net_type == 3
    load(data_folder); ball_radius = 0.001;
end
for sigma = sigma_d
    well_cond_mat = 0;
    sigma_str = ['sigma_d_', strrep(strrep(num2str(sigma), '.', ''), '-', '')];
    fprintf('*****CRLB | sigma=%10.10f*****\n', sigma)
    avgs.CRLB.(sigma_str) = 0;

    for nn = 1:num_nets
        if net_type == 2
            load([data_folder, num2str(nn), '.mat']); ball_radius = 0.005;
        end
        [r_cond, CRLB_net] = CRLB_SNL(net, sigma, ball_radius/3);
        if r_cond > 10^(-15)
            well_cond_mat = well_cond_mat + 1;
            avgs.CRLB.(sigma_str) = avgs.CRLB.(sigma_str) + CRLB_net;
        end
        avgs.CRLB.(sigma_str) = (1/well_cond_mat)*avgs.CRLB.(sigma_str);
    end
end


%% Plot Function Values and Bias for smallest sigma_d

figure(1); hold on
figure(2); hold on
sigma_str = ['sigma_d_', strrep(strrep(num2str(sigma_d(1)), '.', ''), '-', '')];
sigma_d_title = ['\textrm{', num2str(sigma_d(end)), '}'];

for mm = 1:length(methods)
    method = methods{mm};

    switch method
        case 'SNAM_cen_s40'
            figure(1); semilogy(avgs.(method).(sigma_str).funV, '--', 'Color', colors(1, :));
            figure(2); semilogy(avgs.(method).(sigma_str).bias, '--', 'Color', colors(1, :));

        case 'SNAM_cen_s80'
            figure(1); semilogy(avgs.(method).(sigma_str).funV, 'Color', colors(1, :));
            figure(2); semilogy(avgs.(method).(sigma_str).bias, 'Color', colors(1, :));

        case 'SNAM_dist_s40'
            figure(1); semilogy(avgs.(method).(sigma_str).funV, '--', 'Color', colors(2, :));
            figure(2); semilogy(avgs.(method).(sigma_str).bias, '--', 'Color', colors(2, :));

        case 'SNAM_dist_s80'
            figure(1); semilogy(avgs.(method).(sigma_str).funV, 'Color', colors(2, :));
            figure(2); semilogy(avgs.(method).(sigma_str).bias, 'Color', colors(2, :));

        case 'AMFC'
            figure(1); semilogy(avgs.(method).(sigma_str).funV, 'Color', colors(3, :));
            figure(2); semilogy(avgs.(method).(sigma_str).bias, 'Color', colors(3, :));

        case 'DC'
            figure(1); semilogy(avgs.(method).(sigma_str).funV, 'Color', colors(5, :));
            figure(2); semilogy(avgs.(method).(sigma_str).bias, 'Color', colors(5, :));

        case 'SGO'
            figure(1); semilogy(avgs.(method).(sigma_str).funV, 'Color', colors(6, :));
            figure(2); semilogy(avgs.(method).(sigma_str).bias, 'Color', colors(6, :));

        case 'SF'
            figure(1); semilogy(avgs.(method).(sigma_str).funV, 'Color', colors(7, :));
            figure(2); semilogy(avgs.(method).(sigma_str).bias, 'Color', colors(7, :));

        case 'SOCP'
            figure(1); semilogy(avgs.(method).(sigma_str).funV*ones(1, max_iter + 1), 'Color', colors(9, :));
            figure(2); semilogy(avgs.(method).(sigma_str).bias*ones(1, max_iter + 1), 'Color', colors(9, :));
    end
end

figure(1); set(gca, 'YScale', 'log')
xlim([0 max_iter + 1]); xlabel('$\textrm{Iterations}$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\textrm{ML\ Function\ Value}$', 'Interpreter', 'latex', 'FontSize', 14)
title(['$\textrm{Average\ ML\ Function\ Value\ for\ }\sigma_{d}=', sigma_d_title, '$'], 'Interpreter', 'latex', 'FontSize', 14)
legend(names_to_legend); hold off

figure(2); set(gca, 'YScale', 'log')
xlim([0 max_iter + 1]); xlabel('$\textrm{Iterations}$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\textrm{ML\ Function\ Value}$', 'Interpreter', 'latex', 'FontSize', 14)
title(['$\textrm{Average\ Squared\ Norm\ of\ Bias\ for\ }\sigma_{d}=', sigma_d_title, '$'], 'Interpreter', 'latex', 'FontSize', 14)
legend(names_to_legend); hold off


%% Plot RMSE and CRLB

figure(3); hold on
for mm = 1:length(methods)
    method = methods{mm};
    RMSE_method = zeros(size(sigma_d));

    for sig = 1:length(sigma_d)
        sigma_str = ['sigma_d_', strrep(strrep(num2str(sigma_d(sig)), '.', ''), '-', '')];
        RMSE_method(sig) = avgs.(method).(sigma_str).RMSE;
    end
    plot(RMSE_method);
end

CRLB_vals = zeros(size(sigma_d));
for sig = 1:length(sigma_d)
    sigma_str = ['sigma_d_', strrep(strrep(num2str(sigma_d(sig)), '.', ''), '-', '')];
    CRLB_vals(sig) = avgs.CRLB.(sigma_str);
end
plot(CRLB_vals, 'Color', 'black');
legend([names_to_legend, 'CRLB']); hold off



%% Auxiliary Functions

% Create str of methods
function methods = method_SNAM(methods, str, s_SNAM)
ind = find(contains(methods, str));
methods = setdiff(methods, methods{ind});
new_meth = cell(1, length(s_SNAM));
for ss = 1:length(s_SNAM)
    s_val = s_SNAM(ss);
    new_meth{ss} = [str, '_s', num2str(s_val)];
end
methods = [methods new_meth];
end

% Calculate CRLB
function [r_cond, CRLB] = CRLB_SNL(net, sigma_dist, sigma_a)
X = net.Matrices.X_real;
[n, K] = size(X);
N = K - net.anchors;
FIM = zeros(n*K);

for i = 1:K
    for j = net.node{1, i}.neighbors
        ij = X(:, i) - X(:, j);
        normij = norm(ij);
        mij = (sigma_dist^(-2))*(ij*ij')./normij;
        FIM(n*(i-1) + 1:n*i, n*(j-1) + 1:n*j) = -mij;  % off-diagonal
        FIM(n*(j-1) + 1:n*j, n*(i-1) + 1:n*i) = -mij;  % off-diagonal

        FIM(n*(i-1) + 1:n*i, n*(i-1) + 1:n*i) = FIM(n*(i-1) + 1:n*i, n*(i-1) + 1:n*i) + mij;  % diagonal
    end
    if i > N  % if anchor
        FIM(n*(i-1) + 1:n*i, n*(i-1) + 1:n*i) = FIM(n*(i-1) + 1:n*i, n*(i-1) + 1:n*i) + (sigma_a^(-2))*eye(n);
    end
end

r_cond = rcond(FIM);
CRLB = sqrt(trace(inv(FIM)));
end
