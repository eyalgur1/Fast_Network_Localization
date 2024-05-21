clear; format compact
rng(1);
formattedDateTime = "2024_03_31_16_56_34";  % string(char(datetime('now', 'Format', 'yyyy_MM_dd_HH_mm_ss')));
%maxNumCompThreads(1);

%% Hyper-parameters

net_type = 1;  % 1 - Benchmark 1000, 2 - Random 1000, 3 - Random 10000
num_nets = 1;  % only for net_type = 2
methods = {'DC'};
max_iter = 10000;
real_arr = 1:1;
sigma_d_intervals = 1;
sigma_d_lb = 0.02;  % sets sigma_d to lb% to ub% of shortest edge
sigma_d_ub = 0.5;
x0_bound = 0.01;
regul = 0.0001;

% Hyper-parameters for SNAM
s_SNAM = [40]; r_SNAM = 1000;


% set noises and ball uncertainty radius
if net_type == 1
    sigma_d = logspace(-3, -1.5, 4);
    %sigma_d = sqrt(logspace(log10(0.001*sigma_d_lb), log10(0.001*sigma_d_ub), sigma_d_intervals));
    ball_radius = 0.005; sigma_a = ball_radius/3;
elseif net_type == 2
    sigma_d = 0.01;%logspace(-3, -1, 5);
    %sigma_d = logspace(log10(0.000431*sigma_d_lb), log10(0.000431*sigma_d_ub), sigma_d_intervals);
    %sigma_d = logspace(-5, -1.5, sigma_d_intervals);
    ball_radius = 0.005; sigma_a = ball_radius/3;
elseif net_type == 3
    %sigma_d = logspace(log10(0.000431*sigma_d_lb), log10(0.000431*sigma_d_ub), sigma_d_intervals);
    sigma_d = logspace(-4, -4, sigma_d_intervals);
    ball_radius = 0.001; sigma_a = ball_radius/3;
end


real_struct = struct;
for nn = 1:num_nets
    if net_type == 1
        net_original = load('datasets/net1000bench.mat');
        net_str = 'net1000bench';
    elseif net_type == 2
        net_original = load(['datasets/net', num2str(nn), '.mat']);
        net_str = ['net', num2str(nn)];
    elseif net_type == 3
        net_original = load('datasets/net10000.mat');
        net_str = 'net10000';
    end

    for sig_d = sigma_d
        sigma_d_str = ['sigma_d_', strrep(strrep(num2str(sig_d), '.', ''), '-', '')];

        for rr = real_arr
            real_str = ['real', num2str(rr)];
            [noised_net, sp] = create_realization(net_original, sig_d, sigma_a, ball_radius, x0_bound);
            real_struct.(net_str).(sigma_d_str).(real_str).net_noise = noised_net;
            real_struct.(net_str).(sigma_d_str).(real_str).x0 = sp;
        end
    end
end

%% Run Algorithms

out = struct;
if net_type == 1
    num_nets = 1;
end
num_real = length(real_arr);  % number of realizations

for nn = 1:num_nets
    if net_type == 1
        net_original = load('datasets/net1000bench.mat');
        net_str = 'net1000bench';
    elseif net_type == 2
        net_original = load(['datasets/net', num2str(nn), '.mat']);
        net_str = ['net', num2str(nn)];
    elseif net_type == 3
        net_original = load('datasets/net10000.mat');
        net_str = 'net10000';
    end

    for sig_d = sigma_d
        sigma_d_str = ['sigma_d_', strrep(strrep(num2str(sig_d), '.', ''), '-', '')];

        for rr = real_arr
            real_str = ['real', num2str(rr)];
            fprintf('*****Net type=%1d | Net number=%1d | sigma=%10.10f | Realization=%1d*****\n', net_type, nn, sig_d, rr)
            net_noise = real_struct.(net_str).(sigma_d_str).(real_str).net_noise;
            x0 = real_struct.(net_str).(sigma_d_str).(real_str).x0;

            for mm = methods
                met = mm{1};

                switch met
                    case 'SNAM_cen'
                        for ss = 1:length(s_SNAM)
                            s_str = ['s', num2str(s_SNAM(ss))];
                            o = alg_SNAM(net_noise, net_original, x0, r_SNAM, s_SNAM(ss), max_iter, 1);
                            out.(met).(s_str).(sigma_d_str).(net_str).realization{rr, 1} = o;
                        end

                    case 'SNAM_dist'
                        for ss = 1:length(s_SNAM)
                            s_str = ['s', num2str(s_SNAM(ss))];
                            o = alg_SNAM(net_noise, net_original, x0, r_SNAM, s_SNAM(ss), max_iter, 0);
                            out.(met).(s_str).(sigma_d_str).(net_str).realization{rr, 1} = o;
                        end

                    case 'AMFC'
                        o = alg_AMFC(net_original, net_noise, x0, max_iter);
                        out.(met).(sigma_d_str).(net_str).realization{rr, 1} = o;

                    case 'AMFD'
                        o = alg_AMFD(net_noise, x0, max_iter);
                        out.(met).(sigma_d_str).(net_str).realization{rr, 1} = o;

                    case 'SGO'
                        o = alg_SGO(net_original, net_noise, x0, max_iter);
                        out.(met).(sigma_d_str).(net_str).realization{rr, 1} = o;

                    case 'SF'
                        o = alg_SF(net_original, net_noise, x0, max_iter);
                        out.(met).(sigma_d_str).(net_str).realization{rr, 1} = o;

                    case 'SOCP'
                        o = alg_SOCP(net_original, net_noise);
                        out.(met).(sigma_d_str).(net_str).realization{rr, 1} = o;

                    case 'DC'
                        o = alg_DC(net_original, net_noise, x0, max_iter, regul);
                        out.(met).(sigma_d_str).(net_str).realization{rr, 1} = o;
                end
            end
        end
    end
end


% Save the outputs for each method separately
if net_type == 1
    net_str_save = 'net1000bench';
elseif net_type == 2
    net_str_save = 'net1000rand';
elseif net_type == 3
    net_str_save = 'net10000rand';
end
outFolder = "output/" + formattedDateTime + "/" + net_str_save;
mkdir(outFolder)

for mm = methods
    met = mm{1};

    switch met
        case {'SNAM_cen', 'SNAM_dist'}
            for ss = 1:length(s_SNAM)
                s_str = ['s', num2str(s_SNAM(ss))];
                O = out.(met).(s_str);
                %outFolder_met = outFolder + "/" + met + "_" + s_str; mkdir(outFolder_met)
                %save(outFolder_met + "/" + met + "_" + s_str + ".mat", 'O')
            end

        otherwise
            O = out.(met);
            %outFolder_met = outFolder + "/" + met + "_real_" + num2str(min(real_arr)) + "_to_" + num2str(max(real_arr)) ; mkdir(outFolder_met)
            %save(outFolder_met + "/" + met + ".mat", 'O')
    end
end

fprintf('Run Completed!')