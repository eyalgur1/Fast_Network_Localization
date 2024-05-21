%% Plot Average Function Values and Average Norm of Bias

num_nets = 1;
methods = {'SNAM_cen','SNAM_dist','SGO','SOCP','DC'};
max_iter = 10000;
s_SNAM = [40,80,120];
num_real = 25;
net_type = 1;

%%

sigma_d_intervals = 6;
sigma_d_lb = 0.02;  % sets sigma_d to lb% to ub% of shortest edge
sigma_d_ub = 0.5;
if net_type == 1
    %sigma_d = logspace(log10(0.001*sigma_d_lb), log10(0.001*sigma_d_ub), sigma_d_intervals);
    sigma_d = logspace(-3, -1.5, 4);
elseif net_type == 2
    %sigma_d = logspace(log10(0.000431*sigma_d_lb), log10(0.000431*sigma_d_ub), sigma_d_intervals);
    sigma_d = logspace(-5, -1.5, sigma_d_intervals);
elseif net_type == 3
    %sigma_d = logspace(log10(0.000431*sigma_d_lb), log10(0.000431*sigma_d_ub), sigma_d_intervals);
    sigma_d = logspace(-5, -1.5, sigma_d_intervals);
end

if net_type == 1
    %net_str_save = 'net1000bench\';
    outFolder = "output\" + "Benchmark1000_complete_data\";% + "\" + net_str_save;
    outFolder = "output\2024_01_04_16_25_34\net1000bench\";
elseif net_type == 2
    %net_str_save = 'net1000rand\';
    outFolder = "output\" + "Rand1000_complete_data\";% + "\" + net_str_save;
elseif net_type == 3
    %net_str_save = 'net1000rand\';
    outFolder = "output\" + "Rand10000_complete_data\";% + "\" + net_str_save;
end
%%
%Set structure for averages with 0s
avgs = struct;
avgs.net_type = net_type;
for mm = 1:length(methods)
    met = methods{mm};

    for vv = 1:length(sigma_d)
        sig_d = sigma_d(vv);
        sigma_d_str = ['sigma_d_', strrep(strrep(num2str(sig_d), '.', ''), '-', '')];

        switch met
            case {'SNAM_cen', 'SNAM_dist'}
                for ss = 1:length(s_SNAM)
                    avgs.(met).(['s', num2str(s_SNAM(ss))]).(sigma_d_str).avg_fun_val = 0;
                    avgs.(met).(['s', num2str(s_SNAM(ss))]).(sigma_d_str).avg_bias = 0;
                    avgs.(met).(['s', num2str(s_SNAM(ss))]).(sigma_d_str).location = 0;
                end

            otherwise
                avgs.(met).(sigma_d_str).avg_fun_val = 0;
                avgs.(met).(sigma_d_str).avg_bias = 0;
                avgs.(met).(sigma_d_str).location = 0;
        end

    end
end

%
%Calculate the averages
if net_type == 1 || net_type == 3
    avg_factor = 1/(num_real);
    num_nets = 1;
else
    avg_factor = 1/(num_real+num_nets);
end

for nn = 1:num_nets
    if net_type == 1
        net_str = 'net1000bench';
    elseif net_type == 2
        net_str = ['net', num2str(nn)];
    elseif net_type == 3
        net_str = 'net10000';
    end

    for vv = 1:length(sigma_d)
        sig_d = sigma_d(vv);
        sigmad_d_str = ['sigma_d_', strrep(strrep(num2str(sig_d), '.', ''), '-', '')];

        for rr = 1:num_real

            for mm = 1:length(methods)
                met = methods{mm};
                fprintf(['net=%1d | sigma=%1d | real=%1d | ', met, '\n'], nn, vv, rr)

                switch met
                    case {'SNAM_cen', 'SNAM_dist'}
                        for ss = 1:length(s_SNAM)
                            s_str = ['s', num2str(s_SNAM(ss))];
                            %                             if s_SNAM(ss) == 40 && contains(met, 'cen')
                            %                                 result = result_SNAM_cen_s40;
                            %                             elseif s_SNAM(ss) == 40 && contains(met, 'dist')
                            %                                 result = result_SNAM_dist_s40;
                            %                             elseif s_SNAM(ss) == 80 && contains(met, 'cen')
                            %                                 result = result_SNAM_cen_s80;
                            %                             elseif s_SNAM(ss) == 80 && contains(met, 'dist')
                            %                                 result = result_SNAM_dist_s80;
                            %                             end
                            result = load(outFolder + met + "_" + s_str + "\" + met + "_" + s_str +".mat");
                            avgs.(met).(s_str).(sigmad_d_str).avg_fun_val = avgs.(met).(s_str).(sigmad_d_str).avg_fun_val + avg_factor*result.O.(sigmad_d_str).(net_str).realization{rr, 1}.fun_val;
                            avgs.(met).(s_str).(sigmad_d_str).avg_bias = avgs.(met).(s_str).(sigmad_d_str).avg_bias + avg_factor*result.O.(sigmad_d_str).(net_str).realization{rr, 1}.norm_bias;
                            avgs.(met).(s_str).(sigmad_d_str).location = avgs.(met).(s_str).(sigmad_d_str).location + avg_factor*result.O.(sigmad_d_str).(net_str).realization{rr, 1}.location;
                        end

                    otherwise
                        result = load(outFolder + met + "\" + met + ".mat");
                        %                         if contains(met, 'AMFC')
                        %                             result = result_AMFC;
                        %                         elseif contains(met, 'SGO')
                        %                             result = result_SGO;
                        %                         end
                        avgs.(met).(sigmad_d_str).avg_fun_val = avgs.(met).(sigmad_d_str).avg_fun_val + avg_factor*result.O.(sigmad_d_str).(net_str).realization{rr, 1}.fun_val;
                        avgs.(met).(sigmad_d_str).avg_bias = avgs.(met).(sigmad_d_str).avg_bias + avg_factor*result.O.(sigmad_d_str).(net_str).realization{rr, 1}.norm_bias;
                        avgs.(met).(sigmad_d_str).location = avgs.(met).(sigmad_d_str).location + avg_factor*result.O.(sigmad_d_str).(net_str).realization{rr, 1}.location;
                end
            end
        end
    end
end

%%
% Set legend for plots
names_to_legend = {};
for mm = 1:length(methods)
    met = methods{mm};
    switch met
        case {'SNAM_cen', 'SNAM_dist'}
            for ss = 1:length(s_SNAM)
                switch met
                    case 'SNAM_cen'
                        names_to_legend{end + 1} = ['$\textrm{SNAM,\ centralized,\ }s=', num2str(s_SNAM(ss)), '$'];
                    case 'SNAM_dist'
                        names_to_legend{end + 1} = ['$\textrm{SNAM,\ distributed,\ }s=', num2str(s_SNAM(ss)), '$'];
                end
            end

        otherwise
            names_to_legend{end + 1} = ['$\textrm{', met,'}$'];
    end
end

%%
% Plot function values and bias
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultLegendInterpreter' ,'latex');
set(groot, 'defaultLineLineWidth',1)
set(groot, 'defaultAxesFontSize',14)
colors = distinguishable_colors(10);

if net_type == 1
    net_str_save = 'Benchmark1000_complete_data\';
elseif net_type == 2
    net_str_save = 'Rand1000_complete_data\';
else
    net_str_save = 'Rand10000_complete_data\';
end
load(['output\', net_str_save, 'avgs.mat'])
outFigFolder = "output\" + net_str_save + "\figs";
% outFigFolder = "output\" + formattedDateTime + "\" + net_str_save + "\figs";
%mkdir(outFigFolder)
%%
for vv = 1:length(sigma_d)
    sig_d = sigma_d(vv);
    figure(vv); hold on; figure(100 + vv); hold on
    sigma_d_str = ['sigma_d_', strrep(strrep(num2str(sig_d), '.', ''), '-', '')];
    sigma_d_title = ['\textrm{', num2str(sig_d), '}'];

    for mm = 1:length(methods)
        met = methods{mm};

        switch met
            case 'SNAM_cen'
                for ss = 1:length(s_SNAM)
                    if s_SNAM(ss) == 40
                        figure(vv); semilogy(avgs.(met).(['s', num2str(s_SNAM(ss))]).(sigma_d_str).avg_fun_val, '--' , 'Color', colors(1, :));
                        figure(100 + vv); semilogy(avgs.(met).(['s', num2str(s_SNAM(ss))]).(sigma_d_str).avg_bias, '--', 'Color', colors(1, :));
                    elseif s_SNAM(ss) == 80
                        figure(vv); semilogy(avgs.(met).(['s', num2str(s_SNAM(ss))]).(sigma_d_str).avg_fun_val, '*','Color', colors(1, :));
                        figure(100 + vv); semilogy(avgs.(met).(['s', num2str(s_SNAM(ss))]).(sigma_d_str).avg_bias, '*','Color', colors(1, :));
                    elseif s_SNAM(ss) == 120
                        figure(vv); semilogy(avgs.(met).(['s', num2str(s_SNAM(ss))]).(sigma_d_str).avg_fun_val, 'Color', colors(1, :));
                        figure(100 + vv); semilogy(avgs.(met).(['s', num2str(s_SNAM(ss))]).(sigma_d_str).avg_bias, 'Color', colors(1, :));
                    end
                end

            case 'SNAM_dist'
                for ss = 1:length(s_SNAM)
                    if s_SNAM(ss) == 40
                        figure(vv); semilogy(avgs.(met).(['s', num2str(s_SNAM(ss))]).(sigma_d_str).avg_fun_val, '--', 'Color', colors(2, :));
                        figure(100 + vv); semilogy(avgs.(met).(['s', num2str(s_SNAM(ss))]).(sigma_d_str).avg_bias, '--', 'Color', colors(2, :));
                    elseif s_SNAM(ss) == 80
                        figure(vv); semilogy(avgs.(met).(['s', num2str(s_SNAM(ss))]).(sigma_d_str).avg_fun_val, '*','Color', colors(2, :));
                        figure(100 + vv); semilogy(avgs.(met).(['s', num2str(s_SNAM(ss))]).(sigma_d_str).avg_bias, '*','Color', colors(2, :));
                    elseif s_SNAM(ss) == 120
                        figure(vv); semilogy(avgs.(met).(['s', num2str(s_SNAM(ss))]).(sigma_d_str).avg_fun_val, 'Color', colors(2, :));
                        figure(100 + vv); semilogy(avgs.(met).(['s', num2str(s_SNAM(ss))]).(sigma_d_str).avg_bias, 'Color', colors(2, :));
                    end
                end

            case 'AMFC'
                figure(vv); semilogy(avgs.(met).(sigma_d_str).avg_fun_val, 'Color', colors(3, :));
                figure(100 + vv); semilogy(avgs.(met).(sigma_d_str).avg_bias, 'Color', colors(3, :));

            case 'AMFD'
                figure(vv); semilogy(avgs.(met).(sigma_d_str).avg_fun_val, 'Color', colors(4, :));
                figure(100 + vv); semilogy(avgs.(met).(sigma_d_str).avg_bias, 'Color', colors(4, :));

            case 'SGO'
                figure(vv); semilogy(avgs.(met).(sigma_d_str).avg_fun_val, 'Color', colors(5, :));
                figure(100 + vv); semilogy(avgs.(met).(sigma_d_str).avg_bias, 'Color', colors(5, :));

            case 'SF'
                figure(vv); semilogy(avgs.(met).(sigma_d_str).avg_fun_val, 'Color', colors(6, :));
                figure(100 + vv); semilogy(avgs.(met).(sigma_d_str).avg_bias, 'Color', colors(6, :));

            case 'SOCP'
                figure(vv); semilogy(avgs.(met).(sigma_d_str).avg_fun_val*ones(1, max_iter + 1), 'Color', colors(7, :));
                figure(100 + vv); semilogy(avgs.(met).(sigma_d_str).avg_bias*ones(1, max_iter + 1), 'Color', colors(7, :));

            case 'DC'
                figure(vv); semilogy(avgs.(met).(sigma_d_str).avg_fun_val*ones(1, max_iter + 1), 'Color', colors(8, :));
                figure(100 + vv); semilogy(avgs.(met).(sigma_d_str).avg_bias*ones(1, max_iter + 1), 'Color', colors(8, :));

        end
    end
    figure(vv); set(gca, 'YScale', 'log')
    xlim([0 max_iter + 1]); xlabel('$\textrm{Iterations}$', 'Interpreter', 'latex', 'FontSize', 14)
    ylabel('$\textrm{ML\ Function\ Value}$', 'Interpreter', 'latex', 'FontSize', 14)
    title(['$\textrm{Average\ ML\ Function\ Value\ for\ }\sigma_{d}=', sigma_d_title, '$'], 'Interpreter', 'latex', 'FontSize', 14)
    legend(names_to_legend)
    %savefig(outFigFolder + "\funv_" + net_str + num2str(sig_d) + ".fig")
    hold off;

    figure(100 + vv); set(gca, 'YScale', 'log')
    xlim([0 max_iter + 1]); xlabel('$\textrm{Iterations}$', 'Interpreter', 'latex', 'FontSize', 14)
    ylabel('$\textrm{ML\ Norm\ of\ Bias}$', 'Interpreter', 'latex', 'FontSize', 14)
    title(['$\textrm{Average\ Norm\ of\ Bias\ for\ }\sigma_{d}=', sigma_d_title, '$'], 'Interpreter', 'latex', 'FontSize', 14)
    legend(names_to_legend)
    %savefig(outFigFolder + "\bias_" + net_str + "_" + num2str(sig_d) + ".fig")
    hold off;
end
