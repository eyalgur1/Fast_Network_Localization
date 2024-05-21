sigma_d_intervals = 6;
%sigma_d_lb = 0.02;  % sets sigma_d to lb% to ub% of shortest edge
%sigma_d_ub = 0.5;
sigma_d = logspace(-3, -1, 5);
reals = [1:10:30;10:10:30];

for j=1:length(reals)
    ind = reals(:, j);
    folder_str = ['SF_real_',num2str(ind(1)),'_to_',num2str(ind(2))];
    load([folder_str,'\SF.mat']);

    for sig_d = sigma_d
        sigma_d_str = ['sigma_d_', strrep(strrep(num2str(sig_d), '.', ''), '-', '')];

        %for net = 1:30
            net_str = 'net1000bench';%['net',num2str(net)];
            for r = ind(1):ind(2)
                o.(sigma_d_str).(net_str).realization{r,1} = O.(sigma_d_str).(net_str).realization{r,1};
            end
        %end
    end
end