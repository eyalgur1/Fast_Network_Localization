sigm = {'0001','00031623','001','0031623'};
xreal=net.Matrices.X_real;
for ss = 1:length(sigm)
    sb = 0;
    sig = ['sigma_d_',sigm{ss}];
    for rl = 1:25
        sb = sb + (1/25)*norm(O.(sig).net1000bench.realization{rl}.location-xreal)^2;
    end
end
sb