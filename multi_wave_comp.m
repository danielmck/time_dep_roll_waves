% Compares the flux of a range of different simulations
flux_vals = zeros(50,301);
for i=11:50
    lambda=10*i;
    dirname = ['four_eqn_var_rho/results_ssh/lambda_',num2str(lambda),'_tau0_0_theta_12_2500'];
    dat=hs.Load(dirname);
    [t_vals, flux_ave] = average_waves(dat);
    flux_vals(i,1:size(flux_ave,2)) = flux_ave;
end