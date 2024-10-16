function slow_comp
[phi_c,rho_f,rho_p,~,eta_f,g] = get_params_water();
theta = 5;
dirname_st = 'four_eqn_var_rho/results/';
n_runs = 14;
npts = 1001;
flux_vals = zeros(n_runs,npts);
for i = 1:n_runs
    lambda = 20+(i-1)*10;
    exact_dir = ['lambda_',num2str(lambda),'_tau0_0_theta_5_4000'];
    dirname = [dirname_st,exact_dir];
    dat=hs.Load(dirname);
    for j=1:min(npts,size(dat,2))
        final = dat(j);
        final_y = permute(final.data,[3,1,2]);
        final_grid = final.xGrid;
        wave_len = final_grid(end);
        final_h = final_y(1,:);
        final_hu = final_y(2,:);
        final_u = final_hu./final_h;
        final_hphi = final_y(3,:);
        final_phi = final_hphi./final_h;
        final_pbh = final_y(4,:);
        rho = rho_p*final_phi+rho_f*(1-final_phi);
        chi = (rho_f+3*rho)./rho/4;
        final_pb = final_pbh./final_h+rho.*g.*cosd(theta).*chi.*final_h;
        final_pe = final_pb - final.params.rhof.*g.*cosd(theta).*final_h;
        LR = (final_pe./(rho.*g.*cosd(theta).*final_h));
        num_static = sum(final_u<1e-5);
        flux = num_static/size(final_grid,2);
%         flux = min(num_static);
%         if size(flux,2) == 0
%             flux = NaN;
%         end
        for k=2:size(final_grid,2)
            flux = flux + (LR(k)+LR(k-1))/2*(final_grid(k)-final_grid(k-1))/wave_len;
        end
        flux_vals(i,j) = flux;
    end
    
end
save("four_eqn_var_rho/flux_decrease_5deg_LR_ave_fine.txt","flux_vals","-ascii")

end