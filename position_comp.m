function position_comp
[phi_c,rho_f,rho_p,~,eta_f,g] = get_params_water();
theta = 5;
dirname_st = 'ive_rep_inflow/results/tau0_0_theta_31_9500';
n_runs = 14;
npts = 500;
flux_vals = zeros(1,npts);
dat = hs.Load(dirname_st);
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
%     rho = rho_p*final_phi+rho_f*(1-final_phi);
%     chi = (rho_f+3*rho)./rho/4;
%     final_pb = final_pbh./final_h+rho.*g.*cosd(theta).*chi.*final_h;
%     final_pe = final_pb - final.params.rhof.*g.*cosd(theta).*final_h;
%     LR = (final_pe./(rho.*g.*cosd(theta).*final_h));
%     num_static = sum(final_u<1e-5);
%     flux = num_static/size(final_grid,2);
%         flux = min(num_static);
%         if size(flux,2) == 0
%             flux = NaN;
%         end
%     for k=2:size(final_grid,2)
%         flux = flux + (LR(k)+LR(k-1))/2*(final_grid(k)-final_grid(k-1))/wave_len;
%     end
    flux_vals(1,j) = final_h(700);
end
save("ive_rep_inflow/height_val_2m.txt","flux_vals","-ascii")

end