function [t_vals, ave] = average_waves(dat)
    ntimes = size(dat,2);
    ave = zeros(1,ntimes);
    t_vals = zeros(1,ntimes);
    for i = 1:ntimes
        final = dat(i);
        t_vals(i) = final.time;
        final_grid = final.xGrid;
        lambda = final_grid(end);
        npts = size(final_grid,2);
        g = final.params.g;
        theta = final.params.theta;
        final_y = permute(final.data,[3,1,2]);
        final_h = final_y(1,:);
        final_hu = final_y(2,:);
        final_u = final_hu./final_h;
        final_hphi = final_y(3,:);
        final_pbh = final_y(4,:);
        final_phi = final_hphi./final_h;
        final_rho = final.params.rhog*final_phi+final.params.rhof*(1-final_phi);
        final_flux = final_rho.*final_u.*final_h;
        final_chi = (final.params.rhof+3*final_rho)./final_rho/4;
        final_pb = final_pbh./final_h+final_rho.*g.*cosd(theta).*final_chi.*final_h;
        final_pp = final_h.*final_rho.*g.*cosd(theta)-final_pb;
        final_pe = final_pb - final.params.rhof.*g.*cosd(theta).*final_h;
        final_psi1 = final_phi./final_rho./final_h.*final_pe;
        sin_ave = 0;
        for j = 2:npts
            sin_ave = sin_ave + (final_flux(j)+final_flux(j-1))/2*(final_grid(j)-final_grid(j-1))/lambda;
        end
        ave(i) = sin_ave;
    end
end