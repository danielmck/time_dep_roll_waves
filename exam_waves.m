% dirname = 'channel_roll_wave/results/tau0_0_theta_12_2000';
dirname = 'four_eqn_model/results/tau0_0_theta_12_4001';
dat=hs.Load(dirname);
% hs.Plot(dat,1,[0.0,0.02])

final = dat(end);
theta = final.params.theta;
h0 = 0.1;
lambda = final.xSize;
% rho = final.params.rho;
% chi = final.params.chi;
% g = final.params.g;



if isfield(final.params,"tau0")
    tau0=final.params.tau0;
else
    tau0=0;
end

final_grid = final.xGrid;
final_y = permute(final.data,[3,1,2]);
final_h = final_y(1,:);


[~,ind] = min(final_h);
final_h = horzcat(final_h(ind:end),final_h(1:ind-1));
final_grid = mod((horzcat(final_grid(ind:end),final_grid(1:ind-1))-final_grid(ind)),lambda);
final_y = horzcat(final_y(:,ind:end),final_y(:,1:ind-1));

final_hu = final_y(2,:);
final_u = final_hu./final_h;

if size(final_y,1) > 2
    final_hphi = final_y(3,:);
    final_pbh = final_y(4,:);

    final_phi = final_hphi./final_h;

    if isfield(final.params,"rho")
        rho=final.params.rho;
    else
        rho=final.params.rhog*final_phi+final.params.rhof*(1-final_phi);
    end

    if isfield(final.params,"chi")
        chi=final.params.chi;
    else
        chi=(final.params.rhof+3*rho)./rho/4;
    end
    final_pb = final_pbh./final_h+rho.*g.*cosd(theta).*chi.*final_h;
    final_pp = final_h.*rho.*g.*cosd(theta)-final_pb;
    final_pe = final_pb - final.params.rhof.*g.*cosd(theta).*final_h;
end

ode_denom = final_h.^3*g*cosd(theta)-final_hu.^2;
%%
% [phi_c,rho_f,rho_p,rho,eta_f,g] = get_params_water();
% P = (rho-rho_f)/rho;
% [Fr,eq_Iv] = crit_Iv_tau0_h(theta, rho_p, rho_f, eta_f, h0, tau0,0);
% u0 = Fr*sqrt(g*cosd(theta)*h0);
% 
% a=eq_Iv*dmudIv_fn(eq_Iv);
% b=-eq_Iv/2*dmudIv_fn(eq_Iv);
% 
% omega = 1;
% omega_dl = omega/u0*h0;
% g_omega = 4*omega_dl^2/Fr^2-(a-b)^2*P^2/Fr^4;
% determ = sqrt(-g_omega+sqrt(g_omega^2+16*P^2*omega_dl^2/Fr^4*(a/Fr^2-b)^2));
% k_i_plus = ((a-b)*P/2/(Fr^2-1)+1/(2*sqrt(2)*(Fr^2-1))*determ)/h0;
% k_i_minus = ((a-b)*P/2/(Fr^2-1)-1/(2*sqrt(2)*(Fr^2-1))*determ)/h0;
% % k_i_plus = -eq_Iv*P*(2+1/Fr)/(2*Fr*(1+1/Fr))/h0;
% % k_i_minus = -eq_Iv*P*(2-1/Fr)/(2*Fr*(1-1/Fr))/h0;
% hold on
% plot(final_grid,abs(final_h-0.1))
% % plot(final_grid,h0+0.001*h0*exp(k_i_minus*final_grid))
% plot(final_grid,0.001*h0*exp(k_i_plus*final_grid))
% ylim([0,0.2])
%%

hold on
% SetPaperSize(8,8)
plot(final_grid,final_h,"DisplayName","Variable $\rho$ model")
% plot(final_grid,final.params.rhof.*g.*cosd(theta).*final_h)
ylabel("$h$ ($m$)")
% ylabel("$u$ ($ms^{-1}$)")
% ylabel("$\phi$")
% ylabel("$p_b$ ($Pa$)")
xlabel("$\xi$ (m)")
% % ylim([0,1.8])
legend("Location","south")
title("$\theta = "+num2str(theta)+"^{\circ}$, $\tau_0 = "+num2str(tau0)+"$Pa, $t="+num2str(final.time)+"$s")
% exp_graph(gcf,"rho_vary_comp_"+num2str(theta)+"deg_tau0_"+num2str(tau0)+"_h.pdf")