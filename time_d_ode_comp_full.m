dirname = 'four_eqn_var_rho/results/tau0_0_theta_12_2500';
dat=hs.Load(dirname);
final = dat(end);
theta = final.params.theta;
lambda = final.xSize;
% rho = final.params.rho;
% chi = final.params.chi;
rho_p = final.params.rhog;
rho_f = final.params.rhof;
g = final.params.g;
d = final.params.d;
alpha = final.params.alpha;

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
final_hphi = final_y(3,:);
final_pbh = final_y(4,:);

final_u = final_hu./final_h;
final_phi = final_hphi./final_h;
rho = rho_p*final_phi+rho_f*(1-final_phi);
chi = (3*rho+rho_f)./rho;
final_pb = final_pbh./final_h+rho.*g*cosd(theta).*chi.*final_h;

final_flux = get_flux(final_grid,final_h,final_u);
[phi_c,rho_f,rho_p,rho,eta_f,g] = get_params_water();
[Fr_eq,Iv_eq] = crit_Iv_tau0_h(theta, rho_p, rho_f, eta_f, h0, tau0,0, true);
u0 = Fr_eq*sqrt(g*cosd(theta)*h0);
base_flux = h0*u0;
flux_ratio = final_flux/base_flux;
% hs.Plot(dat,1,100);

remake_ode = false;
pres_h_file = [dirname, '/ode_pres_h.txt'];
if isfile(pres_h_file) && ~remake_ode
    in_pres_h = load(pres_h_file);
    xi_pres_h = in_pres_h(1,:);
    y_pres_h = in_pres_h(2:end,:);
    param_file = load([dirname, '/ode_params.txt']);
    
    Fr_pres_h = param_file(1);
    theta_pres_h = param_file(2);
    d_pres_h = param_file(4);
    alpha_pres_h = param_file(3);
    tau0_pres_h = param_file(5);
    [h0_pres_h,Iv_pres_h] = crit_Iv_tau0(theta_pres_h, rho_p, rho_f, eta_f, Fr_pres_h, tau0_pres_h,false,true);
    u0_pres_h = Fr_pres_h*sqrt(g*cosd(theta_pres_h)*h0_pres_h);
    
    u_w_pres_h = param_file(6)*u0_pres_h;
    lambda_pres_h = param_file(7);
    crit_xi = param_file(8);
    xi_pres_h = horzcat(xi_pres_h(xi_pres_h<1)*crit_xi,crit_xi+(xi_pres_h(xi_pres_h>=1)-1)*(lambda/h0-crit_xi));
    xi_pres_h = xi_pres_h*h0_pres_h;
else
    [xi_pres_h,y_pres_h] = time_dep_recreate_full(h0,theta,lambda,tau0,d,alpha,1,1);
    out_pres_h = vertcat(xi_pres_h,y_pres_h);
    save(pres_h_file,"out_pres_h","-ascii")
end

% h_comp = y_comp(3,:);
% u_w_comp = y_comp(1,1);
% Q1_comp = y_comp(2,1);
% u_comp = u_w_comp - Q1_comp./h_comp;
% flux_comp_ode = y_comp(4,:);
% flux_comp_2 = get_flux(xi_comp,h_comp,u_comp);
% 
% h_flux = y_flux(3,:);
% u_w_flux = y_flux(1,1);
% Q1_flux = y_flux(2,1);
% u_flux = u_w_flux - Q1_flux./h_flux;

h_pres_h = y_pres_h(2,:)*h0;
% u_w_pres_h = y_pres_h(1,1);
Q1_pres_h = y_pres_h(1,:)*h0*u0;
u_pres_h = (u_w_pres_h - Q1_pres_h./h_pres_h);
phi_pres_h = h0*u0*y_pres_h(4,:)./Q1_pres_h;
rho_pres_h = rho_p*phi_pres_h+rho_f*(1-phi_pres_h);
chi_pres_h = (3*rho_pres_h+rho_f)./rho_pres_h;
pb_pres_h = rho_f*g*cosd(theta)*h0*y_pres_h(5,:) + rho_pres_h.*g*cosd(theta).*chi_pres_h.*h_pres_h;

%%

% vis_dirname = 'four_eqn_model_vis/results/tau0_0_theta_12_1000';
% vis_dat=hs.Load(vis_dirname);
% 
% vis_final = vis_dat(101);
% vis_lambda = vis_final.xSize;
% 
% vis_final_grid = vis_final.xGrid;
% vis_final_y = permute(vis_final.data,[3,1,2]);
% 
% vis_final_h = vis_final_y(1,:);
% [~,ind] = min(vis_final_h);
% vis_final_h = horzcat(vis_final_h(ind:end),vis_final_h(1:ind-1));
% vis_final_grid = mod((horzcat(vis_final_grid(ind:end),vis_final_grid(1:ind-1))-vis_final_grid(ind)),vis_lambda);
% vis_final_y = horzcat(vis_final_y(:,ind:end),vis_final_y(:,1:ind-1));
%%

hold on
% SetPaperSize(8,8)
% time_d_u_comp = interp1(final_grid,final_u,xi_comp);
% plot(xi_comp,time_d_u_comp-u_comp,"DisplayName","ODE solution")
plot(final_grid,final_pb,"DisplayName","Time dependent, $t="+num2str(final.time)+"$s")
% plot(vis_final_grid,vis_final_h,"DisplayName","Viscous Time dep. solution, $t="+num2str(final.time)+"$s")
plot(xi_pres_h,pb_pres_h,"DisplayName","ODE")
% plot(xi_flux,u_flux,"DisplayName","Flux Matched ODE solution")
% ylabel("$h$ ($m$)")
% ylabel("$u$ ($ms^{-1}$)")
% ylabel("$\phi$")
ylabel("$p_b$ ($Pa$)")
xlabel("$x$ (m)")
% % ylim([0,1.8])
% legend("Location","best")
title("$\theta = "+num2str(theta)+"^{\circ}$, $\tau_0 = "+num2str(tau0)+"$Pa")
% exp_graph(gcf,"full_wave_comp_"+num2str(theta)+"deg_tau_"+num2str(tau0)+"_pb.pdf")