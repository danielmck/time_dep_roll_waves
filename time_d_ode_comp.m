dirname = 'four_eqn_var_rho/results_ssh/tau0_0_theta_12_6000';
dat=hs.Load(dirname);
h0 = 0.00729731;%;0.0084677798
n_times = [1,3,11,51];%26,51,301,151,201
C = viridis(size(n_times,2));
hold on
for i = 1:size(n_times,2)
    final = dat(n_times(i));
    g = final.params.g;
    theta = final.params.theta;
    lambda = final.xSize;

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
    
    if final.nDims > 3
        final_hphi = final_y(3,:);
        final_phi = final_hphi./final_h;
        rho = final_phi*final.params.rhog+(1-final_phi)*final.params.rhof;
        if isfield(final.params,"chi")
            chi=final.params.chi;
        else
            chi=(final.params.rhof+3*rho)./rho/4;
        end
        final_pbh = final_y(4,:);
        final_pb = final_pbh./final_h+rho.*g.*cosd(theta).*chi.*final_h;
        final_pp = final_h.*rho.*g.*cosd(theta)-final_pb;
        final_pe = final_pb - final.params.rhof.*g.*cosd(theta).*final_h;
    else
        final_phi = final.params.phim;
        final_pb = final.params.rhof.*g.*cosd(theta).*final_h;
        final_pp = (rho-final.params.rhof).*g.*cosd(theta).*final_h;
        final_pe = zeros(size(final_h));
    end

    final_flux = get_flux(final_grid,final_h,final_u);
    [phi_c,rho_f,rho_p,rho,eta_f,g] = get_params_water();
    [Fr_eq,Iv_eq] = crit_Iv_tau0_h(theta, rho_p, rho_f, eta_f, h0, tau0,0,true);
    u0 = Fr_eq*sqrt(g*cosd(theta)*h0);
    base_flux = h0*u0;
    flux_ratio = final_flux/base_flux;
    
    plot(final_grid,final_pb,"DisplayName","$t="+num2str(final.time)+"$s",'color',C(i,:))
end
% hs.Plot(dat,1,100);

%%
% dat_len = size(struct2table(dat),1);
% 
% max_pos = zeros(dat_len,1);
% wave_speed = zeros(dat_len-1,1);
% for i = 1:dat_len
%     curr = dat(i);
%     curr_y = permute(curr.data,[3,1,2]);
%     curr_h = curr_y(1,:);
%     [max_h,ind] = max(curr_h);
%     max_pos(i) = final.xGrid(ind);
%     if i > 1
%         del_t = dat(i).time-dat(i-1).time;
%         wave_speed(i-1) = mod(max_pos(i) - max_pos(i-1),lambda)/del_t;
%     end
% end
% plot(wave_speed)
%%

remake_ode = false;
ode_file = [dirname, '/ode_orig.txt'];
flux_file = [dirname, '/ode_check.txt'];
% ode_file = [dirname, '/ode_pres_h.txt'];
if isfile(ode_file) && ~remake_ode
    in_comp = load(ode_file);
    xi_comp = in_comp(1,:);
    y_comp = in_comp(2:end,:);
else   
    [xi_comp,y_comp] = time_dep_recreate(h0,theta,lambda,tau0);
    out_comp = vertcat(xi_comp,y_comp);
    save(ode_file,"out_comp","-ascii")
end
% if isfile(flux_file) && ~remake_ode
%     in_flux = load(flux_file);
%     xi_flux = in_flux(1,:);
%     y_flux = in_flux(2:end,:);
% else
%     [xi_flux,y_flux] = time_dep_recreate(h0,theta,lambda,tau0,flux_ratio);
%     out_flux = vertcat(xi_flux,y_flux);
%     save(flux_file,"out_flux","-ascii")
% end
% if isfile(pres_h_file) && ~remake_ode
%     in_pres_h = load(pres_h_file);
%     xi_pres_h = in_pres_h(1,:);
%     y_pres_h = in_pres_h(2:end,:);
% else
%     [xi_pres_h,y_pres_h] = time_dep_recreate(h0,theta,lambda,tau0,1,1);
%     out_pres_h = vertcat(xi_pres_h,y_pres_h);
%     save(pres_h_file,"out_pres_h","-ascii")
% end

% params = load([dirname, '/ode_params.txt']);

h_comp = y_comp(4,:).*h0;
u_w_comp = y_comp(1,1).*u0;
Q1_comp = y_comp(2,:).*u0*h0;
xi_comp = y_comp(3,1)*xi_comp*h0;
u_comp = u_w_comp - Q1_comp./h_comp;
flux_comp_ode = y_comp(5,:);
flux_comp_2 = get_flux(xi_comp,h_comp,u_comp);
if size(y_comp,1) > 5
    phi_comp = y_comp(6,:).*u0*h0./Q1_comp;
    rho_comp = rho_p*phi_comp+rho_f*(1-phi_comp);
    chi_comp = (rho_f+3*rho_comp)./(4*rho_comp);
    pb_comp = y_comp(7,:).*rho_f*g*cosd(theta)*h0 + rho_comp.*g*cosd(theta).*chi_comp.*h_comp;
end

% h_flux = y_flux(4,:).*h0;
% u_w_flux = y_flux(1,1).*u0;
% Q1_flux = y_flux(2,:).*u0*h0;
% xi_flux = y_flux(3,1)*xi_flux*h0;
% u_flux = u_w_flux - Q1_flux./h_flux;
% flux_flux_ode = y_flux(5,:);
% flux_flux_2 = get_flux(xi_flux,h_flux,u_flux);
% if size(y_flux,1) > 5
%     phi_flux = y_flux(6,:).*u0*h0./Q1_flux;
%     rho_flux = rho_p*phi_flux+rho_f*(1-phi_flux);
%     chi_flux = (rho_f+3*rho_flux)./(4*rho_flux);
%     pb_flux = y_flux(7,:).*rho_f*g*cosd(theta)*h0 + rho_flux.*g*cosd(theta).*chi_flux.*h_flux;
% end

% h_flux = y_flux(3,:);
% u_w_flux = y_flux(1,1);
% Q1_flux = y_flux(2,1);
% u_flux = u_w_flux - Q1_flux./h_flux;

hold on
SetPaperSize(7.6,7.6)
% time_d_u_comp = interp1(final_grid,final_u,xi_comp);
% plot(final_grid,final_h,"DisplayName","ODE solution")

plot(xi_comp,pb_comp,'r--',"DisplayName","ODE solution")
% plot(xi_flux,u_flux,"DisplayName","Flux Matched ODE solution")
% ylabel("$h$ ($m$)")
% ylabel("$u$ ($ms^{-1}$)")
% ylabel("$\phi$")
ylabel("$p_b$ ($Pa$)")
xlabel("$\xi$ ($m$)")
% ylim([0,1.8])
% legend("Location","best")
% ax = gca;
% ax.YAxis.Exponent = 0;
title("$Fr = 0.8$, $\theta = "+num2str(theta)+"^{\circ}$, $\tau_0 = "+num2str(tau0)+"Pa$")
exp_graph(gcf,"report_time_conv_full_pb.pdf")