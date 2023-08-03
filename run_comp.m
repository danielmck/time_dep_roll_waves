main_name = 'four_eqn_var_rho/results_ssh/'; %';finalTheta_8.2_initTheta_10_tau0_0_3500

% hs.Plot(dat,1,[0.0,0.2])

run_names = {'tau0_0_theta_12_6000','limiter_MinMod2_tau0_0_theta_12_6000'}; %
times = 76*ones(1,size(run_names,2));
C = viridis(size(run_names,2));
plot_names = ["WENO", "MinMod2"]; %["750 points","1500 points","3000 points","6000 points"];
hold on
for i = 1:size(run_names,2)
    dat=hs.Load([main_name,run_names{i}]);
    t_now = times(i);
    final = dat(t_now);
    theta = final.params.theta;
    lambda = final.xSize;
    % rho = final.params.rho;
    % chi = final.params.chi;
    g = final.params.g;
    rho_p = final.params.rhog;
    rho_f = final.params.rhof;
    % eta_f = final.params.etaf;
    phi_c = final.params.phim;

    if isfield(final.params,"d")
        d=final.params.d;
    else
        d=1e-4;
    end

    if isfield(final.params,"alpha")
        alpha=final.params.alpha;
    else
        alpha=1e-5;
    end
    if isfield(final.params,"etaf")
        eta_f=final.params.etaf;
    else
        eta_f=1.0013e-3;
    end

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

    if final.nDims > 2
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
    %     ppval = rho.*g.*cosd(theta).*final_h-final_pb;
    else
        rho=final.params.rho;
        final_phi = final.params.phim;
        final_pb = final.params.rhof.*g.*cosd(theta).*final_h;
        final_pp = (rho-final.params.rhof).*g.*cosd(theta).*final_h;
        final_pe = zeros(size(final_h));
    end
    plot(final_grid,final_h,"DisplayName",plot_names(i),"color",C(i,:))
end

%%
pres_h_file = '/ode_pres_h.txt';
in_pres_h = load([main_name,run_names{i}, pres_h_file]);
xi_pres_h = in_pres_h(1,:);
y_pres_h = in_pres_h(2:end,:);
param_file = load([main_name, run_names{i},'/ode_params.txt']); %

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
xi_pres_h = horzcat(xi_pres_h(xi_pres_h<1)*crit_xi,crit_xi+(xi_pres_h(xi_pres_h>=1)-1)*(lambda/h0_pres_h-crit_xi));
xi_pres_h = xi_pres_h*h0_pres_h;

h_pres_h = y_pres_h(2,:)*h0_pres_h;
% u_w_pres_h = y_pres_h(1,1);
Q1_pres_h = y_pres_h(1,:)*h0_pres_h*u0_pres_h;
u_pres_h = (u_w_pres_h - Q1_pres_h./h_pres_h);
phi_pres_h = h0_pres_h*u0_pres_h*y_pres_h(4,:)./Q1_pres_h;
rho_pres_h = rho_p*phi_pres_h+rho_f*(1-phi_pres_h);
chi_pres_h = (3*rho_pres_h+rho_f)./rho_pres_h/4;
pb_pres_h = rho_f*g*cosd(theta)*h0_pres_h*y_pres_h(5,:) + rho_pres_h.*g*cosd(theta).*chi_pres_h.*h_pres_h;
pp_pres_h = rho_pres_h.*g.*cosd(theta).*h_pres_h-pb_pres_h;

plot(xi_pres_h, h_pres_h, "DisplayName", "ODE Solution" ,'color','r')

%%
SetPaperSize(10,10)
ylabel("$h$ ($m$)")
% ylabel("$u$ ($ms^{-1}$)")
% ylabel("$\phi$")
% ylabel("$p_b$ ($Pa$)")
xlabel("$x$ (m)")
xlim([0.01,0.012])
ylim([5.83e-3,5.88e-3])
legend("Location","best")
title("$\theta = "+num2str(theta)+"^{\circ}$, $\tau_0 = "+num2str(tau0)+"$Pa, $t="+num2str(final.time)+"$s") %,
exp_graph(gcf,"limiter_comp_h_ode_zoom.pdf")