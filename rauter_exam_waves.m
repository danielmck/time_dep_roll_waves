dirname = 'rauter_closure/results/tau0_0_theta_12_3000'; %';finalTheta_8.2_initTheta_10_tau0_0_3500
dat=hs.Load(dirname);
% hs.Plot(dat,1,[0.0,0.05])

final = dat(end);
theta = final.params.theta;
lambda = final.xSize;
% rho = final.params.rho;
% chi = final.params.chi;
g = final.params.g;
rho_p = final.params.rhog;
rho_f = final.params.rhof;
% eta_f = final.params.etaf;
phi_c = final.params.phim;
phi_rlp = final.params.phi_rlp;
phi_rcp = final.params.phi_rcp;
a = final.params.a;

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
    final_phi = final_hphi./final_h;

    if isfield(final.params,"rho")
        rho=final.params.rho;
    else
        rho=final.params.rhog*final_phi+final.params.rhof*(1-final_phi);
    end
%     ppval = rho.*g.*cosd(theta).*final_h-final_pb;
else
    rho=final.params.rho;
    final_phi = final.params.phim;
    final_pb = final.params.rhof.*g.*cosd(theta).*final_h;
    final_pp = (rho-final.params.rhof).*g.*cosd(theta).*final_h;
    final_pe = zeros(size(final_h));
end
iv_phi = (phi_c./final_phi-1).^2;
pp_contact = max(a.*(final_phi-phi_rlp)./(phi_rcp-final_phi),0);
pp_shear = 3*final_u*eta_f./iv_phi;
final_pp = pp_contact+pp_shear;
final_pb = rho.*g.*cosd(theta).*final_h-final_pp;
final_pe = final_pb-rho_f.*g.*cosd(theta).*final_h;
P = (rho-rho_f)./rho;


% if (absu == 0)
%     signU = 1;
% else 
%     signU = uval/absu;
% end

beta = 150*phi_c^2*eta_f/((1-phi_c)^3*d^2);

iv = (3*eta_f).*(abs(final_u)./final_h)./final_pp;

iv_use = max(iv,1e-7);
ppval_use = max(1e-7,final_pp);
pb_use = rho.*g.*cosd(theta).*final_h-ppval_use;

mubf = ppval_use.*mu_Iv_fn(iv_use);

D = -2/beta./final_h.*(pb_use-rho_f*g*cosd(theta)*final_h);

absFriction = (1./rho).*(-mubf - tau0 + (rho-rho_f).*D.*final_u);
force_bal = final_h*g*sind(theta)+absFriction;

%%

hold on
% figure(2)
% SetPaperSize(8,8)
plot(final_grid,final_h,"DisplayName","Wave Profile")
% plot(final_grid(end),pb_max,"x","DisplayName","Shock condition maximum")
% plot(final_grid,final.params.rhof.*g.*cosd(theta).*final_h)
% ylabel("$h$ ($m$)")
% ylabel("$u$ ($ms^{-1}$)")
% ylabel("$\phi$")
ylabel("$p_b$ ($Pa$)")
xlabel("$\xi$ (m)")
% ylim([30,205])
legend("Location","best")
title("$\theta = "+num2str(theta)+"^{\circ}$, $\tau_0 = "+num2str(tau0)+"$Pa, $t="+num2str(final.time)+"$s")
% exp_graph(gcf,"rho_con_"+num2str(theta)+"deg_tau0_"+num2str(tau0)+"_shock_pb.pdf")