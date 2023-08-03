td_dirname = 'rauter_closure/results/tau0_0_theta_12_3000'; %';finalTheta_8.2_initTheta_10_tau0_0_3500
dat=hs.Load(td_dirname);
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
final_P = (rho-rho_f)./rho;


% if (absu == 0)
%     signU = 1;
% else 
%     signU = uval/absu;
% end

beta = 150*phi_c^2*eta_f/((1-phi_c)^3*d^2);

td_iv = (3*eta_f).*(abs(final_u)./final_h)./final_pp;

td_iv_use = max(td_iv,1e-7);
td_ppval_use = max(1e-7,final_pp);
td_pb_use = rho.*g.*cosd(theta).*final_h-td_ppval_use;

td_mubf = td_ppval_use.*mu_Iv_fn(td_iv_use);

td_D = -2/beta./final_h.*(td_pb_use-rho_f*g*cosd(theta)*final_h);

td_absFriction = (1./rho).*(-td_mubf - tau0 + (rho-rho_f).*td_D.*final_u);
td_force_bal = final_h*g*sind(theta)+td_absFriction;

%%
ode_name = "time_dep_match.txt";
dirname = "~/Documents/MATLAB/1D_System/xVariation/Nonlinear_Waves_Iv/Rauter_Waves/Results/";
master_file = load(dirname+ode_name);
xi_ode = master_file(1,:);
y_ode = master_file(2:end,:);


record = readtable(dirname+'full_record.csv');
in_table = strcmp(record.Name, ode_name);
ode_wave_type = record.wave_type(in_table);
ode_pres_h = (ode_wave_type == "full_pres_h");
ode_theta = record.theta(in_table); 
ode_lambda = record.lambda(in_table);
ode_Fr = record.Fr(in_table);
ode_nu = record.nu(in_table);
ode_tau0 = record.tau0(in_table);
%         if strcmp(wave_type,'full')
ode_d = record.d(in_table);
ode_a = record.a(in_table);
ode_u_w = record.u_w(in_table);
phi_rlp = 0.53;


[h0_ode, crit_Iv] = crit_Iv_rauter(ode_theta, rho_p, rho_f, eta_f, ode_Fr, ode_a, phi_rlp, phi_rcp, ode_tau0);
u_eq = ode_Fr*sqrt(g*cosd(ode_theta)*h0_ode);
pp_eq_ode = 3*eta_f*u_eq/crit_Iv/h0_ode;
phi_eq_ode = pp_eq_ode/g/cosd(ode_theta)/h0_ode/(rho_p-rho_f);
pb_eq_ode = rho_f*g*cosd(ode_theta)*h0_ode;
ode_u_w = ode_u_w*u_eq;

xi_ode = xi_ode*h0_ode;
Q1_ode = y_ode(1,:)*h0_ode*u_eq;
h_ode = y_ode(2,:)*h0_ode;
u_ode = ode_u_w - Q1_ode./h_ode;
phi_ode = y_ode(5,:)./y_ode(1,:);
%%
hold on
plot(final_grid,final_h);
plot(xi_ode,h_ode);

