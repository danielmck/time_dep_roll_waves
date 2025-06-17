% Debugging code recreating the terms that have been used in the model
dirname = 'four_eqn_var_rho/results/lambda_300_theta_12_3000'; 
dat=hs.Load(dirname);
% hs.Plot(dat,1,[0.0,3.0])50

first = dat(1);
second = dat(2);
theta = first.params.theta;
lambda = first.xSize;
% rho = first.params.rho;
% chi = first.params.chi;
g = first.params.g;
rho_p = first.params.rhog;
rho_f = first.params.rhof;
% eta_f = first.params.etaf;
phi_c = first.params.phim;

if isfield(first.params,"d")
    d=first.params.d;
else
    d=1e-4;
end

if isfield(first.params,"alpha")
    alpha=first.params.alpha;
else
    alpha=1e-5;
end
if isfield(first.params,"etaf")
    eta_f=first.params.etaf;
else
    eta_f=1.0013e-3;
end

if isfield(first.params,"tau0")
    tau0=first.params.tau0;
else
    tau0=0;
end

first_grid = first.xGrid;
npts = size(first_grid,2);
first_y = permute(first.data,[3,1,2]);
first_h = first_y(1,:);
first_hu = first_y(2,:);
first_u = first_hu./first_h;
first_hu2 = first_hu.*first_u;
first_hphi = first_y(3,:);
 
first_phi = first_hphi./first_h;

first_hphiu = first_hphi.*first_u;

first_rho=first.params.rhog*first_phi+first.params.rhof*(1-first_phi);
first_chi=(first.params.rhof+3*first_rho)./first_rho/4;
first_pbh = first_y(4,:);
first_pb = first_pbh./first_h+first_rho.*g.*cosd(theta).*first_chi.*first_h;
first_pp = first_h.*first_rho.*g.*cosd(theta)-first_pb;
first_pe = first_pb - first.params.rhof.*g.*cosd(theta).*first_h;

second_grid = second.xGrid;
second_y = permute(second.data,[3,1,2]);
second_h = second_y(1,:);

second_hu = second_y(2,:);
second_u = second_hu./second_h;
second_hu2 = second_hu.*second_u;

second_hphi = second_y(3,:);
second_hphiu = second_hphi.*second_u;

second_phi = second_hphi./second_h;
second_rho=first.params.rhog*second_phi+first.params.rhof*(1-second_phi);
second_chi=(first.params.rhof+3*second_rho)./second_rho/4;

second_pbh = second_y(4,:);
second_pb = second_pbh./second_h+second_rho.*g.*cosd(theta).*second_chi.*second_h;
second_pp = second_h.*second_rho.*g.*cosd(theta)-second_pb;
second_pe = second_pb - second.params.rhof.*g.*cosd(theta).*second_h;

dhdt = (second_h-first_h)/0.05;
dhudt = (second_hu-first_hu)/0.05;
dhphidt = (second_hphi-first_hphi)/0.05;
dpbhdt = (second_pbh-first_pbh)/0.05;

first_dhdx = [diff(first_h) first_h(end)-first_h(1)]/(first_grid(2)-first_grid(1));
second_dhdx = [diff(second_h) second_h(end)-second_h(1)]/(second_grid(2)-second_grid(1));
first_dhudx = [diff(first_hu) first_hu(end)-first_hu(1)]/(first_grid(2)-first_grid(1));
second_dhudx = [diff(second_hu) second_hu(end)-second_hu(1)]/(second_grid(2)-second_grid(1));
first_dhphiudx = [diff(first_hphiu) first_hphiu(end)-first_hphiu(1)]/(first_grid(2)-first_grid(1));
second_dhphiudx = [diff(second_hphiu) second_hphiu(end)-second_hphiu(1)]/(second_grid(2)-second_grid(1));
first_dhu2dx = [diff(first_hu2) first_h(end)-first_h(1)]/(first_grid(2)-first_grid(1));
second_dhu2dx = [diff(second_hu2) second_h(end)-second_h(1)]/(second_grid(2)-second_grid(1));

first_P = (first_rho-rho_f)./first_rho;
second_P = (second_rho-rho_f)./second_rho;

beta = 150*phi_c^2*eta_f/((1-phi_c)^3*d^2);

first_zeta = 3./(2.*alpha.*first_h) + rho_f*g*cosd(theta).*first_P/4;
second_zeta = 3./(2.*alpha.*second_h) + rho_f*g*cosd(theta).*second_P/4;

first_iv = (3*eta_f).*(abs(first_u)./first_h)./first_pp;
first_iv = max(first_iv,1e-8);
second_iv = (3*eta_f).*(abs(second_u)./second_h)./second_pp;
second_iv = max(second_iv,1e-8);

first_mubf = first_pp.*mu_Iv_fn(first_iv);
second_mubf = second_pp.*mu_Iv_fn(second_iv);

first_D = -2/beta./first_h.*(first_pb-rho_f*g*cosd(theta)*first_h);
second_D = -2/beta./first_h.*(second_pb-rho_f*g*cosd(theta)*second_h);

first_absFriction = (1./first_rho).*(-first_mubf - tau0 + (first_rho-rho_f).*first_D.*first_u);
first_force_bal = first_h*g*sind(theta)+first_absFriction;

second_absFriction = (1./second_rho).*(-second_mubf - tau0 + (second_rho-rho_f).*second_D.*second_u);
second_force_bal = first_h*g*sind(theta)+second_absFriction;

first_tanpsi = first_phi - phi_c./(1+sqrt(first_iv));
first_dilatancy = 9/2./alpha./first_h.*first_u.*first_tanpsi;

second_tanpsi = second_phi - phi_c./(1+sqrt(second_iv));
second_dilatancy = 9/2./alpha./second_h.*second_u.*second_tanpsi;

first_psi1 = first_D.*first_P;
first_psi5 = first_zeta.*first_D - first_dilatancy;

second_psi1 = second_D.*second_P;
second_psi5 = second_zeta.*second_D - second_dilatancy;

first_source_h = first_psi1;
first_source_hphi = -first_D.*first_phi.*rho_f./first_rho;
first_source_pbh = (first_psi5-g*cosd(theta)*rho_f/4.*first_psi1).*first_h+(first_pb-first_rho*g*cosd(theta).*first_chi.*first_h).*first_psi1;

second_source_h = second_psi1;
second_source_hphi = -second_D.*second_phi.*rho_f./second_rho;
second_source_pbh = (second_psi5-g*cosd(theta)*rho_f/4.*second_psi1).*second_h+(second_pb-second_rho*g*cosd(theta).*second_chi.*second_h).*second_psi1;

first_dhdt_eqn = -first_dhudx+first_P.*first_D;
second_dhdt_eqn = -second_dhudx+second_P.*second_D;

first_dhudt_eqn = -first_dhu2dx-g*cos(theta).*first_h.*first_dhdx+first_force_bal;
second_dhudt_eqn = -first_dhu2dx-g*cos(theta).*second_h.*second_dhdx+second_force_bal;

hold on
% plot(first_grid,dhudt)
plot(first_grid,first_h)
% plot(first_grid,second_dhudt_eqn)