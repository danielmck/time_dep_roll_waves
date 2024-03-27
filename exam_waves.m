% dirname = 'channel_roll_wave_vis/results/tau0_0_theta_12_1500'; %change_t_10_finalTheta_7_initTheta_10_2000
dirname = 'four_eqn_var_rho/results/lambda_200_tau0_0_theta_12_4000'; %'theta_12_1500'; %lambda_15_lambda_300_
dat=hs.Load(dirname);
% hs.Plot(dat,3,[0.0,0.3])

final = dat(100);
theta = final.params.theta;
lambda = final.xSize;
% rho = final.params.rho;
% chi = final.params.chi;
g = final.params.g;
rho_p = final.params.rhog;
rho_f = final.params.rhof;
% eta_f = final.params.etaf;
phi_c = final.params.phim;

init = dat(1);
init_y = permute(init.data,[3,1,2]);
h0 = (min(init_y(1,:))+max(init_y(1,:)))/2;

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
[Fr,eq_Iv] = crit_Iv_tau0_h(theta, rho_p, rho_f, eta_f, h0, tau0,0,true);
u_eq = Fr*sqrt(g*cosd(theta)*h0);
phi_eq = phi_c/(1+sqrt(eq_Iv));
rho_eq = rho_p*phi_eq+rho_f*(1-phi_eq);

final_grid = final.xGrid;
npts = size(final_grid,2);
final_y = permute(final.data,[3,1,2]);
final_h = final_y(1,:);

if ~all(dirname(18:23)=='inflow') && ~all(dirname(19:24)=='inflow')
%     ind = sum(final_grid<12.6);
    [~,ind] = min(final_h);
    final_h = horzcat(final_h(ind:end),final_h(1:ind-1));
    final_grid = mod((horzcat(final_grid(ind:end),final_grid(1:ind-1))-final_grid(ind)),lambda);
    final_y = horzcat(final_y(:,ind:end),final_y(:,1:ind-1));
end

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
    if final.nDims > 3
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
%     ppval = rho.*g.*cosd(theta).*final_h-final_pb;
else
    rho=final.params.rho;
    final_phi = final.params.phim;
    final_pb = final.params.rhof.*g.*cosd(theta).*final_h;
    final_pp = (rho-final.params.rhof).*g.*cosd(theta).*final_h;
    final_pe = zeros(size(final_h));
end
final_pf = rho_f.*g.*cosd(theta).*final_h+final_pe;
final_ptot = rho.*final_h.*g.*cosd(theta);

P = (rho-rho_f)./rho;


% if (absu == 0)
%     signU = 1;
% else 
%     signU = uval/absu;
% end

beta = 150*phi_c^2*eta_f/((1-phi_c)^3*d^2);

zeta = 3./(2.*alpha.*final_h) + rho_f*g*cosd(theta).*P/4;

iv = (3*eta_f).*(abs(final_u)./final_h)./final_pp;

iv_use = max(iv,1e-7);
ppval_use = max(1e-7,final_pp);
pb_use = rho.*g.*cosd(theta).*final_h-ppval_use;

mubf = ppval_use.*mu_Iv_fn(iv_use);

D = -2/beta./final_h.*(pb_use-rho_f*g*cosd(theta)*final_h);

absFriction = (1./rho).*(-mubf - tau0 + (rho-rho_f).*D.*final_u);
force_bal = final_h*g*sind(theta)+absFriction;

if final.nDims > 2
    tanpsi = final_phi - phi_c./(1+sqrt(iv_use));
    dilatancy = 9/2./alpha./final_h.*final_u.*tanpsi;

    psi1 = D.*P;
    psi5 = zeta.*D - dilatancy;
    source_h = psi1;
    source_hphi = -D.*final_phi.*rho_f./rho;

    source_pbh = (psi5-g*cosd(theta)*rho_f/4.*psi1).*final_h+(pb_use-rho*g*cosd(theta).*chi.*final_h).*psi1;
end

ode_denom = final_h.^3*g*cosd(theta)-final_hu.^2;

max_pos = zeros(size(dat,2),1);
tvals = zeros(size(dat,2),1);
lambda = dat(1).xGrid(end);
for i=1:size(dat,2)
    curr = dat(i);
    curr_y = permute(curr.data,[3,1,2]);
    curr_xi = curr.xGrid;
    [h_max, h_loc] = max(curr_y(1,:));
    tvals(i) = curr.time;
    max_pos(i) = curr_xi(h_loc);
end
u_w_time = (mod(diff(max_pos),lambda))./diff(tvals);

if final.nDims > 3
    t_val = sum(final_grid<35)+1;
    sin_ave = 0;
    for j = 2:npts
        sin_ave = sin_ave + (rho(j)*final_h(j)+rho(j-1)*final_h(j-1))/2*(final_grid(j)-final_grid(j-1))/lambda;
    end
end
u_w = u_w_time(end);
% 
Q1 = final_h*u_w-final_hu;

final_m = zeros(1,npts);
int = 0;
for i=2:npts
    int = int+1/2*(rho(i)*final_h(i)+rho(i-1)*final_h(i-1))*(final_grid(i)-final_grid(i-1))/rho_eq/lambda;
    final_m(i) = int;
end
dhdxi = diff(final_h)./diff(final_grid);
dhdxi(npts)=dhdxi(1);
% 
% [h_min, min_ind] = min(final_h);
% hu_min = final_hu(min_ind);
% Q1_min = u_w*h_min-hu_min;
% pb_min = final_pb(min_ind);
% h_max = (-h_min + sqrt(h_min.^2+8*Q1_min^2./(h_min*g*cosd(theta))))/2;
% if final.nDims > 2
%     pb_max = pb_min+rho(end).*g.*cosd(theta).*chi(end).*(h_max-h_min);
% end
% min_len = 10;
% num_unstab = zeros(1,npts);
% for l=1:npts
%     test = Fr_stab_noneq(final_h(l),final_u(l),final_phi(l),final_pb(l),theta, rho_p, rho_f, d, eta_f, alpha,tau0);
%     num_unstab(l) = test(1);
%     min_len = min(min_len,test(2));
% end
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
% theta_change = 8.2;
% pbh_exp = (final_pbh-rho.*g.*cosd(theta_change).*chi.*final_h).*final_h;
in_range = (final_grid>-1); % & final_grid>0.4);
n_wave = 1;
crop_out = [];
for k=1:n_wave
    crop_out = horzcat(crop_out,final_h(in_range));
end
save("four_eqn_var_rho/time_d_load_h.txt","crop_out","-ascii")
crop_out = [];
for k=1:n_wave
    crop_out = horzcat(crop_out,final_hu(in_range));
end
% crop_out = final_hu(in_range);
save("four_eqn_var_rho/time_d_load_hu.txt","crop_out","-ascii")
crop_out = [];
for k=1:n_wave
    crop_out = horzcat(crop_out,final_hphi(in_range));
end
% crop_out = final_hphi(in_range);
save("four_eqn_var_rho/time_d_load_hphi.txt","crop_out","-ascii")
crop_out = [];
for k=1:n_wave
    crop_out = horzcat(crop_out,final_pbh(in_range));
end
% crop_out = final_pbh(in_range);
save("four_eqn_var_rho/time_d_load_pbh.txt","crop_out","-ascii")
%%

hold on
% figure(2)
% SetPaperSize(10,10)

[h_max,max_ind] = max(final_h);
min_ind=1;
% max_ind = size(final_h,2);
eff_fric_coeff = mubf./(rho.*final_h)/g/cosd(theta);

% plot(final_grid(min_ind:end),final_pb(min_ind:end))

iv_list = linspace(1e-6,2e-3,200);
phi_list = phi_c./(1+sqrt(iv_list));
rho_list = rho_p*phi_list+(1-phi_list)*rho_f;
% plot(iv_list,mu_Iv_fn(iv_list).*(rho_list-rho_f)./rho_list,"DisplayName","$\mu(I_v)$ with buoyancy")
% plot(iv(min_ind:max_ind),eff_fric_coeff,"DisplayName","Wave profile")
% dirx = [-4.5,-3,-2,1,1,-2,1]*1e-5;
% diry = [0,-1,-1,0,0,-1,0]*1e-2;
% for k = [1,3,5,6,7]
%     ind = sum(final_grid(1:max_ind)<(k-1)*0.1)+1;
%     scatter(iv(ind),eff_fric_coeff(ind),'kx')
%     txt = ['$\xi=$',num2str((k-1)*0.1),'m'];
%     text(iv(ind)+dirx(k), eff_fric_coeff(ind)+diry(k),txt,'Interpreter','latex')
% end

% plot(final_grid,u_eq*ones(size(final_grid)),"DisplayName","Uniform equilibrium")
subplot(2,2,1)
plot(final_grid,final_h,"DisplayName","Wave Profile")
hold on
xlabel("$\xi$ (m)","interpreter","latex")
ylabel("$h$ (m)","interpreter","latex")
subplot(2,2,2)
plot(final_grid,final_u,"DisplayName","Wave Profile")
hold on
ylabel("$u$","interpreter","latex")
xlabel("$\xi$ (m)","interpreter","latex")
subplot(2,2,3)
plot(final_grid,final_phi,"DisplayName","Wave Profile")
hold on
xlabel("$\xi$ (m)","interpreter","latex")
ylabel("$\phi$","interpreter","latex")
subplot(2,2,4)
plot(final_grid,final_pb,"DisplayName","Wave Profile")
% hold on
ylabel("$p_b$","interpreter","latex")
xlabel("$\xi$ (m)","interpreter","latex")

% plot(final_grid,source_pbh./final_h,"DisplayName","Wave Profile")
% plot(final_grid(end),pb_max,"x","DisplayName","Shock condition maximum")
% plot(final_grid,final.params.rhof.*g.*cosd(theta).*final_h)
% ylabel("$h$ ($m$)")
% ylabel("$\bar{u}$ (ms$^{-1}$)")
% xlabel("$I_v$")
% ylabel("Effective friction coeff.")
% ylabel("$\phi$")
% ylabel("$p_b$ (Pa)")
% ylabel("Liquefaction ratio $\frac{p_f}{p_{tot}}$") %
xlabel("$\xi$ (m)")
% xlim([1.6e-3,1.9e-3])
% ylim([30,205])
% legend("Location","best")
% title("Periodic waveform with $h_0 = 6.1$mm")
% title("$\theta = "+num2str(theta)+"^{\circ}$, $\tau_0 = "+num2str(tau0)+"$Pa, $t="+num2str(final.time)+"$s")
% exp_graph(gcf,"lambda_100_iv_fric_comp.pdf")