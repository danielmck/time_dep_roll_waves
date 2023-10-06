function comp_wave_unif
    dirname = 'four_eqn_var_rho/results/change_t_30_finalTheta_7.6_initTheta_10_4000'; %';finalTheta_8_initTheta_10_tau0_0_2100 tau0_0_theta_8_6000
    dat=hs.Load(dirname);
    [phi_c,rho_f,rho_p,~,eta_f,g] = get_params_water();
    change_t = 30;
    finalTheta = 7.6;
    initTheta = 10;
    
    n_times = size(dat,2);
    flux_vals = zeros(1,n_times);
    u_vals = zeros(1,n_times);
    wave_t = zeros(1,n_times);
    for i = 1:n_times
        curr = dat(i);
        ave_u = 0;
        ave_flux = 0;
        xGrid = curr.xGrid;
        data = permute(curr.data,[3,1,2]);
        lambda = curr.xSize;
        wave_h = data(1,:);
        wave_u = data(2,:)./wave_h;
        wave_phi = data(3,:)./wave_h;
        wave_rho = rho_p.*wave_phi+rho_f.*(1-wave_phi);
        for j=2:curr.xPoints
            ave_u = ave_u + (wave_u(j)+wave_u(j-1))*(xGrid(j)-xGrid(j-1))/2;
            ave_flux = ave_flux + (wave_rho(j).*wave_u(j).*wave_h(j)+wave_rho(j-1).*wave_u(j-1).*wave_h(j-1))/2*(xGrid(j)-xGrid(j-1));
        end
        ave_u = ave_u/lambda;
        ave_flux = ave_flux/lambda;
        wave_t(i) = curr.time;
        u_vals(i) = ave_u;
        flux_vals(i) = ave_flux;
    end
    theta_vec = initTheta + (finalTheta-initTheta)*wave_t/change_t;
    theta_vec(wave_t>change_t) = finalTheta;
    da_root = "~/Documents/MATLAB/1D_System/Iverson_DA/Results/";
    da_file = da_root+"Ive_da_7.6deg_10init_change_30_dim.txt";
    da_data = load(da_file);
    t_uni = da_data(:,1);
    h_uni = da_data(:,2);
    phi_uni = da_data(:,3);
    u_uni = da_data(:,4);
    pb_uni = da_data(:,5);
    rho_uni = rho_p.*phi_uni+rho_f.*(1-phi_uni);
    flux_uni = rho_uni.*h_uni.*u_uni;
    
    hold on
    SetPaperSize(10,10)
    yyaxis left
    plot(wave_t,flux_vals,"DisplayName","Wave average")
    plot(t_uni,flux_uni,"DisplayName","Uniform flow")
    ylim([0,250])
    
    xlim([0,200])
    yyaxis right
    plot(wave_t,theta_vec,"DisplayName","Slope Angle")
    ylabel("Slope angle ($^\circ$)")
    xlabel("Time (s)")
    yyaxis left
    ylabel("Flux (kg m$^2$ s$^{-1}$)")
    legend();
    exp_graph(gcf,"slope_dec_wave_unif_comp.pdf")
end