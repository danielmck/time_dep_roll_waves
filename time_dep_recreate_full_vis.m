function [xi_final,y_final] = time_dep_recreate_full_vis(h0_dim,theta,lambda_dim,tau0,d,alpha,rel_flux,pres_h)
    [phi_c,rho_f,rho_p,rho,eta_f,g] = get_params_water();
    P = (rho-rho_f)/rho;
    
    [Fr_eq,~] = crit_Iv_tau0_h(theta, rho_p, rho_f, eta_f, h0_dim, tau0,0);
    lambda = lambda_dim/h0_dim;
    u_eq_dim = Fr_eq*sqrt(g*cosd(theta)*h0_dim);
    pb_dim = rho_f*g*cosd(theta)*h0_dim;
    stat_dist = 0;
    
    if ~exist("rel_flux","var")
        rel_flux=1;
    end
    
    if ~exist("pres_h","var")
        pres_h=false;
    end
    
%     Res_dir = "~/Documents/MATLAB/1D_System/xVariation/Nonlinear_Waves_Iv/Results/";
    [xi_dl,y_dl] = bvp_full_from_master([Fr_eq,theta,lambda,2.5e-5,alpha,d,tau0,rel_flux,pres_h]);
    
    xi_final = (stat_dist/lambda+xi_dl*(1-stat_dist/lambda))*h0_dim;
    y_final = [u_eq_dim,u_eq_dim*h0_dim,h0_dim,1,u_eq_dim*h0_dim,u_eq_dim*h0_dim,pb_dim]'.*y_dl;
end
    
    