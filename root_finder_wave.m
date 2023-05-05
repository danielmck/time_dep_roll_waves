function root_finder_wave

    master_name = "non_vis_convert.txt";
    master_file = load("../Results/"+master_name);
    master_xi = master_file(1,:);
    master_y = master_file(2:end,:);
    record = readtable('../Results/wave_record.csv');

    in_table = strcmp(record.Name, master_name);
    wave_type = record.wave_type(in_table);
    theta = record.theta(in_table);
    lambda_in = record.lambda(in_table);
    Fr = record.Fr(in_table);
    tau0 = record.tau0(in_table);
    alpha = record.alpha(in_table);
    d = record.d(in_table);
    
    
    rho = rho_p*phi_c+rho_f*(1-phi_c);
    P = (rho-rho_f)/rho;
    chi = (rho_f+3*rho)/(4*rho);
    
    [h0, crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0);
    u_eq = Fr*sqrt(g*cosd(theta)*h0);
    phi_eq = phi_c/(1+sqrt(crit_Iv));
    
    crit_pb = rho_f*g*cosd(theta)*h0;
    
    p_tot = rho*g*cosd(theta);
    
    z_scale = h0;
    v_scale = u_eq;
    p_scale = crit_pb;
    t_scale = z_scale/v_scale;

    tau0_dl = tau0/p_scale;
    eta_f_dl = eta_f/(p_scale*t_scale);
    alpha_dl = alpha*p_scale;
    g_dl = g*t_scale/v_scale; 
    
    u_eq_dl = u_eq/v_scale;
    
    rho_f_dl = rho_f*v_scale^2/p_scale;
    rho_p_dl = rho_p*v_scale^2/p_scale;
    d_dl = d/z_scale;
    
    p_tot_grad_dl = p_tot/p_scale*z_scale;
    rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
    
    pres_h = 0;
    h_b_pert = 1e-3;
    h_min = 0.8;
    
    beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);
    
    
    
    function [Q_diff,phi_diff,pb_diff] = project_from_crit(h_c,phi_c,pb_c)
        D_c = -2/beta_dl/h_c*(pb_c-h_c);
        pp_c = p_tot_grad_dl*h_c-pb_c;
        crit_mu_val = (tand(theta)-tau0_dl*rho_f/rho/h+P*D_c*h_c)*(p_tot_grad_dl*h_c)/pp_c;
        Iv_c = 5e-3;
        while (abs(crit_mu_val-mu_Iv_fn(Iv_c))>1e-8)
            Iv_c = Iv_c - crit_mu_val-mu_Iv_fn(Iv_c)/dmudIv_fn(Iv_c);
        end
        u_c = Iv_c/3/eta_f_dl*(h*pp_c);
        u_w = u_c+h_c^(1/2)/Fr;
        Q1_c = h^(3/2)/Fr;
        m_c = 0.5;
        h_max = (-h_min + sqrt(h_min.^2+8*u_w^2./(h_min/Fr^2)))/2;
        opts_l = odeset('Events',@h_stop_l);
        opts_r = odeset('Events',@h_stop_r);
        
        [t_r,y_r] = ode45(@viscous_syst,[0,500],[Q1_c,h_c,m_c,phi_c*Q1_c,pb_c+g_dl*cosd(theta)*rho_dl*chi*h_c]',opts_r);
        
        [t_l,y_l] = ode45(@viscous_syst,[0,-500],[Q1_c,h_c,m_c,phi_c*Q1_c,pb_c+g_dl*cosd(theta)*rho_dl*chi*h_c]',opts_l);
        
        Q_diff = y_r(1,end)-y_l(1,1);
        m_diff = y_r(3,end)-y_l(3,1)-rf;
        phi_diff = y_r(4,end)-y_l(4,1);
        pb_diff = y_r(5,end)-y_l(5,1);
%         total_res = abs(Q_diff)+abs(m_diff)+abs(phi_diff)+abs(pb_diff);
        
        function [value,isterminal,direction] = h_stop_r(t,y)
            value = y(2)-h_max;     % Detect height = 0
            isterminal = 1;   % Stop the integration
            direction = 0;
        end
        
        function [value,isterminal,direction] = h_stop_l(t,y)
            value = y(2)-h_min;     % Detect height = 0
            isterminal = 1;   % Stop the integration
            direction = 0;
        end
    
        function dydxi = viscous_syst(xi,y)
            % The system from the viroulet paper with a different friction
            % law. Need to solve for u_w and Q1 too.
            Q1 = y(1);
            h = y(2);
            u = (-Q1 + h.*u_w)./h;
            m = y(3);
            phi = y(4)./Q1;
            pb = y(5) + g_dl*cosd(theta)*rho_dl*chi.*h;


            zeta = 3/(2*alpha_dl*h) + P/4;
            p_p = p_tot_grad_dl*h-pb;
            Iv = 3*eta_f_dl*abs(u)/h/p_p;
            D = -2/beta_dl/h*(pb-h);

            denom = (h.^3/Fr^2-Q1.^2);
            fric_coeff = p_p/(p_tot_grad_dl*h)*mu_Iv_fn(Iv);
            fb_val = tand(theta)-fric_coeff-tau0_dl*rho_f/rho/h;
            dhdxi = 1/Fr^2.*h.^3.*(fb_val+P*D*h)./denom;

            R_w3 = -phi*rho_f_dl/rho_dl*D;
            R_w4 = (-P*chi+zeta)*D - 2*3/alpha_dl/h*u*(phi - phi_c./(1+sqrt(Iv)));

            dQdxi = -P*D;
            dmdxi = h/(lambda)*u^(1-pres_h);
            dy6dxi = -R_w3;
            dy7dxi = R_w4/(u-u_w);

            dydxi = [dQdxi,dhdxi,dmdxi,dy6dxi,dy7dxi];
            if abs(denom)<denom_eps
                A_pp_coeff = -mu_Iv_fn(Iv)./(p_tot_grad_dl.*h)+p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*Iv./p_p;
                A_pb_coeff = -P.*2./beta_dl-A_pp_coeff;
                A_h_coeff = p_p./(p_tot_grad_dl.*h.^2).*mu_Iv_fn(Iv)+p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*3*eta_f_dl.*abs(u)./h.^2./p_p+tau0_dl.*rho_f./rho./h.^2-P.*2./beta_dl+A_pp_coeff.*rho_dl.*g_dl.*cosd(theta);

                full_h_coeff = (h.^3./Fr^2.*(g_dl*cosd(theta)*rho_dl.*chi.*A_pb_coeff+A_h_coeff-p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*Iv./u.*Q1./h.^2)+3.*h.^2./Fr^2.*(fb_val+P.*D.*h));
                numer_h_term = full_h_coeff.*dhdxi;
                numer_other_term = h.^3./Fr^2.*(dy7dxi.*A_pb_coeff+p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*Iv./u.*dQdxi./h.^2);

                h_quad_roots = roots([3*h.^2./Fr^2, -full_h_coeff-2.*Q1.*dQdxi, -numer_other_term]);
                if imag(h_quad_roots(1))==0
                    dydxi(2) = min(h_quad_roots(h_quad_roots>0));
                else
                    dydxi(2) = 0;
                end
            end
        end
    end
end