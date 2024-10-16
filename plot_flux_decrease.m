function plot_flux_decrease
    [phi_c,rho_f,rho_p,~,eta_f,g] = get_params_water();
    SetPaperSize(15,10)
    
    flux_vals = load("four_eqn_var_rho/flux_decrease_5deg_LR_ave_fine.txt");
%     flux_vals(1,1) = NaN;
    flux_init = flux_vals(:,1);
%     flux_vals = max(flux_vals,0);
%     flux_vals = flux_vals./flux_init;
%     flux_vals(flux_vals==0)=NaN;
%     flux_vals = flux_vals*1000;
    theta=5;
    h0 = 0.0061;
    [Fr, crit_Iv] = crit_Iv_tau0_h(12, 2500, 1000, 0.0010016, h0, 0.0, false, false);
    u_eq = Fr*sqrt(cosd(12)*g*h0);
    flux_eq = u_eq*h0;
%     flux_vals = flux_vals/flux_eq;
    
    t_vals = linspace(0,10,1001);
    lambda_vals = linspace(20,150,14);
    
    da_vals = load('four_eqn_var_rho/results/Ive_da_5deg_12init_change_0.txt');
    da_time = da_vals(:,1);
    da_h = da_vals(:,2);
    da_u = da_vals(:,4);
    da_phi = da_vals(:,3);
    da_pb = da_vals(:,5);
    da_pe = da_pb-rho_f*g*cosd(theta)*da_h;
    da_rho = rho_p*da_phi+(1-da_phi)*rho_f;
    da_LR = da_pe./(da_rho.*g*cosd(theta).*da_h);
    
    
    depo_ind = sum(da_u>1e-5);
    depo_time = da_time(depo_ind);
    hold on
    C = viridis(4);
    plot(t_vals(1:size(flux_vals,2)),flux_vals(4,:), "DisplayName", "$\lambda/h_0 = 50$","color",C(1,:))
    plot(t_vals(1:size(flux_vals,2)),flux_vals(9,:), "DisplayName", "$\lambda/h_0 = 100$","color",C(2,:))
    plot(t_vals(1:size(flux_vals,2)),flux_vals(14,:), "DisplayName", "$\lambda/h_0 = 150$","color",C(3,:))
    plot(da_time, da_LR, 'r--',"DisplayName", "Uniform flow")
%     contourf(t_vals(1:size(flux_vals,2)),lambda_vals,flux_vals,12,'LineStyle','none');
    
%     contour(t_vals(1:size(flux_vals,2)),lambda_vals,flux_vals,[1 1]*h0*1000,'k','linewidth',2);
%     text(1.1,90,{'Contour',' at $h_0$'},'Color','k','Interpreter','latex')
    
%     hcb = colorbar();
%     hcb.Label.String = "Basal excess fluid pressure (Pa)";
%     hcb.Label.String = "Average flux (m$^2$ s$^{-1}$)";
%     hcb.Label.String = "Flux relative to initial value";
%     hcb.Label.String = "Flux relative to uniform flow at $\theta_i$";
%     hcb.Label.String = "Deposited fraction of wavelength";
%     hcb.Label.String = "Lowest flowing point (mm)";
%     ylabel('Average flux (m$^2$ s$^{-1}$)')
%     ylabel('Average basal excess fluid pressure (Pa)')
    ylabel('$\frac{\lambda}{h_0}$');
    xlabel('$t$ (s)')
    title('Full model, $h_0 = 6.1$mm, $\theta_i = 12^{\circ}$, $\theta_f = 5^{\circ}$')
    
%     X = [0.11, 0.445];
%     Y = [0.06, 0.06];
%     annotation('doublearrow',X,Y);
%     X = [0.12, 0.775];
%     Y = [0.06, 0.06];
%     annotation('doublearrow',X,Y);
% %     annotation('textbox',[.2 .2 .4 .3],'String','Linear decreasing $\theta$','EdgeColor','none')
% %     dim = [0.15 0.0 0.01 0.1];
%     annotation('textbox',[0.5 0.015, 0.3 0.044979],'String','Linear decreasing $\theta$','FontSize',9,'EdgeColor','none','Interpreter','latex');
% %     annotation('textbox',[0.58661 0.015, 0.1518 0.044979],'String','$\theta = 5^{\circ}$','FontSize',9,'EdgeColor','none','Interpreter','latex');
% %     annotation('textbox',[0.25 0.085 0.13 0.02910416],'String','$\theta = \theta_c$','FontSize',9,'EdgeColor','none','Interpreter','latex');
%     plot(20*(12-8.506)/7*ones(19,1),lambda_vals,'w--','HandleVisibility','off')
%     plot(depo_time*ones(14,1),lambda_vals,'w--','HandleVisibility','off')
%     text(depo_time-0.65,30,{'Uniform flow', 'deposits'},'Color','w','Interpreter','latex')
%     text(20*(12-8.506)/7+0.2,25,'$\theta=\theta_c$','Color','w','Interpreter','latex')
    xlim([0,3.0]);
%     caxis([0.0, 1.0])
%     ylim([0,1e-3]);
%     ax = gca;
%     ax.YAxis.Exponent = 0;
%     legend()
     %
%     exp_graph(gcf,"Ive_flux_decrease_5deg_dep_frac_fine.pdf")
    
end