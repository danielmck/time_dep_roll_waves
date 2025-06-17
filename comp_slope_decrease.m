% Code for comparing the periodic box case when the slope angle is varied in time. The angle is decreased from an initial
% value theta to finalTheta over a time of change_t
n_run = 11;
slows = zeros(n_run,n_run);
for i=1:n_run
    finalTheta = 7+1.5/(n_run-1)*(i-1);
    for j=1:n_run
        change_t = 5+25/(n_run-1)*(j-1);
        dirname = ['four_eqn_var_rho/results/change_t_',num2str(change_t),'_finalTheta_',num2str(finalTheta),'_initTheta_10_800'];
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

        final_grid = final.xGrid;
        final_y = permute(final.data,[3,1,2]);
        final_h = final_y(1,:);

        if ~all(dirname(18:23)=='inflow')
            [~,ind] = min(final_h);
            final_h = horzcat(final_h(ind:end),final_h(1:ind-1));
            final_grid = mod((horzcat(final_grid(ind:end),final_grid(1:ind-1))-final_grid(ind)),lambda);
            final_y = horzcat(final_y(:,ind:end),final_y(:,1:ind-1));
        end

        final_hu = final_y(2,:);
        final_u = final_hu./final_h;

        initial = dat(1);
        init_grid = initial.xGrid;
        init_y = permute(initial.data,[3,1,2]);
        init_h = init_y(1,:);

        if ~all(dirname(18:23)=='inflow')
            [~,ind] = min(final_h);
            init_h = horzcat(init_h(ind:end),init_h(1:ind-1));
            init_grid = mod((horzcat(init_grid(ind:end),init_grid(1:ind-1))-init_grid(ind)),lambda);
            init_y = horzcat(init_y(:,ind:end),init_y(:,1:ind-1));
        end
        init_hu = init_y(2,:);
        init_u = init_hu./init_h;

        init_u_ave = 0;
        final_u_ave = 0;
        for k=2:size(final_grid,2)
            init_u_ave = init_u_ave+(init_u(k)+init_u(k-1))/2*(init_grid(k)-init_grid(k-1))/lambda;
            final_u_ave = final_u_ave+(final_u(k)+final_u(k-1))/2*(final_grid(k)-final_grid(k-1))/lambda;
        end
    end
    slows(i,j) = (final_u_ave>init_u_ave);
end
%%
SetPaperSize(12,10)
% yyaxis left
contourf(5+25/(n_run-1)*(1:n_run),7+1.5/(n_run-1)*(1:n_run),slows,1)
xlabel("Length of slope decrease period (s)")
ylabel("Final slope angle ($^\circ$)", 'color','k')
c = colorbar('Ticks',[0,0.5],'TickLabels',{'Decreasing','Increasing'});
c.Label.String = "Change in u";
c.Label.Position(1) = c.Label.Position(1)-1.5;
exp_graph(gcf,"slope_angle_decrease.pdf")
