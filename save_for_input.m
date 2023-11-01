function save_for_input(dirname)
    dat=hs.Load(dirname);
    final = dat(end);
    final_y = permute(final.data,[3,1,2]);
    final_h = final_y(1,:);
    final_hu = final_y(2,:);
    final_hphi = final_y(3,:);
    final_pbh = final_y(4,:);

    save('four_eqn_var_rho/time_d_load_h.txt','final_h',"-ascii")
    save("four_eqn_var_rho/time_d_load_hu.txt","final_hu","-ascii")
    save("four_eqn_var_rho/time_d_load_hphi.txt","final_hphi","-ascii")
    save("four_eqn_var_rho/time_d_load_pbh.txt","final_pbh","-ascii")
end