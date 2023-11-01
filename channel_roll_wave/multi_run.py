import os
n_run = 19;
os.system('make')
# os.system('pwd')
for x in range(1):
    finalTheta = 8.0 + x*1.5/(n_run-1)
    for y in range(1,n_run):
        lambda_in = 20+(y-1)*180/(n_run-1)
        # print('.\\four_eqn_model '+str(finalTheta)+' '+str(change_t))
        # if not os.path.isdir("results/change_t_"+str(change_t)+"_finalTheta_"+str(finalTheta)+"_initTheta_10_800"):
        # if (x != 2):
        os.system('./channel_roll_wave '+str(lambda_in))
