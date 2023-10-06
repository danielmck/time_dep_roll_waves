import os
n_run = 6;
os.system('make')
# os.system('pwd')
for x in range(1):
    finalTheta = 8.0 + x*1.5/(n_run-1)
    for y in range(1,n_run):
        change_t = 5.0 + y*25/(n_run-1)
        # print('.\\four_eqn_model '+str(finalTheta)+' '+str(change_t))
        # if not os.path.isdir("results/change_t_"+str(change_t)+"_finalTheta_"+str(finalTheta)+"_initTheta_10_800"):
        # if (x != 2):
        os.system('./four_eqn_model '+str(finalTheta)+' '+str(change_t))
