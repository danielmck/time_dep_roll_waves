import os
import time
# time.sleep(21600)
n_run = 19
os.system('make')
# os.system('pwd')
for x in range(1):
    finalTheta = 8.0 + x*1.5/(n_run-1)
    for y in range(n_run):
        lambda_in = 20+y*180/(n_run-1)
        # print('.\\four_eqn_model '+str(finalTheta)+' '+str(change_t))
        # if not os.path.isdir("results/change_t_"+str(change_t)+"_finalTheta_"+str(finalTheta)+"_initTheta_10_800"):
        # if (x != 2):
        os.system("matlab -nodisplay -nosplash -nodesktop -r \"addpath ../../matlab;save_for_input('results/lambda_"+str(int(lambda_in))+"_tau0_0_theta_12_4000');exit;\"")
        os.system('./channel_roll_wave '+str(lambda_in))
