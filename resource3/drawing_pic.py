import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def draw_block(r_R,c_R,beta,save_path,epoch,eta,T):
    name = save_path + str(epoch)
    r_R,c_R,beta
    fig = plt.figure(figsize=(10, 7))
    ax  = fig.add_subplot(111)
    ax.plot(r_R,beta,label='beta',color='r')
    ax.legend(loc='left')
    ax2 = ax.twinx()
    ax2.plot(r_R,c_R,label='chord')
    print('max c:',np.max(c_R),'max beta:',np.max(beta))
    ax2.legend(loc='right')

    plt.title("eta = {:.3f},T={:.3f}".format(eta,T))
    plt.savefig(name)
    plt.close()

def draw_process(process_eta,base_eta,process_T,base_T,fig_path):
    fig = plt.figure(figsize=(10, 7))
    ax  = fig.add_subplot(111)
    ax.plot(process_eta,label='eta',color = 'r')
    ax.axhline(y = base_eta,color = 'r')
    ax.legend(loc='left')
    ax2 = ax.twinx()
    ax2.plot(process_T,label = 'T',color = 'b')
    ax2.axhline(y = base_T,color = 'b')
    ax2.legend(loc='right')
    plt.savefig(fig_path+'_process')
    plt.close()

