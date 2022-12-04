from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.problems import get_problem
from pymoo.operators.crossover.pntx import TwoPointCrossover
from pymoo.operators.mutation.bitflip import BitflipMutation
from pymoo.operators.sampling.rnd import BinaryRandomSampling,FloatRandomSampling
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter
from pymoo.core.problem import Problem
import pymoo.gradient.toolbox as anp

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd
import scipy.interpolate as interpolate
from scipy.signal import savgol_filter
from bisect import bisect_left
from BEMT_program.solver import Solver

fig_path = 'process/'

def draw_block(r_R,c_R,beta,save_path,epoch,eta,T):
    name = save_path + str(epoch)
    r_R,c_R,beta
    fig = plt.figure(figsize=(10, 7))
    ax  = fig.add_subplot(111)
    ax.plot(r_R,beta,label='beta',color='r')
    ax.legend(loc='left')
    ax2 = ax.twinx()
    ax2.plot(r_R,c_R,label='c/R')
    print('max c:',np.max(c_R),'max beta:',np.max(beta))
    ax2.legend(loc='right')

    plt.title("eta = {:.3f},T={:.3f}".format(eta,T))
    plt.savefig(name)
    plt.close()

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
 
    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

BEM_path = 'example_OPTIMIZE/prop.ini'
s = Solver(BEM_path)
R = s.rotor.diameter
sections = s.rotor.sections
foil_name = []
radius   = []
pitch    = []
chord    = []

v_inf    = s.v_inf
rpm      = s.rpm
for i in sections:
    foil_name.append(i.airfoil.name)
    radius.append(i.radius)
    pitch.append(i.pitch)
    chord.append(i.chord/R)

num_sec = len(foil_name)
radius,pitch,c_R = np.array(radius),np.array(pitch),np.array(chord)


k = 4
lenth = 13
tail  = np.array([0,0,0,0,0])

t_pitch,c_pitch,_  = interpolate.splrep(radius,pitch,s=0,k=k)
c_pitch = c_pitch[:lenth]
t_cR,c_cR,_  = interpolate.splrep(radius,c_R,s=0,k=k)          
c_cR    = c_cR[:lenth]

para_one = np.max(c_pitch)/np.max(c_cR)
c_cR0 = np.concatenate((c_cR,tail))

# c_cR0 = c_cR + np.random.randn(len(c_cR))*0.01

# print(c_cR,c_cR0)

# exit()


with HiddenPrints():
    df, sections = s.run_sweep('v_inf', 1, v_inf, v_inf)
    base_eta     = df['eta'][0]
    base_T       = df["T"][0]
draw_block(radius,c_R*R,pitch*180/np.pi,'test/','base',base_eta,base_T)

# print(c_cR,c_pitch)
# exit()


process_eta = []
process_T   = []

min_T   = base_T*0.9
min_eta = base_eta*0.9*100

local_cl    = [(0.05,1.2),(0.10,1.4)]
local_pitch = [(0.05,30),(0.10,20)]
local_cR    = [(0.05,0.03),(0.10,0.03)]


degree = (0.8,3)
def f3(t, x):
    i = bisect_left(t, x)
    if t[i] - x > 0.5:
        i-=1
    return i
    
class SphereWithConstraint(Problem):

    def __init__(self, n_var=2*lenth, n_obj=2,n_ieq_constr = 2):
        self.l  = 0
        self.radius = radius
        super().__init__(n_var=n_var, n_obj=n_obj, n_ieq_constr=n_ieq_constr, xl=-0.1, xu=0.1)

    def _evaluate(self, x, out, *args, **kwargs):
        self.l += 1
        radius = self.radius
        print('current loop:',self.l)
        tail0     = np.expand_dims(tail,0).repeat(x.shape[0],axis=0)
        t_pitch0  = np.expand_dims(t_pitch,0).repeat(x.shape[0],axis=0)
        t_cR0     = np.expand_dims(t_cR,0).repeat(x.shape[0],axis=0)
        pitch0   = np.expand_dims(c_pitch,0).repeat(x.shape[0],axis=0)
        c_R0     = np.expand_dims(c_cR,0).repeat(x.shape[0],axis=0)
        
        pitch0  = np.add(x[:,:lenth],pitch0)
        c_R0     = np.add(x[:,lenth:]*para_one,c_R0)
        
        pitch0   = np.concatenate((pitch0,tail0),1)
        c_R0     = np.concatenate((c_R0,tail0),1)
        
        pitch_x = []
        c_Rx    = []
    
        for i in range(len(t_pitch0)):
            t1 = t_pitch0[i]
            t2 = t_cR0[i]
            p1 = pitch0[i]
            p2 = c_R0[i]

            pitch_x.append(interpolate.splev(radius,(t1,p1,k))*180/np.pi)
            c_Rx.append(interpolate.splev(radius,(t2,p2,k))*R)
        pitch_x = np.array(pitch_x)
        c_Rx    = np.array(c_Rx)
        
        # print(pitch0*180/np.pi,pitch_x[0])
        eta = []
        T   = []
        g_cl_list    = []
        g_pitch_list = []
        g_cR_list    = []
        
        for i in range(len(pitch_x)):
        
            c1 = c_Rx[i]
            p1 = pitch_x[i]
            
            window_lenth = int(len(c1)*degree[0])
            if window_lenth%2 ==0:
                window_lenth = window_lenth -1
            smooth_k = int(degree[1])

            c1 = savgol_filter(c1,window_lenth,smooth_k)
            p1 = savgol_filter(p1,window_lenth,smooth_k)

            s1       = Solver(BEM_path,c1=c1,p1=p1)
            with HiddenPrints():
                df, sections = s1.run_sweep('v_inf', 1, v_inf, v_inf)
            if df['eta'][0]>1:
                df['eta'][0] = 100
            else:
                df['eta'][0] = df['eta'][0]*100
            
            
            va     = sections[0].values

            radius = va[:,0]
            chord  = va[:,1]
            pitch  = va[:,2]*180/np.pi
            cl     = va[:,3]

            po_cl_list    = []
            po_pitch_list = []
            po_cR_list    = []
            for _ in local_cl:
                po_cl_list.append(f3(radius,_[0]))
            for _ in local_pitch:
                po_pitch_list.append(f3(radius,_[0]))
            for _ in local_cR:
                po_cR_list.append(f3(radius,_[0]))
            
            a = 0
            for _ in range(len(po_cl_list)):
                a += (cl[po_cl_list[_]] - local_cl[_][1] )**2
            g_cl_list.append(a)

            a = 0
            for _ in range(len(po_pitch_list)):
                a += (pitch[po_pitch_list[_]] - local_pitch[_][1] )**2
            g_pitch_list.append(a)

            a = 0
            for _ in range(len(po_cR_list)):
                a +=(chord[po_cR_list[_]] - local_cR[_][1])**2
            g_cR_list.append(a)
            

            eta.append(df['eta'][0])
            T.append(df['T'][0])
        eta = np.array(eta)
        T   = np.array(T)

        process_eta.append(np.max(eta)/100)
        process_T.append(np.max(T))
        
        #plot figure
        fig = plt.figure(figsize=(10, 7))
        ax  = fig.add_subplot(111)
        ax.plot(process_eta,label='eta',color = 'r')
        ax.axhline(y = base_eta)
        ax.legend(loc='left')
        ax2 = ax.twinx()
        ax2.plot(process_T,label = 'T',color = 'b')
        ax2.axhline(y = base_T)
        ax2.legend(loc='right')
        plt.savefig(fig_path+'process')
        plt.close()

        f1  = -1*eta 
        f2  = -1*T

        g1  = min_T - T
        g2  = min_eta - eta
        
        g1 =np.array(g1)
    
        np.save('test',x)
        out["F"] = np.column_stack([f1, f2])
        a = np.column_stack([f1, f2])
        b = np.column_stack([a,a])
        print(b.shape)
        out["G"] = np.column_stack([g1, g2])
        # out["G"] = 0.1 - out["F"]
        exit()
 

# x为基因型，F为表现型最小化函数，G为不等式约束<=0


problem =SphereWithConstraint()
algorithm = NSGA2(pop_size=5,
                  sampling=FloatRandomSampling(),
                #   crossover=TwoPointCrossover(),
                #   mutation=BitflipMutation(),
                  eliminate_duplicates=True)
res = minimize(problem,
               algorithm,
               ('n_gen', 10),
               seed=1,
            #    save_history = True,
            #    verbose=True)
)


def get_BEM(r_R,pitch,c_R):

    with HiddenPrints():
    
        s = Solver(BEM_path,r1=r_R,p1=pitch,c1=c_R)
        df, sections = s.run_sweep('v_inf', 1, v_inf, v_inf)
        base_eta     = df['eta'][0]
        base_T       = df["T"][0]
    
    return base_eta,base_T


X = np.load('test.npy')
for i in range(len(X)):
    pitch0  = np.add(X[i][:lenth],c_pitch)
    c_R0     = np.add(X[i][lenth:],c_cR)
    
    pitch   = np.concatenate((pitch0,tail))
    c_R     = np.concatenate((c_R0,tail))
    
    pitch_x  = interpolate.splev(radius,(t_pitch,pitch,k))*180/np.pi
    c_Rx     = interpolate.splev(radius,(t_cR,c_R,k))*R

    window_lenth = int(len(c_Rx)*degree[0])
    if window_lenth%2 ==0:
        window_lenth = window_lenth -1
    smooth_k = int(degree[1])
    c_Rx = savgol_filter(c_Rx,window_lenth,smooth_k)
    pitch_x = savgol_filter(pitch_x,window_lenth,smooth_k)
    eta,T = get_BEM(radius,pitch_x,c_Rx)
    draw_block(radius,c_Rx,pitch_x,'test/',i,eta,T)

