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

from BEMT_program.solver import Solver


def draw_block(r_R,c_R,beta,save_path,epoch):
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
draw_block(radius,c_R*R,pitch*180/np.pi,'test/','base')
with HiddenPrints():
    df, sections = s.run_sweep('v_inf', 1, v_inf, v_inf)
    base_eta     = df['eta'][0]
    base_T       = df["T"][0]
# print(df['eta'][0])
# print(np.add(pitch,c_R))
# print(c_R*R,pitch)
# print(pitch,c_R)
# s1 = Solver(BEM_path,c1=c_R*R,p1=pitch*180/np.pi)
# df, sections = s1.run_sweep('v_inf', 1, v_inf, v_inf)



class SphereWithConstraint(Problem):

    def __init__(self, n_var=2*num_sec, n_obj=2,n_ieq_constr = 0):
        self.l  = 0
        self.l0 = 0
        super().__init__(n_var=n_var, n_obj=n_obj, n_ieq_constr=n_ieq_constr, xl=-0.5, xu=0.5)

    def _evaluate(self, x, out, *args, **kwargs):
        self.l += 1
      
        print('current loop:',self.l)
        pitch0   = np.expand_dims(pitch,0).repeat(x.shape[0],axis=0)
        c_R0     = np.expand_dims(c_R,0).repeat(x.shape[0],axis=0)
        
        pitch_x  = np.add(x[:,:num_sec],pitch0)*180/np.pi
        c_x      = np.add(x[:,:num_sec],c_R0)*R
        # print(pitch0*180/np.pi,pitch_x[0])
        eta = []
        T   = []
        for i in range(len(pitch_x)):
        
            c1 = c_x[i]
            p1 = pitch_x[i]
            
            s1       = Solver(BEM_path,c1=c1,p1=p1)
            with HiddenPrints():
                df, sections = s1.run_sweep('v_inf', 1, v_inf, v_inf)
            eta.append(df['eta'][0])
            T.append(df['T'][0])
        eta = np.array(eta)
        T   = np.array(T)

        f1  = -1*eta 
        f2  = -1*T
        np.save('test',x)
        out["F"] = np.column_stack([f1, f2])
        # out["G"] = 0.1 - out["F"]

# x为基因型，F为表现型最小化函数，G为不等式约束<=0

problem =SphereWithConstraint()
algorithm = NSGA2(pop_size=50,
                  sampling=FloatRandomSampling(),
                #   crossover=TwoPointCrossover(),
                #   mutation=BitflipMutation(),
                  eliminate_duplicates=True)
res = minimize(problem,
               algorithm,
               ('n_gen', 100),
               seed=1,
               verbose=False)
X,F = res.opt.get("X","F")
np.save(X,'test')

X = np.load('test.npy')
print(X.shape)
for i in range(len(X)):
    pitch0 = X[i][:num_sec]
    c0     = X[i][num_sec:]
    pitchx    = np.add(pitch,pitch0)*180/np.pi
    cx        = np.add(c0,c_R)*R
    draw_block(radius,cx,pitchx,'test/',i)



