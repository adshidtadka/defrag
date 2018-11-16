# -*- coding: utf-8 -*-

from __future__ import print_function
from math import fabs
import cplex
import sys
import math
import datetime

f = open('./../result/ssrr_lno_output.txt')

#number of signals
line = f.readline()
line = line.split("=")
S    = int(line[1])

#number of requests
line = f.readline()
line = line.split("=")
absP = int(line[1])

#number of slots
line = f.readline()
line = line.split("=")
absF = int(line[1])

#number of links
line = f.readline()
line = line.split("=")
absE = int(line[1])

#number of nodes
line = f.readline()
line = line.split("=")
absV = int(line[1])

#number of allowed transitions
line = f.readline()
line = line.split("=")
T = int(line[1])
# print("T = ", T)
# print("absV = ", absV)

line = f.readline()#空行
line = f.readline()#param : s_p k_p_init n_p f_p_init := 

s_p = []
k_p_init = []
n_p = []
f_p_init = []

line = f.readline()
while line[0] != ";":
  line = line.split()
  s_p.append(int(line[0]))      #signal index
  k_p_init.append(int(line[1])) #prim or back
  n_p.append(int(line[2]))      #lp_size
  f_p_init.append(int(line[3])) #spec_ind or bp_ind
  line = f.readline()

# print("s_p = ", s_p)
# print("k_p_init = ", k_p_init)
# print("n_p = ", n_p)
# print("f_p_init = ", f_p_init)

line = f.readline()#空行
line = f.readline()#E := 

E = []
line = f.readline()
while line[0] != ";":
    line = line.split()
    E.append([int(line[0]), int(line[1])])
    line = f.readline()

# print("E = ", E)

line = f.readline()#空行
line = f.readline()#r_p_i_j_init := 

#initial state of link (i,j) of path p
r_p_i_j_init = []
for a in range(absP):
    r_p_i_j_init.append([])
    for b in range(absV):
        r_p_i_j_init[a].append([])
        for g in range(absV):
            r_p_i_j_init[a][b].append(0)

line = f.readline()
while line[0] != ";":
  line = line.split()
  r_p_i_j_init[int(line[0])][int(line[1])][int(line[2])] = 1
  line = f.readline()
# print("r_p_i_j_init = ", r_p_i_j_init)

line = f.readline()#空行
line = f.readline()#o_p d_p := 
o_p = []
d_p = []

line = f.readline()
while line[0] != ";":
    line = line.split()
    o_p.append(int(line[0]))
    o_p.append(int(line[0]))
    d_p.append(int(line[1]))
    d_p.append(int(line[1]))
    line = f.readline()
# print("o_p = ", o_p)
# print("d_p = ", d_p)

f.close()

# epsilon

epsilon = 0.01

#problem variables

x_p_f_t = []
y_p_f_t = []
z_p_f_t_i_j = []
w_p_f_t_i_j = []
r_p_t_i_j = []
m_p_t = []
U_t = []
h_p_t = []
k_p_t = []
l_p_t = []






################################################################





def setupproblem(c):

    c.objective.set_sense(c.objective.sense.minimize)

    # ログを表示しない
    c.set_log_stream(None)
    c.set_error_stream(None)
    c.set_warning_stream(None)
    c.set_results_stream(None)

    # 時間制限
    #c.parameters.timelimit.set(5)

    # assignment variables: x_p_f_t
    allx_p_f_t = []
    for a in range(absP):
        x_p_f_t.append([])
        for b in range(absF):
            x_p_f_t[a].append([])
            for g in range(T+1): #from 0 to T
                varname = "x_" + str(a) + "_" + str(b) + "_" + str(g)
                allx_p_f_t.append(varname)
                x_p_f_t[a][b].append(varname)
    c.variables.add(names=allx_p_f_t, lb=[0] * len(allx_p_f_t),
                    ub=[1] * len(allx_p_f_t),
                    types=["B"] * len(allx_p_f_t))


    # assignment variables: y_p_f_t

    ally_p_f_t = []

    for a in range(absP):

        y_p_f_t.append([])

        for b in range(absF):

            y_p_f_t[a].append([])

            for g in range(T+1): #from 0 to T

                varname = "y_" + str(a) + "_" + str(b) + "_" + str(g)

                ally_p_f_t.append(varname)

                y_p_f_t[a][b].append(varname)

    c.variables.add(names=ally_p_f_t, lb=[0] * len(ally_p_f_t),

                    ub=[1] * len(ally_p_f_t),

                    types=["B"] * len(ally_p_f_t))


    # assignment variables: z_p_f_t_i_j
    allz_p_f_t_i_j = []
    for a in range(absP):
        z_p_f_t_i_j.append([])
        for b in range(absF):
            z_p_f_t_i_j[a].append([])
            for g in range(T+1):
                z_p_f_t_i_j[a][b].append([])
                for d in range(absV):
                    z_p_f_t_i_j[a][b][g].append([])
                    for e in range(absV):
                        varname = "z_" + str(a) + "_" + str(b) + "_" + str(g) + "_" + str(d) + "_" + str(e)
                        allz_p_f_t_i_j.append(varname)
                        z_p_f_t_i_j[a][b][g][d].append(varname)
    c.variables.add(names=allz_p_f_t_i_j, lb=[0] * len(allz_p_f_t_i_j),
                    ub=[1] * len(allz_p_f_t_i_j),
                    types=["B"] * len(allz_p_f_t_i_j))

    # assignment variables: w_p_f_t_i_j
    allw_p_f_t_i_j = []
    for a in range(absP):
        w_p_f_t_i_j.append([])
        for b in range(absF):
            w_p_f_t_i_j[a].append([])
            for g in range(T+1):
                w_p_f_t_i_j[a][b].append([])
                for d in range(absV):
                    w_p_f_t_i_j[a][b][g].append([])
                    for e in range(absV):
                        varname = "w_" + str(a) + "_" + str(b) + "_" + str(g) + "_" + str(d) + "_" + str(e)
                        allw_p_f_t_i_j.append(varname)
                        w_p_f_t_i_j[a][b][g][d].append(varname)
    c.variables.add(names=allw_p_f_t_i_j, lb=[0] * len(allw_p_f_t_i_j),
                    ub=[1] * len(allw_p_f_t_i_j),
                    types=["B"] * len(allw_p_f_t_i_j))

    # assignment variables: r_p_t_i_j
    allr_p_t_i_j= []
    for a in range(absP):
        r_p_t_i_j.append([])
        for b in range(T+1):
            r_p_t_i_j[a].append([])
            for g in range(absV):
                r_p_t_i_j[a][b].append([])
                for d in range(absV):
                    varname = "r_" + str(a) + "_" + str(b) + "_" + str(g) + "_" + str(d)
                    allr_p_t_i_j.append(varname)
                    r_p_t_i_j[a][b][g].append(varname)
    c.variables.add(names=allr_p_t_i_j, lb=[0] * len(allr_p_t_i_j),
                    ub=[1] * len(allr_p_t_i_j),
                    types=["B"] * len(allr_p_t_i_j))


    # assignment variables: U_t

    allU_t = []
    for a in range(T): # U(T) is objective function so I defined later
        varname = "U_" + str(a)
        allU_t.append(varname)
        U_t.append(varname)
    c.variables.add(names=allU_t, lb=[0] * len(allU_t),
                    types=["I"] * len(allU_t))

    # assignment variables: k_p_t
    allk_p_t = []
    for a in range(absP):
        k_p_t.append([])
        for b in range(T+1): #from 0 to T
            varname = "k_" + str(a) + "_" + str(b)
            allk_p_t.append(varname)
            k_p_t[a].append(varname)
    c.variables.add(names=allk_p_t, lb=[0] * len(allk_p_t),
                    ub=[1] * len(allk_p_t),
                    types=["B"] * len(allk_p_t))


    # # assignment variables: l_p_t
    # alll_p_t = []
    # for a in range(absP):
    #     l_p_t.append([])
    #     for b in range(T+1):
    #         varname = "l_" + str(a) + "_" + str(b)
    #         alll_p_t.append(varname)
    #         l_p_t[a].append(varname)
    # c.variables.add(names=alll_p_t, lb=[0] * len(alll_p_t),
    #                 ub=[1] * len(alll_p_t),
    #                 types=["B"] * len(alll_p_t))

    # fujun: the objective function in jornal paper is 

    # defferent with the one in GLPK code. In jornal paper,

    # t is from 1 to T, but it is from 0 to T - 1 in GLPK code.


    # assignment variables: h_p_t
    #fujun: I changed
    allhml_p_t = []
    for a in range(absP):
        h_p_t.append([])
        for b in range(T): #from 0 to T - 1,
            varname = "h_" + str(a) + "_" + str(b)
            allhml_p_t.append(varname)
            h_p_t[a].append(varname)

    # assignment variables: m_p_
    ##from 0 to T - 1, fujun: I changed,
    for a in range(absP):
        m_p_t.append([])
        for b in range(T):
            varname = "m_" + str(a) + "_" + str(b)
            allhml_p_t.append(varname)
            m_p_t[a].append(varname)

    # assignment variables: l_p_t
    for a in range(absP):
        l_p_t.append([])
        for b in range(T):
            varname = "l_" + str(a) + "_" + str(b)
            allhml_p_t.append(varname)
            l_p_t[a].append(varname)


    # assignment objectivefunction (1)
    varname = "U_" + str(T)
    U_t.append(varname)
    c.variables.add(names=[U_t[T]], lb=[0],
                    ub=[cplex.infinity],
                    types=["I"],
                    obj = [1])
    c.variables.add(names=allhml_p_t, lb=[0] * len(allhml_p_t),
                    ub=[1] * len(allhml_p_t),
                    types=["B"] * len(allhml_p_t),
                    obj = [epsilon] * len(allhml_p_t))

    # fujun: assignment variables: h_p_t[p][T] and m_p_t[p][T] and l_p_t[p][T]
    allhml_p_T = []
    for a in range(absP):
        varname = "m_" + str(a) + "_" + str(T)
        allhml_p_T.append(varname)
        m_p_t[a].append(varname)
        varname = "h_" + str(a) + "_" + str(T)
        allhml_p_T.append(varname)
        h_p_t[a].append(varname)
        varname = "l_" + str(a) + "_" + str(b)
        allhml_p_T.append(varname)
        l_p_t[a].append(varname)
    c.variables.add(names=allhml_p_T, lb=[0] * len(allhml_p_T),
                    ub=[1] * len(allhml_p_T),
                    types=["B"] * len(allhml_p_T))      



############################################################################



    #constrains
    # x_{p}^{f}(0)=1(2) Init
    for a in range(absP):
        thevars = [x_p_f_t[a][f_p_init[a]][0]]
        thecoefs = [1]
        c.linear_constraints.add(
            lin_expr=[cplex.SparsePair(thevars, thecoefs)],
            senses=["E"],
            rhs=[1])

    # k_{p}(0)=k_{p}^{\rm init}(3) Init2
    for a in range(absP):
        thevars = [k_p_t[a][0]] # fujun: I changed. KEY POINT
        thecoefs = [1]
        c.linear_constraints.add(
            lin_expr=[cplex.SparsePair(thevars, thecoefs)],
            senses=["E"],
            rhs=[k_p_init[a]]) # fujun: I changed, as k_p_init[a] is not a varible. KEY POINT.


    # r_{p}(0,i,j) = r_{p}^{\rm init}(i,j) OK
    for a in range(absP):
        for g in range(absV):
            for d in range(absV):
                if [g,d] in E:
                    thevars = [r_p_t_i_j[a][0][g][d]]
                    thecoefs = [1]
                    c.linear_constraints.add(
                        lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                        senses=["E"],
                        rhs=[r_p_i_j_init[a][g][d]])
    

    # \sum_{f\in F}x_{p}^{f}(t)=1(4) Alloc
    for a in range(absP):
        for b in range (T+1): #from 0 to T
            thevars = []
            thecoefs = []
            for g in range(absF):
                thevars.append(x_p_f_t[a][g][b])
                thecoefs.append(1)
            c.linear_constraints.add(
                lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                senses=["E"],
                rhs=[1])

    # x_{p}^{f}(t)\leq y_{p}^{f^{\prime}}(t)(5) Cons1
    for a in range(absP):
       for b in range(T+1): #from 0 to T
           for g in range(absF - n_p[a] + 1): #from 0 to F-n(p)
               for d in range(g,g + n_p[a]): # from f to f+n(p)-1
                   thevars = [x_p_f_t[a][g][b], y_p_f_t[a][d][b]]
                   thecoefs = [1, -1]
                   c.linear_constraints.add(
                       lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                       senses=["L"],
                       rhs=[0])

    # x_{p}^{f}(t)=0(6) Cons2
    for a in range(absP):
       for b in range(T+1): #from 0 to T
           for g in range(absF - n_p[a] + 1, absF): #from F-n(p)+1 to F-1
               thevars = [x_p_f_t[a][g][b]]
               thecoefs = [1]
               c.linear_constraints.add(
                   lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                   senses=["E"],
                   rhs=[0])

    # #\sum_{p\in P}z_{p}^{f}(t,i,j) + \sum_{p\in P}z_{p}^{f}(t,j,i) \leq 1 Nondirectional
    # for b in range(T+1): #from 0 to T
    #     for g in range(absV):
    #         for d in range(absV):
    #             if [g,d] in E:
    #                 for e in range(absF):
    #                     thevars = []
    #                     thecoefs = []
    #                     for a in range(absP):
    #                         thevars.extend([z_p_f_t_i_j[a][e][b][g][d],z_p_f_t_i_j[a][e][b][d][g]])
    #                         thecoefs.extend([1,1])
    #                     c.linear_constraints.add(
    #                         lin_expr=[cplex.SparsePair(thevars, thecoefs)],
    #                         senses=["L"],
    #                         rhs=[1])

        #\sum_{p\in P}z_{p}^{f}(t,i,j) \leq 1 Unidirectional
    for b in range(T+1): #from 0 to T
        for g in range(absV):
            for d in range(absV):
                if [g,d] in E:
                    for e in range(absF):
                        thevars = []
                        thecoefs = []
                        for a in range(absP):
                            thevars.append(z_p_f_t_i_j[a][e][b][g][d])
                            thecoefs.append(1)
                        c.linear_constraints.add(
                            lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                            senses=["L"],
                            rhs=[1])


    #z_{p}^{f}(t,i,j) \leq y_{p}^{f}(t)
    for a in range(absP):
        for b in range(T+1):
            for g in range(absV):
                for d in range(absV):
                    if [g,d] in E:
                        for e in range(absF):
                            thevars = [z_p_f_t_i_j[a][e][b][g][d],y_p_f_t[a][e][b]]
                            thecoefs = [1, -1]
                            c.linear_constraints.add(
                                lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                                senses=["L"],
                                rhs=[0])

    # z_{p}^{f}(t,i,j) \leq r_{p}(t,i,j)
    for a in range(absP):
        for b in range(T+1):
            for g in range(absV):
                for d in range(absV):
                    if [g,d] in E:
                        for e in range(absF):
                            thevars = [z_p_f_t_i_j[a][e][b][g][d],r_p_t_i_j[a][b][g][d]]
                            thecoefs = [1, -1]
                            c.linear_constraints.add(
                                lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                                senses=["L"],
                                rhs=[0])

    # z_{p}^{f}(t,i,j) \geq y_{p}^{f}(t) + r_{p}(t,i,j) -1
    for a in range(absP):
        for b in range(T+1):
            for g in range(absV):
                for d in range(absV):
                    if [g,d] in E:
                        for e in range(absF):
                            thevars = [z_p_f_t_i_j[a][e][b][g][d],y_p_f_t[a][e][b],r_p_t_i_j[a][b][g][d]]
                            thecoefs = [1, -1, -1]
                            c.linear_constraints.add(
                                lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                                senses=["G"],
                                rhs=[-1])

    # \sum_{j:(i,j)\in E}r_{p}(t,i,j) - \sum_{j:(i,j)\in E}r_{p}(t,j,i) =1
    for a in range(absP):
       for b in range(T+1):
           for g in range(absV):
               if g == o_p[a]:
                   thevars = []
                   thecoefs = []
                   for d in range(absV):
                       if [g,d] in E:
                           thevars.append(r_p_t_i_j[a][b][g][d])
                           thecoefs.append(1)
                           thevars.append(r_p_t_i_j[a][b][d][g])
                           thecoefs.append(-1)
                   c.linear_constraints.add(
                       lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                       senses=["E"],
                       rhs=[1])

    # \sum_{j:(i,j)\in E}r_{p}(t,i,j) - \sum_{j:(i,j)\in E}r_{p}(t,j,i) =0
    for a in range(absP):
       for b in range(T+1):
           for g in range(absV):
               if g != o_p[a] and g != d_p[a]:
                   thevars = []
                   thecoefs = []
                   for d in range(absV):
                       if [g,d] in E:
                           thevars.append(r_p_t_i_j[a][b][g][d])
                           thecoefs.append(1)
                           thevars.append(r_p_t_i_j[a][b][d][g])
                           thecoefs.append(-1)
                   c.linear_constraints.add(
                       lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                       senses=["E"],
                       rhs=[0])

    # # r_{p}(t,i,j)+r_{p}(t,j,i)+r_{p^{\prime}}(t,i,j)+r_{p^{\prime}}(t,j,i) \leq 1 Nondirectional
    # for a in range(absP):
    #    for aa in range(absP):
    #        if a != aa and s_p[a] == s_p[aa]:
    #            for b in range(T+1):
    #                for g in range(absV):
    #                    for d in range(absV):
    #                     if [g,d] in E:
    #                            thevars = [r_p_t_i_j[a][b][g][d],r_p_t_i_j[a][b][d][g],r_p_t_i_j[aa][b][g][d],r_p_t_i_j[aa][b][d][g]]
    #                            thecoefs = [1,1,1,1]
    #                            c.linear_constraints.add(
    #                                lin_expr=[cplex.SparsePair(thevars, thecoefs)],
    #                                senses=["L"],
    #                                rhs=[1])

    # r_{p}(t,i,j)+r_{p^{\prime}}(t,i,j) \leq 1 Unidirectional
    for a in range(absP):
       for aa in range(absP):
           if a != aa and s_p[a] == s_p[aa]:
               for b in range(T+1):
                   for g in range(absV):
                       for d in range(absV):
                        if [g,d] in E:
                               thevars = [r_p_t_i_j[a][b][g][d], r_p_t_i_j[aa][b][g][d]]
                               thecoefs = [1,1]
                               c.linear_constraints.add(
                                   lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                                   senses=["L"],
                                   rhs=[1])


    # f \times y_{p}^{f}(t) \leq U(t)(8) MaxInd
    for a in range(absP):
       for b in range(T+1): #from 0 to T
           for g in range(absF):
                   thevars = [y_p_f_t[a][g][b], U_t[b]]
                   thecoefs = [g, -1]
                   c.linear_constraints.add(
                       lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                       senses=["L"],
                       rhs=[0])


    # # z_{p}^{f}(t,i,j)+z_{p}^{f}(t,j,i)+w_{p}^{f}(t+1,j,i) \leq 1 Nondirectional
    # for a in range(absP):
    #     for b in range(T):
    #         for g in range(absV):
    #             for d in range(absV):
    #                 if [g,d] in E:
    #                     for e in range(absF):
    #                         thevars = [z_p_f_t_i_j[a][e][b][g][d],z_p_f_t_i_j[a][e][b][d][g], w_p_f_t_i_j[a][e][b+1][g][d],w_p_f_t_i_j[a][e][b+1][d][g]]
    #                         thecoefs = [1,1,1,1]
    #                         c.linear_constraints.add(
    #                             lin_expr=[cplex.SparsePair(thevars, thecoefs)],
    #                             senses=["L"],
    #                             rhs=[1])

    # z_{p}^{f}(t,i,j)+w_{p}^{f}(t+1,i,j) \leq 1 Unidirectional
    for a in range(absP):
        for b in range(T):
            for g in range(absV):
                for d in range(absV):
                    if [g,d] in E:
                        for e in range(absF):
                            thevars = [z_p_f_t_i_j[a][e][b][g][d], w_p_f_t_i_j[a][e][b+1][g][d]]
                            thecoefs = [1,1]
                            c.linear_constraints.add(
                                lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                                senses=["L"],
                                rhs=[1])

    # w_{p}^{f}(t+1,i,j) \leq r_{p}(t,i,j)
    for a in range(absP):
        for b in range(T):
            for g in range(absV):
                for d in range(absV):
                    if [g,d] in E:
                        for e in range(absF):
                            thevars = [w_p_f_t_i_j[a][e][b+1][g][d],r_p_t_i_j[a][b][g][d]]
                            thecoefs = [1,-1]
                            c.linear_constraints.add(
                                lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                                senses=["L"],
                                rhs=[0])

    # w_{p}^{f}(t+1,i,j) \leq \sum_{p^{\prime}\in P:p^{\prime}\neq p}z_{p^{\prime}}^{f}(t+1,i,j)
    for b in range(T):
        for g in range(absV):
            for d in range(absV):
                if [g,d] in E:
                    for e in range(absF):
                        for a in range(absP):
                            thevars = []
                            thecoefs = []
                            thevars.append(w_p_f_t_i_j[a][e][b+1][g][d])
                            thecoefs.append(1)
                            for aa in range(absP):
                                if a != aa:
                                    thevars.append(z_p_f_t_i_j[aa][e][b+1][g][d])
                                    thecoefs.append(-1)
                            c.linear_constraints.add(
                                lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                                senses=["L"],
                                rhs=[0])

    # w_{p}^{f}(t+1,i,j) \geq r_{p}(t,i,j) + \sum_{p^{\prime}\in P:p^{\prime}\neq p}z_{p^{\prime}}^{f}(t+1,i,j) -1
    for b in range(T):
        for g in range(absV):
            for d in range(absV):
                if [g,d] in E:
                    for e in range(absF):
                        for a in range(absP):
                            thevars = []
                            thecoefs = []
                            thevars.extend([w_p_f_t_i_j[a][e][b+1][g][d],r_p_t_i_j[a][b][g][d]])
                            thecoefs.extend([1, -1])
                            for aa in range(absP):
                                if a != aa:
                                    thevars.append(z_p_f_t_i_j[aa][e][b+1][g][d])
                                    thecoefs.append(-1)
                            c.linear_constraints.add(
                                lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                                senses=["G"],
                                rhs=[-1])



    # |x_{p}^{f}(t)-x_{p}^{f}(t+1)| \leq m_{p}(t)(9) Move1,2

    for a in range(absP):

       for b in range(T): #from 0 to T-1

           for g in range(absF):

               thevars = [x_p_f_t[a][g][b], x_p_f_t[a][g][b+1], m_p_t[a][b]]

               thecoefs = [1, -1, 1]

               c.linear_constraints.add(

                   lin_expr=[cplex.SparsePair(thevars, thecoefs)],

                   senses=["G"],

                   rhs=[0])

    for a in range(absP):

        for b in range(T): #from 0 to T-1

           for g in range(absF):

               thevars = [x_p_f_t[a][g][b], x_p_f_t[a][g][b+1], m_p_t[a][b]]

               thecoefs = [1, -1, -1]

               c.linear_constraints.add(

                   lin_expr=[cplex.SparsePair(thevars, thecoefs)],

                   senses=["L"],

                   rhs=[0])



    # |r_{p}(t,i,j)-r_{p}(t+1,i,j)| \leq l_{p}(t)
    for a in range(absP):
        for b in range(T):
            for g in range(absV):
                for d in range(absV):
                    if [g,d] in E:
                        for e in range(absF):
                            thevars = [r_p_t_i_j[a][b][g][d],r_p_t_i_j[a][b+1][g][d],l_p_t[a][b]]
                            thecoefs = [1, -1, 1]
                            c.linear_constraints.add(
                                lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                                senses=["G"],
                                rhs=[0])
    for a in range(absP):
        for b in range(T):
            for g in range(absV):
                for d in range(absV):
                    if [g,d] in E:
                        for e in range(absF):
                            thevars = [r_p_t_i_j[a][b][g][d],r_p_t_i_j[a][b+1][g][d],l_p_t[a][b]]
                            thecoefs = [1, -1, -1]
                            c.linear_constraints.add(
                                lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                                senses=["L"],
                            rhs=[0])
        

    # k_{p}(t)+m_{p}(t) \leq 1(11) fujun: PBmove

    for a in range(absP):

       for b in range(T+1): #from 0 to T

           thevars = [k_p_t[a][b], m_p_t[a][b]]

           thecoefs = [1, 1]

           c.linear_constraints.add(

               lin_expr=[cplex.SparsePair(thevars, thecoefs)],

               senses=["L"],

               rhs=[1])



    # k_{p}(t+1)+m_{p}(t) \leq 1(12) fujun: PBmove2

    for a in range(absP):

       for b in range(T): #from 0 to T-1

           thevars = [k_p_t[a][b+1], m_p_t[a][b]]

           thecoefs = [1, 1]

           c.linear_constraints.add(

               lin_expr=[cplex.SparsePair(thevars, thecoefs)],

               senses=["L"],

               rhs=[1])

    # k_{p}(t)+l_{p}(t) \leq 1
    for a in range(absP):
        for b in range(T+1):
            thevars = [k_p_t[a][b], l_p_t[a][b]]
            thecoefs = [1, 1]
            c.linear_constraints.add(
                lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                senses=["L"],
                rhs=[1])

    # k_{p}(t+1)+l_{p}(t) \leq 1(12) fujun: PBmove2
    for a in range(absP):
       for b in range(T): #from 0 to T-1
           thevars = [k_p_t[a][b+1], l_p_t[a][b]]
           thecoefs = [1, 1]
           c.linear_constraints.add(
               lin_expr=[cplex.SparsePair(thevars, thecoefs)],
               senses=["L"],
               rhs=[1])


    # k_{p}(t)+k_{p{\prime}}(t)=1(13) Toggle1

    for a in range(absP):

       for aa in range(absP):

           if a != aa and s_p[a] == s_p[aa]:

               for b in range(T+1): #from 0 to T

                   thevars = [k_p_t[a][b], k_p_t[aa][b]]

                   thecoefs = [1, 1]

                   c.linear_constraints.add(

                       lin_expr=[cplex.SparsePair(thevars, thecoefs)],

                       senses=["E"],

                       rhs=[1])





    # k_{p}(t)-k_{p}(t+1) \leq h_{p}(t)(14) Toggle2

    for a in range(absP):

       for b in range(T): #from 0 to T-1

           thevars = [k_p_t[a][b], k_p_t[a][b+1], h_p_t[a][b]]

           thecoefs = [1, -1, -1]

           c.linear_constraints.add(

               lin_expr=[cplex.SparsePair(thevars, thecoefs)],

               senses=["L"],

               rhs=[0])



def Rerouting():

    c = cplex.Cplex()


    # print(datetime.datetime.now())
    setupproblem(c)

    # c.write("../result/ssrr_lno.lp") # fujun: you can debug by reading lp file.
    # print(datetime.datetime.now())
    c.solve()
    # print(datetime.datetime.now())
    sol = c.solution

    # solution.get_status() returns an integer code
    # print("Solution status = ", sol.get_status(), ":", end=' ')
    # the following line prints the corresponding string
    # print(sol.status[sol.get_status()])



    if sol.is_primal_feasible():
        f = open('./../result/ssrr_lno_result.txt', 'w')
        h_total = 0
        m_total = 0
        for a in range(absP):
          for b in range(T):
            h_total += int(sol.get_values(h_p_t[a][b]))
            m_total += int(sol.get_values(m_p_t[a][b]))
        f.write("sum_h_p_t =" + str(h_total) + "\n")
        f.write("sum_m_p_t =" + str(m_total) + "\n" + "\n")

        for a in range(absP):
            f.write("r_" + str(a) + "_t_i_j = \n")
            # print("r_" + str(a) + "_t_i_j = \n")
            for g in range(absV):
                for d in range(absV):
                    if round(sol.get_values(r_p_t_i_j[a][T][g][d])) == 1:
                        f.write(str(a) + " " + str(g) + " " + str(d) + " " + str(0) + "\n")
                        # print(str(a) + " " + str(g) + " " + str(d) + " " + str(0))
            f.write(str(0) + " " + str(0) + " " + str(0) + " " + str(1) + "\n")#パスの最後を示す判別子
            # print(str(0) + " " + str(0) + " " + str(0) + " " + str(1))#パスの最後を示す判別子
        f.write("\n")
        # print("\n")

        f.write("x_p_f_T =" + "\n")
        for a in range(absP):
          for b in range(absF):
            # print(x_p_f_t[a][b][T] + " = " + str(sol.get_values(x_p_f_t[a][b][T])))
            if round(sol.get_values(x_p_f_t[a][b][T])) == 1:
              f.write(str(a)+ " " + str(b) + "\n")

        f.close()

        # print("Solution value  = ", sol.get_objective_value())



        # print("U(T) = %d" %

        #         sol.get_values(U_t[T]))

        # print()


        # for a in range(absP):
        #     for d in range(T+1):
        #         print()
        #         for b in range(absV):
        #             unmaped_list = []
        #             for g in range(absV):
        #                 unmaped_list.append(int(sol.get_values(r_p_t_i_j[a][d][b][g])))
        #             maped_list = map(str, unmaped_list)
        #             mojiretsu = ' '.join(maped_list)
        #             print("r_%d(%d,%d) = %s"% (a,d,b, mojiretsu))
        #         print()
        #     print('------------------------------------')

        # for a in range(absP):

        #     unmaped_list = []
        #     for b in range(T+1):
        #         unmaped_list.append(int(sol.get_values(k_p_t[a][b])))
        #     maped_list = map(str, unmaped_list)
        #     mojiretsu = ' '.join(maped_list)
        #     print("k_%d = %s"% (a, mojiretsu))
        # print()

        # for a in range(absP):

        #     unmaped_list = []
        #     for b in range(T):
        #         unmaped_list.append(int(sol.get_values(h_p_t[a][b])))
        #     maped_list = map(str, unmaped_list)
        #     mojiretsu = ' '.join(maped_list)
        #     print("h_%d = %s"% (a, mojiretsu))
        # print()

        # for a in range(absP):
        #     unmaped_list = []
        #     for b in range(T):
        #         unmaped_list.append(int(sol.get_values(m_p_t[a][b])))
        #     maped_list = map(str, unmaped_list)
        #     mojiretsu = ' '.join(maped_list)
        #     print("m_%d = %s"% (a, mojiretsu))
        # print()


        # for a in range(absP):
        #     unmaped_list = []
        #     for b in range(T):
        #         unmaped_list.append(int(sol.get_values(l_p_t[a][b])))
        #     maped_list = map(str, unmaped_list)
        #     mojiretsu = ' '.join(maped_list)
        #     print("l_%d = %s"% (a, mojiretsu))
        # print()



        # for a in range(absP):

        #     unmaped_list = []
        #     for b in range(absF):

        #         unmaped_list.append(int(sol.get_values(x_p_f_t[a][b][T])))
        #     maped_list = map(str, unmaped_list)
        #     mojiretsu = ' '.join(maped_list)
        #     print("x_%d(T) = %s"% (a, mojiretsu))
        # print()

        # for a in range(absP):

        #     unmaped_list = []
        #     for b in range(absF):
        #         unmaped_list.append(int(sol.get_values(y_p_f_t[a][b][T])))
        #     maped_list = map(str, unmaped_list)
        #     mojiretsu = ' '.join(maped_list)
        #     print("y_%d(T) = %s"% (a,mojiretsu))
        # print()


    else:
        print("No solution available.")

if __name__ == "__main__":

    Rerouting()

