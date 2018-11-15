# -*- coding: utf-8 -*-

from __future__ import print_function
from math import fabs
import cplex
import sys
import math
import datetime

f = open('./../result/ssr_lno_output.txt')
# print(datetime.datetime.now())

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

#number of allowed transitions
line = f.readline()
line = line.split("=")
T = int(line[1])
# print("T = ", T)

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
line = f.readline()#set P_e := 

P_e = []
for e in range (absE):
  P_e.append([])

line = f.readline()
while line[0] != ";":
  line = line.split()
  P_e[int(line[1])].append(int(line[0]))
  line = f.readline()

# print("P_e = ", P_e)

f.close()

# epsilon
epsilon = 0.0001

#problem variables
x_p_f_t = []
y_p_f_t = []
m_p_t = []
U_t = []
h_p_t = []
k_p_t = []

################################################################

def setupproblem(c):
    # ログを表示しない
    c.set_log_stream(None)
    c.set_error_stream(None)
    c.set_warning_stream(None)
    c.set_results_stream(None)
    c.objective.set_sense(c.objective.sense.minimize)

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

    # assignment variables: U_t
    allU_t = []
    for a in range(T): # U(T) is objective function so I defined later
        varname = "U_" + str(a)
        allU_t.append(varname)
        U_t.append(varname)
    c.variables.add(names=allU_t, lb=[0] * len(allU_t),
                    types=["I"] * len(allU_t))

    # assignment variables: k_p_
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

    # assignment variables: h_p_t
    # fujun: the objective function in jornal paper is 
    # defferent with the one in GLPK code. In jornal paper,
    # t is from 1 to T, but it is from 0 to T - 1 in GLPK code.
    allhm_p_t = []
    for a in range(absP):
        h_p_t.append([])
        for b in range(T): #from 0 to T - 1, fujun: I changed # sawa:I changed
            varname = "h_" + str(a) + "_" + str(b)
            allhm_p_t.append(varname)
            h_p_t[a].append(varname)

    # assignment variables: m_p_t
    for a in range(absP):
        m_p_t.append([])
        for b in range(T): #from 0 to T - 1, fujun: I changed
    # sawa:I changed
            varname = "m_" + str(a) + "_" + str(b)
            allhm_p_t.append(varname)
            m_p_t[a].append(varname)

    # assignment objectivefunction (1)
    varname = "U_" + str(T)
    U_t.append(varname)
    c.variables.add(names=[U_t[T]], lb=[0],
                    ub=[cplex.infinity],
                    types=["I"],
                    obj = [1])
    c.variables.add(names=allhm_p_t, lb=[0] * len(allhm_p_t),
                    ub=[1] * len(allhm_p_t),
                    types=["B"] * len(allhm_p_t),
                    obj = [epsilon] * len(allhm_p_t))

    # fujun: assignment variables: h_p_t[p][T] and m_p_t[p][T]
    allhm_p_T = []
    for a in range(absP):
        varname = "m_" + str(a) + "_" + str(T)
        allhm_p_T.append(varname)
        m_p_t[a].append(varname)
        varname = "h_" + str(a) + "_" + str(T)
        allhm_p_T.append(varname)
        h_p_t[a].append(varname)
    c.variables.add(names=allhm_p_T, lb=[0] * len(allhm_p_T),
                    ub=[1] * len(allhm_p_T),
                    types=["B"] * len(allhm_p_T),
                    obj = [0] * len(allhm_p_T))       

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

    #\sum_{p\in P_{e}}y_{p}^{f}(t) \leq 1(7) SlotCap
    for b in range(T+1): #from 0 to T
       for g in range(absF):
           for d in range(absE):
               thevars = []
               thecoefs = []
               for a in range(absP):
                   if a in P_e[d]: # fujun: I changed according P_e sawa:I changed again
                       thevars.append(y_p_f_t[a][g][b])
                       thecoefs.append(1)
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

    # y_{p}^{f}(t)+\sum _{p^{\prime}\in P_{e}:p^{\prime}\neq p}y_{p^{\prime}}^{f}(t+1) \leq 1(10) PConf
    for b in range(T): #from 0 to T-1
       for g in range(absF):
          for d in range(absE):
              for a in range(absP):
                  thevars = []
                  thecoefs = []
                  if a in P_e[d]: # a means f
                      thevars.append(y_p_f_t[a][g][b])
                      thecoefs.append(1)
                      for aa in range(absP):
                          if aa in P_e[d] and aa != a: # aa means f'
                              thevars.append(y_p_f_t[aa][g][b+1])
                              thecoefs.append(1)
                      c.linear_constraints.add(
                          lin_expr=[cplex.SparsePair(thevars, thecoefs)],
                          senses=["L"],
                          rhs=[1])

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
    setupproblem(c)
    # c.write("./../result/ssr_lno.lp") # fujun: you can debug by reading lp file.
    c.solve()

    sol = c.solution

    # print("変数の数 = ", c.variables.get_num())    

    # print()

    # # solution.get_status() returns an integer code
    # print("Solution status = ", sol.get_status(), ":", end=' ')
    # # the following line prints the corresponding string
    # print(sol.status[sol.get_status()])

    if sol.is_primal_feasible():
      # print("Solution value  = ", sol.get_objective_value())
      # print()

      f = open('./../result/ssr_lno_result.txt', 'w')
      h_total = 0
      m_total = 0
      for a in range(absP):
        for b in range(T):
          h_total += int(sol.get_values(h_p_t[a][b]))
          m_total += int(sol.get_values(m_p_t[a][b]))
      f.write("sum_h_p_t =" + str(h_total) + "\n")
      f.write("sum_m_p_t =" + str(m_total) + "\n" + "\n")

      f.write("x_p_f_T =" + "\n")
      for a in range(absP):
        for b in range(absF):
          # print(x_p_f_t[a][b][T] + " = " + str(sol.get_values(x_p_f_t[a][b][T])))
          if round(sol.get_values(x_p_f_t[a][b][T])) == 1:
            f.write(str(a)+ " " + str(b) + "\n")

      # print("U(T) = %d" %
      #         sol.get_values(U_t[T]))
      # print()

      ####################################
      #    int(1.0)が0になる場合がある
      ####################################
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

      # for g in range(T+1):
      #   for a in range(absP):
      #       unmaped_list = []
      #       for b in range(absF):
      #           unmaped_list.append(int(sol.get_values(x_p_f_t[a][b][g])))
      #       maped_list = map(str, unmaped_list)
      #       mojiretsu = ' '.join(maped_list)
      #       print("x_%d(%d) = %s"% (a, g, mojiretsu))
      #   print()

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

