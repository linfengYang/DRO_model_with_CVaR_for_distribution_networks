import gurobipy as gp
import numpy as np
from gurobipy import GRB
from scipy.stats import norm
import math
def inte(x):
    return math.exp(-x**2)

def funcpd(beta ):
    Pd = [3980, 3610, 3340, 2990, 2690, 2920, 3260, 3550, 3980, 4290, 4470, 4610, 4730, 4840, 4930, 4710, 4520, 4830,
          5090, 5385, 4980, 4610, 4410, 4290]  # KW
    Qd = [0] * 24
    # 33节点，总pd3715，总qd2300
    for i in range(24): Qd[i] = Pd[i] * 2300 / 3715
    expectedpd = max(Pd)  # 假定最大负荷
    sub = 1
    T = 24
    with open('france_original_rt.txt', 'r') as file:
        lines = file.readlines()
        modified_lines = [line[18:] for line in lines]
        lines = [line.split() for line in modified_lines]
        while [] in lines:
            lines.remove([])
        for i in range(len(lines)):
            lines[i] = lines[i][3]
            lines[i] = lines[i].replace(',', '.')
            lines[i] = float(lines[i])
    history_rt = lines[::2]
    with open('france_da.txt', 'r') as file:
        lines = file.readlines()
        modified_lines = [line[8:] for line in lines]
        lines = [line.split() for line in modified_lines]
        for line in lines:
            for i in range(len(line)):
                line[i] = line[i].replace(',', '.')
                line[i] = float(line[i])
        while [] in lines:
            lines.remove([])

    length = (len(lines) - 24) // 3 + 6
    history = np.zeros((length, 24 * sub))
    for i in range(length):
        for j in range(T):
            history[i, j] = history_rt[i * T + j]
    history = np.array(history)
    def f_pd(Pd,Qd,history,time,beta,delta,flag):
        with open('bus33.txt', 'r') as file:
            f = file.readline()
            bus = []
            while f:
                a = []
                b = f.split()
                for i in b[:-1]:
                    a.append(float(i))
                bus.append(a)
                f = file.readline()
        bus = np.array(bus)
        totalbuspd = sum(bus[:, 2])
        if flag:
            bus=bus[1:2]
        Vbase = 12.66 # 基准电压 kv
        MVAbase = 10  # 基准视在功率 MVA

        T = time
        sub = 1
        ST = sub * T
        num_node = len(bus)  # 节点数
        #扩大负荷，避免精度误差
        for i in range(T):
            Pd[i] *= 1000*totalbuspd/expectedpd
            Qd[i] *= 1000*totalbuspd/expectedpd
            Pd[i] /= (MVAbase * 1e3)
            Qd[i] /= (MVAbase * 1e3)
        sp = sum(bus[:, 2])
        sq = sum(bus[:, 3])
        pd_factor = [i[2] / sp for i in bus]  # 获得每个节点的负荷因子，即负荷占总负荷的比例,
        qd_factor = [i[3] / sq for i in bus]
        if flag:
            pd = [[i * j /3715*100 for i in Pd] for j in pd_factor]  # 获得每个节点各个时段的有功负荷,
            qd = [[i * j /2300*60 for i in Qd] for j in qd_factor]  # 获得每个节点各个时段的无功负荷,
        else:
            pd = [[i * j  for i in Pd] for j in pd_factor]
            qd = [[i * j  for i in Qd] for j in qd_factor]
        pd = np.array(pd)
        qd = np.array(qd)
        history = np.array(history)
        scale = [[pd[i][t] * 0.05 for t in range(T)] for i in range(num_node)] #标准差
        scale = np.array(scale)

        m1 = gp.Model()
        deltapd = m1.addVars(num_node, T, lb=-10, vtype=GRB.CONTINUOUS, name='deltapd')  # 负荷波动
        deltaYoe = m1.addVars(num_node, T, vtype=GRB.CONTINUOUS, name='deltaYoe')  # 负荷高估预期
        deltaYue = m1.addVars(num_node, T, vtype=GRB.CONTINUOUS, name='deltaYue')  # 负荷低估预期
        deltaYuesum = m1.addVars(sub, T, vtype=GRB.CONTINUOUS, name='deltaYuesum')  # 负荷低估之和
        m1.addConstrs(deltapd[i,t]==delta[t] *scale[i,t] for i in range(num_node) for t  in range(T))

        #高估低估预期
        m1.addConstrs(deltaYue[i, t]== scale[i, t] / math.sqrt(2 * math.pi) * math.exp(-(delta[t] *scale[i,t]) ** 2
                      / (2 * scale[i, t] ** 2))- deltapd[i, t]+deltapd[i,t]*norm.cdf(delta[t] *scale[i,t],0,scale[i,t])
                      for i in range(num_node) for t in range(T) if scale[i, t] != 0)
        m1.addConstrs(deltaYoe[i, t]== scale[i, t] / math.sqrt(2 * math.pi) * math.exp(-(delta[t] *scale[i,t]) ** 2
                     / (2 * scale[i, t] ** 2))+deltapd[i,t]*norm.cdf(delta[t] *scale[i,t],0,scale[i,t])
                      for i in range(num_node) for t in range(T) if scale[i, t] != 0)

        # 无负荷需求的点也没有负荷波动
        m1.addConstrs(deltaYoe[i, t] == 0 for i in range(num_node) for t in range(T) if scale[i, t] == 0)
        m1.addConstrs(deltaYue[i, t] == 0 for i in range(num_node) for t in range(T) if scale[i, t] == 0)

        m1.addConstrs(deltaYuesum[s, t] == gp.quicksum(deltaYue[i, t]  for i in range(num_node))/sub
                      for s in range(sub) for t in range(T))

        alpha = m1.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="alpha")  # 在β概率下，保证损失小于α的最小的α
        M,L1 = 366,15
        x_low = np.array([min(history[:, i]) for i in range(ST)])
        x_high = np.array([max(history[:, i]) for i in range(ST)])
        s = [x_low, x_high]
        C_ = []
        for i in range(ST):
            C_.extend([x_low[i] + (x_high[i] - x_low[i]) * l / L1 for l in range(L1 + 1)])
        C_ = np.array(C_)
        gama=[]
        for i in range(ST):
            for l in range(L1+1):
                gama.append( sum([max(history[j][i] - C_[i * (L1 + 1) + l], 0) for j in range(M)]) / M +1 )
                

        # 模糊集矩阵化
        C1 = np.diag([1] * ST)
        C2 = np.diag([-1] * ST)
        C3 = np.diag([0] * ST)
        C4 = [C1, C2]
        [C4.append(C3) for _ in range(2 * (L1 + 1))]
        C = np.concatenate(C4, axis=0)
        k = 0
        for i in range((2 + (L1 + 1)) * ST, (2 + 2 * (L1 + 1)) * ST):
            C[i, k] = 1
            if not (i - (2 + (L1 + 1)) * ST + 1) % (L1 + 1):
                k += 1
        D1 = np.diag([-1] * (L1 + 1) * ST)
        D = np.vstack((np.zeros((2 * ST, (L1 + 1) * ST)), D1, D1))
        d = np.zeros(((L1 + 1) * ST))
        d = np.hstack((d, np.zeros((2 * ST)), d))
        d[:ST] = x_high
        d[ST:2 * ST] = -x_low
        d[(2 + L1 + 1) * ST:] = C_
        G = 2 * (L1 + 1) * ST + 2 * ST

        #模型对偶化
        eta = m1.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="eta")
        omega = m1.addVars((L1 + 1) * ST, vtype=GRB.CONTINUOUS, name="omega")

        pi1 = m1.addVars(G, vtype=GRB.CONTINUOUS, name="pi1")
        m1.addConstr(-gp.quicksum(d[i] * pi1[i] for i in range(G)) + eta - alpha >= 0, name="c29")
        m1.addConstrs((gp.quicksum(C[j, i] * pi1[j] for j in range(G)) == 0 for i in range(ST)),
                      name="c30")
        m1.addConstrs((gp.quicksum(D[j, i] * pi1[j] for j in range(G)) == -omega[i] for i in range((L1 + 1) * ST)),
                      name="c31")

        pi2 = m1.addVars(G, vtype=GRB.CONTINUOUS, name="pi2")
        m1.addConstr(-gp.quicksum(d[i] * pi2[i] for i in range(G)) + eta + (beta * alpha
                      - 100 * gp.quicksum(deltaYoe[i, t] for i in range(num_node) for t in range(T)) ) / (1 - beta) >= 0, name="c16")
        m1.addConstrs((gp.quicksum(C[j, i] * pi2[j] for j in range(G)) ==
                       (deltaYuesum[i // T, i % T]) / (1 - beta) for i in range(ST)), name="c17")
        m1.addConstrs(
            (gp.quicksum(D[j, i] * pi2[j] for j in range(G)) == -omega[i] for i in range((L1 + 1) * ST)), name="c18")

        m1.Params.OutputFlag = 0
        m1.Params.MIPGap = 0
        m1.setObjective( eta + gp.quicksum(omega[i] * gama[i] for i in range((L1) * ST)))
        m1.update()
        m1.optimize()

        deltapd_copy = [[0 for t in range(T)] for i in range(num_node)]
        deltaYoe_copy = [[0 for t in range(T)] for i in range(num_node)]
        deltaYue_copy = [[0 for t in range(T)] for i in range(num_node)]
        deltaYuesum_copy = [[0 for t in range(T)] for i in range(sub)]
        for i in range(num_node):
            for t in range(T):
                deltapd_copy[i][t] = deltapd[i, t].X
                deltaYoe_copy[i][t] = deltaYoe[i, t].X
                deltaYue_copy[i][t] = deltaYue[i, t].X
        for i in range(sub):
            for t in range(T):
                deltaYuesum_copy[i][t] = deltaYuesum[i, t].X
        deltapd_copy = np.array(deltapd_copy)
        deltaYoe_copy = np.array(deltaYoe_copy)
        deltaYue_copy = np.array(deltaYue_copy)
        deltaYuesum_copy = np.array(deltaYuesum_copy)
        obj=m1.objVal
        m1.dispose()
        return deltapd_copy, deltaYoe_copy, deltaYue_copy, deltaYuesum_copy,obj
    TIMEpiece=[0]
    TIMEpiece.extend([1 for i in range(T)])
    b=history[:,1::T]
    TIMEstart=0
    TIMEend=TIMEpiece[0]
    num_node=33
    S_scale = []
    for numt in range(len(TIMEpiece)-1):
        TIMEstart+=TIMEpiece[numt]
        TIMEend+=TIMEpiece[numt+1]
        a = list(history[:, TIMEstart::T])
        for i in range(len(b)):
            for j in range(len(b[0])):
                for n in range(1,TIMEpiece[numt+1]):
                    c = history[:, TIMEstart + n::T]
                    a[i] = list(a[i])
                    a[i].insert(TIMEpiece[numt+1] * j + n, c[i][j])
        S = []
        beta_a=0.1
        if beta==0.95:beta_a=0.3
        if beta == 0.9: beta_a = 0.25
        if beta == 0.85: beta_a = 0.2
        if beta == 0.8: beta_a = 0.15
        DELTA_1 = [beta_a for t in range(T)]
        S_pd=[]
        for n in range(75):
            s = 0
            DELTA1 = DELTA_1[:]
            DELTA1[numt] = DELTA1[numt] + n * (0.02+(0.95-beta)**2) * DELTA_1[numt]
            temp=f_pd(Pd[TIMEstart:TIMEend], Qd[TIMEstart:TIMEend], a, TIMEpiece[numt+1], beta,DELTA1[numt:],1)
            S.append(temp[4])
            S_pd.append(temp[0])
        S_scale.append(S_pd[S.index(min(S))][0][0])
    for t in range(T):
        S_scale[t]=S_scale[t]/0.36954503*Pd[0]/Pd[t]

    deltapd_copy = [[] for i in range(num_node)]
    deltaYoe_copy = [[] for i in range(num_node)]
    deltaYue_copy = [[] for i in range(num_node)]
    deltaYuesum_copy = [[] for i in range(sub)]
    TIMEstart=0
    TIMEend=TIMEpiece[0]
    for numt in range(len(TIMEpiece)-1):
        TIMEstart+=TIMEpiece[numt]
        TIMEend+=TIMEpiece[numt+1]
        a = list(history[:, TIMEstart::T])
        for i in range(len(b)):
            for j in range(len(b[0])):
                for n in range(1,TIMEpiece[numt+1]):
                    c = history[:, TIMEstart + n::T]
                    a[i] = list(a[i])
                    a[i].insert(TIMEpiece[numt+1] * j + n, c[i][j])
        result=f_pd(Pd[TIMEstart:TIMEend], Qd[TIMEstart:TIMEend], a, TIMEpiece[numt+1], beta,S_scale[TIMEstart:],0)

        for i in range(num_node):
            deltapd_copy[i].append(result[0][i]/1000) #除以之前乘的倍数
            deltaYoe_copy[i].append(result[1][i]/1000)
            deltaYue_copy[i].append(result[2][i]/1000)
        for i in range(sub):
            deltaYuesum_copy[i].append(result[3][i]/1000)

    deltapd_copy = np.array(deltapd_copy)
    deltaYoe_copy = np.array(deltaYoe_copy)
    deltaYue_copy = np.array(deltaYue_copy)
    deltaYuesum_copy = np.array(deltaYuesum_copy)
    deltapd=deltapd_copy.reshape(num_node,T)
    deltaYoe=deltaYoe_copy.reshape(num_node,T)
    deltaYue=deltaYue_copy.reshape(num_node,T)
    deltaYuesum=deltaYuesum_copy.reshape(sub,T)
    return deltapd,deltaYoe,deltaYue,deltaYuesum

