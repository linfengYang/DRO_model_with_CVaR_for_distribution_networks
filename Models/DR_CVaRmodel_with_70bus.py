import gurobipy as gp
import numpy as np
from gurobipy import GRB

def funcmain(Pd, Qd, history, time, PWT, deltapd, Yoe, Yue, pdYoe, pdYuesum, beta):
    with open('bus70.txt', 'r') as file:
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
    with open('branch70.txt', 'r') as file:
        f = file.readline()
        branch = []
        while f:
            a = []
            b = f.split()
            for i in b[:-1]:
                a.append(float(i))
            branch.append(a)
            f = file.readline()
    branch = np.array(branch)
    Vbase = 11  # 基准电压
    MVAbase = 1  # 基准视在功率
    Ibase = 90.90  # 基准电流
    branch[:, 2] = branch[:, 2] / ((Vbase) ** 2 / MVAbase)
    branch[:, 3] = branch[:, 3] / ((Vbase) ** 2 / MVAbase)

    WT = 2  # 风电台数
    WT_bus = [13, 58]  # 风电位置
    WTmax = 600 / (MVAbase * 1e3)  # 额定输出
    WTfactor = 0.80
    WT_r = WTmax * WTfactor  # 额定有功

    DG = 2  # DG数量
    DG_bus = [26, 62]  # DG位置
    DGfactor = 0.80
    DGPmax = 1400 / (MVAbase * 1e3)
    DGPmin = 100 / (MVAbase * 1e3)
    DGcost = [[0.575, 159.325, 0], [0.625, 158.875, 0]]

    SVC = 3  # 无功补偿器数量
    SVCmin = 0
    SVCmax = 800  # KVar
    SVCmax /= (MVAbase * 1e3)
    SVC_bus = [22, 47, 50]  # 补偿器位置

    # 支路容量
    Plinemax = 5000 / (MVAbase * 1e3)  # KW
    Qlinemax = 4000 / (MVAbase * 1e3)  # KVar

    vin = 3  # 切入风速
    vr = 8  # 额定风速
    vout = 16  # 切出风速

    T = time
    sub = 2
    sub_bus = [1, 70]
    ST = sub * T
    num_branch = len(branch)  # 分支数量
    num_node = len(bus)  # 节点数量
    MAX = 10
    for i in range(T):
        Pd[i] *= 1.5
        Qd[i] *= 1.5
        Pd[i] /= (MVAbase * 1e3)
        Qd[i] /= (MVAbase * 1e3)
    sp = sum(bus[:, 2])  # 总负荷
    sq = sum(bus[:, 3])
    pd_factor = [i[2] / sp for i in bus]  # 获得每个节点的负荷因子，即负荷占总负荷的比例,
    qd_factor = [i[3] / sq for i in bus]

    pd = [[i * j for i in Pd] for j in pd_factor]  # 获得每个节点各个时段的有功负荷,
    qd = [[i * j for i in Qd] for j in qd_factor]
    pd = np.array(pd)
    qd = np.array(qd)

    Vmax = 1.1
    Vmin = 0.9
    np.random.seed(678)
    IL = 6
    ILbus = [4, 21, 28, 39, 45, 66]  # 可中断负荷bus
    ILrate = 0.3  # 可中断负荷占比
    totalIL = [0] * T  # 总可中断负荷
    for t in range(T):
        temp = 0
        for i in range(len(ILbus)):
            temp += 0.25 * pd[ILbus[i] - 1, t]
        totalIL[t] = temp
    IL_alpha = [[0 for t in range(T)] for i in range(IL)]
    IL_beta = [[0 for t in range(T)] for i in range(IL)]
    for i in range(IL):
        for t in range(T):
            somean = np.random.normal(0, 0.005)
            soscale = np.random.normal(0, 0.002)
            IL_alpha[i][t] = 0.625 * (1 - Pd[t] / (5.385 * 1.5) * 0.1) + np.random.normal(0 + somean * 2,0.02 + soscale * 1.5)
            IL_beta[i][t] = 150 * (1 - Pd[t] / (5.385 * 1.5) * 0.1) + np.random.normal(0 + somean * 1e3, 3 + soscale * 2e2)

    m = gp.Model()
    x = m.addVars(num_branch, vtype=GRB.BINARY, name='x')  # 开关01变量
    DG_x = m.addVars(DG, T, vtype=GRB.BINARY, name='DG_x')  # 表示DG是否开机
    aij = m.addVars(num_branch, vtype=GRB.CONTINUOUS, name='a')
    aji = m.addVars(num_branch, vtype=GRB.CONTINUOUS, name='b')
    V2 = m.addVars(num_node, T, lb=Vmin ** 2, ub=Vmax ** 2, vtype=GRB.CONTINUOUS, name='V2')  # 电压标幺值
    I2 = m.addVars(num_branch, T, vtype=GRB.CONTINUOUS, name='I2')  # 电压标幺值
    Pnode = m.addVars(num_node, T, lb=-Plinemax, vtype=GRB.CONTINUOUS, name='Pnode')  # 节点产生的绝对有功
    Qnode = m.addVars(num_node, T, lb=-Qlinemax, vtype=GRB.CONTINUOUS, name='Qnode')
    Plinehat = m.addVars(num_branch, T, lb=-Vmax * Plinemax, vtype=GRB.CONTINUOUS, name='Pline')  # 线路上的绝对有功
    Qlinehat = m.addVars(num_branch, T, lb=-Vmax * Qlinemax, vtype=GRB.CONTINUOUS, name='Qline')
    PG = m.addVars(num_node, T, vtype=GRB.CONTINUOUS, name='PG')  # 每个节点的有功发电
    QG = m.addVars(num_node, T, vtype=GRB.CONTINUOUS, name='QG')
    lp = m.addVars(DG, T, vtype=GRB.CONTINUOUS, name='lp')  # DG发电成本的线性近似
    alpha = m.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="alpha")  # 在β概率下，保证损失小于α的最小的α
    subP = m.addVars(2 * sub, T, vtype=GRB.CONTINUOUS, name='subP')  # 变电站DA与RT出力

    PIL = m.addVars(IL, T, vtype=GRB.CONTINUOUS, name='PIL')  # 可中断负荷的中断值
    QIL = m.addVars(IL, T, vtype=GRB.CONTINUOUS, name='QIL')
    z_IL = m.addVars(IL, T, vtype=GRB.CONTINUOUS, name='z_IL')  # 可中断负荷成本
    QM = m.addVars(num_node, T, vtype=GRB.CONTINUOUS, name='QM')  # 不考虑负荷波动的节点无功
    PM = m.addVars(num_node, T, vtype=GRB.CONTINUOUS, name='PM')

    lpL = 5  # 线性发电量分段
    DGproduct_l = [[l * (DGPmax) / lpL for l in range(lpL + 1)] for i in range(DG)]  # 发电量分段
    numIL_L = 5  # 可中断负荷分段
    IL_L = [[[l * (ILrate * pd[ILbus[i] - 1, t]) / numIL_L for l in range(numIL_L + 1)] for t in range(T)]
            for i in range(IL)]
    IL_L = np.array(IL_L)

    # IL约束
    m.addConstrs(z_IL[i, t] >= IL_alpha[i][t] * (100 * IL_L[i, t, l] ** 2 + 2 * 100 * IL_L[i, t, l] * (PIL[i, t] - IL_L[i, t, l])) +
                 100 * IL_beta[i][t] * PIL[i, t] for i in range(IL) for t in range(T) for l in range(numIL_L))
    m.addConstrs(PIL[i, t] <= ILrate * pd[ILbus[i] - 1, t] for i in range(IL) for t in range(T))
    m.addConstrs(gp.quicksum(PIL[i, t] for i in range(IL)) <= totalIL[t] for t in range(T))
    m.addConstrs(PIL[i, t] * qd[i, t] == QIL[i, t] * pd[i, t] for i in range(IL) for t in range(T))

    # 线路开关约束
    m.addConstrs((Plinehat[i, t] <= x[i] * MAX for i in range(num_branch) for t in range(T)), name='c1')
    m.addConstrs((Plinehat[i, t] >= -x[i] * MAX for i in range(num_branch) for t in range(T)), name='c1_1')
    m.addConstrs((Qlinehat[i, t] <= x[i] * MAX for i in range(num_branch) for t in range(T)), name='c2')
    m.addConstrs((Qlinehat[i, t] >= -x[i] * MAX for i in range(num_branch) for t in range(T)), name='c2_1')

    # 线路无孤岛
    m.addConstrs(gp.quicksum(aji[j] for j in range(num_branch) if branch[j][1] == i) == 0 for i in sub_bus)
    m.addConstrs(gp.quicksum(aij[j] for j in range(num_branch) if branch[j][0] == i) == 0 for i in sub_bus)
    m.addConstrs(gp.quicksum(aji[j] for j in range(num_branch) if branch[j][1] == i + 1) +
                 gp.quicksum(aij[j] for j in range(num_branch) if branch[j][0] == i + 1) == 1 for i in
                 range(num_node) if (i + 1) not in sub_bus)
    m.addConstrs(aij[i] + aji[i] == x[i] for i in range(num_branch))

    # 二阶锥约束
    m.addConstrs(I2[i, t] <= x[i] * 6 ** 2 for i in range(num_branch) for t in range(T))
    m.addConstrs(
        I2[i, t] * V2[branch[i, 0] - 1, t] >= Plinehat[i, t] ** 2 + Qlinehat[i, t] ** 2 for i in range(num_branch) for t
        in range(T))

    # 节点输入输出平衡
    m.addConstrs((gp.quicksum(Plinehat[k, t] for k in range(num_branch) if branch[k][0] == i + 1) ==
                  gp.quicksum(Plinehat[j, t] - branch[j, 2] * I2[j, t] for j in range(num_branch) if branch[j][1] == i + 1)
                  + Pnode[i, t] for i in range(num_node) for t in range(T)), name='c3')
    m.addConstrs((gp.quicksum(Qlinehat[k, t] for k in range(num_branch) if branch[k][1] == i + 1) ==
                  gp.quicksum( Qlinehat[j, t] - branch[j, 2] * I2[j, t] for j in range(num_branch) if branch[j][0] == i + 1)
                  + Qnode[i, t] for i in range(num_node) for t in range(T)), name='c3')
    m.addConstrs((Plinehat[i, t] <= Plinemax for i in range(num_branch) for t in range(T)), name='c23')
    m.addConstrs(-Plinemax <= Plinehat[i, t] for i in range(num_branch) for t in range(T))
    m.addConstrs((Qlinehat[i, t] <= Qlinemax for i in range(num_branch) for t in range(T)), name='c24')
    m.addConstrs(Qlinehat[i, t] >= -Qlinemax for i in range(num_branch) for t in range(T))

    deltapd_pd = []
    deltaqd_qd = []
    for i in range(num_node):
        for t in range(T):
            if pd[i, t]:deltapd_pd.append(deltapd[i, t] / pd[i, t])
            else: deltapd_pd.append(0)
            if qd[i, t]:deltaqd_qd.append(0.25 * deltapd[i, t] / qd[i, t])
            else:deltaqd_qd.append(0)

    # 节点产生有功无功约束
    m.addConstrs(
        (Pnode[i, t] == PM[i, t] - (1 + deltapd_pd[i * T + t]) * pd[i, t] for i in range(num_node)
         for t in range(T)), name='c6')
    m.addConstrs(
        (Qnode[i, t] == QM[i, t] - (1 + deltaqd_qd[i * T + t]) * qd[i, t] for i in range(num_node)
         for t in range(T)), name='c6_q')
    m.addConstrs(PM[ILbus[i] - 1, t] == PIL[i, t] for i in range(IL) for t in range(T))
    m.addConstrs(QM[ILbus[i] - 1, t] == QIL[i, t] for i in range(IL) for t in range(T))
    m.addConstrs(PM[i, t] == 0 for i in range(num_node) if (i + 1) not in DG_bus and (i + 1) not in WT_bus
                 and (i + 1) not in sub_bus and (i + 1) not in ILbus for t in range(T))
    m.addConstrs(QM[i, t] == 0 for i in range(num_node) if (i + 1) not in DG_bus and (i + 1) not in WT_bus
                 and (i + 1) not in sub_bus and (i + 1) not in SVC_bus and (i + 1) not in ILbus for t in range(T))
    m.addConstrs((PM[(WT_bus[i] - 1), t] == PWT[i, t] for i in range(WT) for t in range(T)), name='c10')
    m.addConstrs((QM[(WT_bus[i] - 1), t] == 0.25 * PWT[i, t] for i in range(WT) for t in range(T)), name='c11')
    m.addConstrs((PM[(sub_bus[i] - 1), t] == PG[(sub_bus[i] - 1), t] for i in range(sub) for t in range(T)), name='c10')
    m.addConstrs((PG[(sub_bus[i] - 1), t] == subP[i, t] + subP[sub + i, t] for i in range(sub) for t in range(T)),name='c10')
    m.addConstrs((QM[(sub_bus[i] - 1), t] == QG[(sub_bus[i] - 1), t] for i in range(sub) for t in range(T)), name='c10')
    m.addConstrs(QM[(SVC_bus[i] - 1), t] <= SVCmax for i in range(SVC) for t in range(T))
    m.addConstrs((PM[(DG_bus[i] - 1), t] <= DGPmax * DG_x[i, t] for i in range(DG) for t in range(T)), name='c10')
    m.addConstrs((PM[(DG_bus[i] - 1), t] >= DGPmin * DG_x[i, t] for i in range(DG) for t in range(T)), name='c10')
    m.addConstrs((QM[(DG_bus[i] - 1), t] == 0.25 * PM[(DG_bus[i] - 1), t] for i in range(DG) for t in range(T)), name='c11')

    # 电压约束
    m.addConstrs((V2[(sub_bus[i] - 1), t] == 1 for i in range(sub) for t in range(T)), name='c28')
    m.addConstrs(V2[branch[i, 0] - 1, t] - V2[branch[i, 1] - 1, t] <= (1 - x[i]) * MAX + 2 * (
                    branch[i, 2] * Plinehat[i, t] + branch[i, 3] * Qlinehat[i, t])
                 - (branch[i, 2] ** 2 + branch[i, 3] ** 2) * I2[i, t] for i in range(num_branch) for t in range(T))
    m.addConstrs(V2[branch[i, 0] - 1, t] - V2[branch[i, 1] - 1, t] >= -(1 - x[i]) * MAX + 2 * (
                    branch[i, 2] * Plinehat[i, t] + branch[i, 3] * Qlinehat[i, t])
                 - (branch[i, 2] ** 2 + branch[i, 3] ** 2) * I2[i, t] for i in range(num_branch) for t in range(T))
    m.addConstrs((((2 * DGcost[i][0] * DGproduct_l[i][l] + DGcost[i][1]) * PM[(DG_bus[i] - 1), t] + DGcost[i][2] -
                   DGcost[i][0] * (DGproduct_l[i][l] ** 2)) <= lp[i, t] for i in range(DG) for t in range(T)
                  for l in range(lpL + 1)), name="c11")

    history = np.array(history)
    ST *= 2
    M, L1 = 366, 30
    x_low = np.array([min(history[:, i]) for i in range(ST)])
    x_high = np.array([max(history[:, i]) for i in range(ST)])
    s, C_ = [x_low, x_high], []
    for i in range(ST):
        C_.extend([x_low[i] + (x_high[i] - x_low[i]) * l / (L1) for l in range(1, L1 + 1)])
    C_ = np.array(C_)
    gama = [sum([max(history[j][i] - C_[i * (L1) + l], 0) for j in range(M)]) / M + 1 for i in range(ST)
            for l in range(L1)]

    # 模糊集矩阵化
    C1 = np.diag([1] * ST)
    C2 = np.diag([-1] * ST)
    C3 = np.diag([0] * ST)
    C4 = [C1, C2]
    [C4.append(C3) for _ in range(2 * (L1))]
    C = np.concatenate(C4, axis=0)
    k = 0
    for i in range((2 + (L1)) * ST, (2 + 2 * (L1)) * ST):
        C[i, k] = 1
        if not (i - (2 + (L1)) * ST + 1) % (L1):
            k += 1

    D1 = np.diag([-1] * (L1) * ST)
    D = np.vstack((np.zeros((2 * ST, (L1) * ST)), D1, D1))
    d = np.zeros(((L1) * ST))
    d = np.hstack((d, np.zeros((2 * ST)), d))
    d[:ST] = x_high
    d[ST:2 * ST] = -x_low
    d[(2 + (L1)) * ST:] = C_
    G = 2 * (L1) * ST + 2 * ST

    # 对偶化模型
    eta = m.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="eta")
    omega = m.addVars((L1) * ST, vtype=GRB.CONTINUOUS, name="omega")
    epi = m.addVars(ST, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="epi")
    mean = np.mean(history, axis=0)

    pi1 = m.addVars(G, vtype=GRB.CONTINUOUS, name="pi1")
    m.addConstr(-gp.quicksum(d[i] * pi1[i] for i in range(G)) + eta - alpha >= 0, name="c12")
    m.addConstrs((gp.quicksum(C[j, i] * pi1[j] for j in range(G)) == -epi[i] for i in range(ST)), name="c13")
    m.addConstrs((gp.quicksum(D[j, i] * pi1[j] for j in range(G)) == -omega[i] for i in range((L1) * ST)), name="c14")

    pi2 = m.addVars(G, vtype=GRB.CONTINUOUS, name="pi2")
    m.addConstr(-gp.quicksum(d[i] * pi2[i] for i in range(G)) + eta + (beta * alpha
                - 0.01 * gp.quicksum(z_IL[i, t] for i in range(IL) for t in range(T))
                - gp.quicksum(lp[i, t] for i in range(DG) for t in range(T)) - 100 * gp.quicksum(
                branch[i, 2] * (I2[i, t]) for i in range(num_branch) for t in range(T))
                - 80 * sum([Yoe[i, t] for i in range(WT) for t in range(T)]) - 100 * sum(
                [Yue[i, t] for i in range(WT) for t in range(T)]) - 100 * sum(
                [pdYoe[i, t] for i in range(num_node) for t in range(T)])) / (1 - beta) >= 0, name="c16")
    m.addConstrs((gp.quicksum(C[j, i] * pi2[j] for j in range(G)) ==
                  subP[i // T, i % T] / (1 - beta) - epi[i] for i in range(ST // 2)), name="c17")
    m.addConstrs((gp.quicksum(C[j, i] * pi2[j] for j in range(G)) ==
                  (subP[i // T, i % T] + pdYuesum[(i - ST // 2) // T, (i - ST // 2) % T]) / (1 - beta) - epi[i]
                  for i in range(ST // 2, ST)), name="c17")
    m.addConstrs((gp.quicksum(D[j, i] * pi2[j] for j in range(G)) == -omega[i] for i in range((L1) * ST)), name="c18")

    m.Params.MIPGAP = 0
    m.Params.OutputFlag = 0
    m.setObjective(eta + gp.quicksum(omega[i] * gama[i] for i in range((L1) * ST)) + gp.quicksum(epi[i] * mean[i] for i in range(ST)))
    m.update()
    m.optimize()
    a = m.objVal
    m.dispose()
    return a