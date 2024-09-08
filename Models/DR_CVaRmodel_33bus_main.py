import numpy as np
import datetime
from DR_CVaRmodel_power_load_with33bus import funcpd
from DR_CVaRmodel_WT_output_with33bus import funcwt
from DR_CVaRmodel_with_33bus import funcmain
Pd = [3980, 3610, 3340, 2990, 2690, 2920, 3260, 3550, 3980, 4290, 4470, 4610, 4730, 4840, 4930, 4710, 4520, 4830,
      5090, 5385, 4980, 4610, 4410, 4290]  # 假定各时段有功负荷波动
Qd = [0]*24
for i in range(24):
    Qd[i] = Pd[i] *2300/3715 # 按33-bus计算无功相对于有功的值
expectedpd=max(Pd) # 假定的最大负荷
sub = 1
T = 24
ST = sub * T
history = []  # 历史数据
M = 366  # 历史数据的数量
np.random.seed(678) # 随机种子
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

history = np.zeros((length, 24 * 2 * sub))
for i in range(length - 24):
    for j in range(T):
        history[i, j] = lines[j + (i // 8) * 24][8 - 1 - i % 8]
for i in range(length - 24, length):
    for j in range(T):
        history[i, j] = lines[j + (i // 8) * 24][6 - 1 - i % 6]
for i in range(length):
    for j in range(T):
        history[i, j + sub * T] = history_rt[i * T + j]
history=np.array(history)

beta = 0.75
result = funcpd(beta)  # pd波动及预期误差值
result1 = funcwt()  # 风电预期出力
deltapd = result[0]
deltaYoe = result[1]
deltaYue = result[2]
deltaYuesum = result[3]
PWT = result1[0]
PWT=np.array(PWT)
deltapd=np.array(deltapd)
deltaYoe = np.array(deltaYoe)
deltaYuesum = np.array(deltaYuesum)
Yoe = np.array(result1[1])
Yue = np.array(result1[2])
TIMEpiece=[0,7,3,11,3]
b=history[:,1::T]
TIMEstart=0
TIMEend=TIMEpiece[0]
obj=0
start=datetime.datetime.now().timestamp()
for numt in range(len(TIMEpiece)-1):
    TIMEstart+=TIMEpiece[numt]
    TIMEend+=TIMEpiece[numt+1]
    a = list(history[:, TIMEstart::T])
    new_PWT = PWT[:, TIMEstart:TIMEend]
    new_deltapd=deltapd[:,TIMEstart:TIMEend]
    new_Yoe=Yoe[:,TIMEstart:TIMEend]
    new_Yue=Yue[:,TIMEstart:TIMEend]
    new_pdYoe=deltaYoe[:,TIMEstart:TIMEend]
    new_pdYuesum=deltaYuesum[:,TIMEstart:TIMEend]
    for i in range(len(b)):
        for j in range(len(b[0])):
            for n in range(1, TIMEpiece[numt+1]):
                c = history[:, TIMEstart + n::T]
                a[i] = list(a[i])
                a[i].insert(TIMEpiece[numt+1] * j + n, c[i][j])
    x_1=funcmain(Pd[TIMEstart:TIMEend], Qd[TIMEstart:TIMEend], a, TIMEpiece[numt+1], new_PWT,new_deltapd,new_Yoe,new_Yue,new_pdYoe,new_pdYuesum,beta)
    obj+=x_1
print('objVal: ',obj)
end=datetime.datetime.now().timestamp()
print('time: ',end-start)
