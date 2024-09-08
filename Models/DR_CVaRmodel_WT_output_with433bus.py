from scipy.special import gamma, gammaincc
import math
def funcwt():
    WT_k=[1.572,1.554,1.578,1.648,1.805,1.936,1.991,1.907,1.920,1.736,1.709,1.695,
                      1.695,1.643,1.730,1.621,1.594,1.469,1.454,1.470,1.474,1.493,1.475,1.516] # weibull公式的24小时k参数
    WT_c=[7.225,8.007,8.534,8.954,9.242,9.897,10.047,10.232,9.964,9.241,8.256,7.541,
          7.455,7.651,7.688,7.656,7.347,6.874,6.664,6.589,6.573,6.304,6.270,6.596] # weibull公式的24小时c参数
    MVAbase=10
    T=24
    WT = 6
    WTmax = 800 / (MVAbase * 1e3)
    WTfactor = 0.80
    WT_r = WTmax * WTfactor
    vin,vr,vout = 3,8,16
    PWT1=[0]*24
    result=[]
    for m in range(T):#考虑24时段
        S = []
        MINN = 1E6 # 最小目标值
        PWT_1=0 # 最小目标值对应的PWT
        for n in range(401):
            s = 0
            for h in range(1):
                PWT = PWT1[m:][:]
                PWT[h]=PWT[h]+n*0.00016  # 令PWT变化
                PWT[h]=min(PWT[h],WT_r)
                PWT[h]=max(PWT[h],0)#调整最大最小上限
                t=m
                s+=80*(PWT[h]*(1-math.exp(-(vin/WT_c[t])**WT_k[t] )+math.exp(-(vout/WT_c[t])**WT_k[t] ))+(WT_r*vin/(vr-vin)
                          +PWT[h])*(math.exp(-(vin/WT_c[t])**WT_k[t]) -math.exp(-((vin+(vr-vin)*PWT[h]/WT_r)/WT_c[t])**WT_k[t]))
                         +WT_r*WT_c[t]/(vr-vin)*((gammaincc(1 + 1 / WT_k[t], ((vin+(vr-vin)*PWT[h]/WT_r) / WT_c[t]) ** WT_k[t])
                         * gamma(1 + 1 / WT_k[t]) -gammaincc(1 + 1 / WT_k[t], (vin / WT_c[t]) ** WT_k[t]) *gamma(1 + 1 / WT_k[t])) ))
                s+=100*( (WT_r - PWT[h]) * (math.exp(-vr ** WT_k[t] / (WT_c[t] ** WT_k[t]))
                          - math.exp(-vout ** WT_k[t] /(WT_c[t] ** WT_k[t]))) + (WT_r * vin / (vr - vin) + PWT[h]) * (
                          math.exp(-vr ** WT_k[t] / (WT_c[t] ** WT_k[t])) - math.exp(-((vin + (vr - vin) * PWT[h] / WT_r)
                          / WT_c[t]) ** WT_k[t]))+ WT_r * WT_c[t] / (vr - vin) * (gammaincc(1 + 1 / WT_k[t],
                          ((vin + (vr - vin) * PWT[h] / WT_r) / WT_c[t]) ** WT_k[t]) * gamma(1 + 1 / WT_k[t])
                          - gammaincc(1 + 1 / WT_k[t], (vr / WT_c[t]) ** WT_k[t]) * gamma(1 + 1 / WT_k[t])) )
            #如果当前解小于已知最小目标值，则将其记录
            if MINN>s:
                PWT_1=PWT[0]
                MINN=s
            S.append(s)
        result.append(round(PWT_1,4))
    temp1 = []
    j = 1
    while j <= WT:
        temp1.append(result)
        j += 1
    yoe=[[0 ]*T for i in range(WT)]
    yue = [[0]*T for i in range(WT)]
    for i in range(WT):
        for t in range(T):
            yoe[i][t]=temp1[i][t]* ( 1 - math.exp(-(vin / WT_c[t]) ** WT_k[t]) + math.exp(-(vout / WT_c[t]) ** WT_k[t]))+(
                    WT_r*vin/(vr-vin) +temp1[i][t])*(math.exp(-(vin / WT_c[t]) ** WT_k[t]) -math.exp(
                -((vin + (vr - vin) * temp1[i][t] / WT_r) / WT_c[t]) ** WT_k[t])) +WT_r * WT_c[t] / (vr - vin) * (
                (gammaincc(1 + 1 / WT_k[t], ((vin + (vr - vin) * temp1[i][t] / WT_r) / WT_c[t]) ** WT_k[t])
                 * gamma(1 + 1 / WT_k[t]) - gammaincc(1 + 1 / WT_k[t], (vin / WT_c[t]) ** WT_k[t]) * gamma(1 + 1 / WT_k[t])) )
            yue[i][t]=(WT_r - temp1[i][t]) * (math.exp(-vr ** WT_k[t] / (WT_c[t] ** WT_k[t]))
                          - math.exp(-vout ** WT_k[t] /(WT_c[t] ** WT_k[t]))) + (WT_r * vin / (vr - vin) + temp1[i][t]) * (
                          math.exp(-vr ** WT_k[t] / (WT_c[t] ** WT_k[t])) - math.exp(-((vin + (vr - vin) * temp1[i][t] / WT_r)
                          / WT_c[t]) ** WT_k[t]))+ WT_r * WT_c[t] / (vr - vin) * (gammaincc(1 + 1 / WT_k[t],
                          ((vin + (vr - vin) * temp1[i][t] / WT_r) / WT_c[t]) ** WT_k[t]) * gamma(1 + 1 / WT_k[t])
                          - gammaincc(1 + 1 / WT_k[t], (vr / WT_c[t]) ** WT_k[t]) * gamma(1 + 1 / WT_k[t]))
    return temp1,yoe,yue
