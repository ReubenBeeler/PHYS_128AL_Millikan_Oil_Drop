from math import sqrt

def find_i(val, vals):
    if not isinstance(vals, list): vals = list(vals)
    for i in range(len(vals) - 1):
        if (val - vals[i] > 0) != (val - vals[i+1] > 0):
            return i
    raise Exception("can't find i!")

def lintearpolate(x, δx, xs, δxs, ys, δys):
    i = find_i(x, xs)
    x_i, x_f = xs[i:i+2]
    δx_i, δx_f = δxs[i:i+2]
    y_i, y_f = ys[i:i+2]
    δy_i, δy_f = δys[i:i+2]
    
    m = (y_f - y_i)/(x_f - x_i)
    ϵ = m*(x - x_i)
    y = ϵ + y_i

    δϵ1 = (δy_f**2 + δy_i**2)/(y_f - y_i)**2
    δϵ2 = (δx_f**2 + δx_i**2)/(x_f - x_i)**2
    δϵ3 = (δx**2 + δx_i**2)/(x - x_i)**2
    δϵ = abs(ϵ) * sqrt(δϵ1 + δϵ2 + δϵ3)
    δy = sqrt(δϵ**2 + δy_i**2)

    return (y, δy)

RT_table_Rs  = [3.239e6 , 3.118e6 , 3.004e6 , 2.897e6 , 2.795e6 , 2.700e6 , 2.610e6 , 2.526e6 , 2.446e6 , 2.371e6 , 2.300e6 , 2.233e6 , 2.169e6 , 2.110e6 , 2.053e6 , 2.000e6 , 1.950e6 , 1.902e6 , 1.857e6 , 1.815e6 , 1.774e6 , 1.736e6 , 1.700e6 , 1.666e6 , 1.634e6 , 1.603e6 , 1.574e6 , 1.547e6 , 1.521e6 , 1.496e6 ]
RT_table_δRs = [0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6, 0.0005e6]
RT_table_Ts  = [  10  ,   11  ,   12  ,   13  ,   14  ,   15  ,   16  ,   17  ,   18  ,   19  ,   20  ,   21  ,   22  ,   23  ,   24  ,   25  ,   26  ,   27  ,   28  ,   29  ,   30  ,   31  ,   32  ,   33  ,   34  ,   35  ,   36  ,   37  ,   38  ,   39  ]
RT_table_δTs = [    .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1,     .1]
def RT_table_lintearpolate(R, δR):
    return lintearpolate(R, δR, RT_table_Rs, RT_table_δRs, RT_table_Ts, RT_table_δTs)

Tη_table_Ts  = [     15  ,     32  ]
Tη_table_δTs = [       .1,       .1]
Tη_table_ηs  = [1.8000e-5, 1.8810e-5]
Tη_table_δηs = [ .0004e-5,  .0004e-5]
def Tη_table_lintearpolate(T, δT):
    return lintearpolate(T, δT, Tη_table_Ts, Tη_table_δTs, Tη_table_ηs, Tη_table_δηs)