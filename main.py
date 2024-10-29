#!/usr/bin/env -S python

from math import sqrt, pi

from matplotlib import pyplot as plt

from table import *
from measurements import *

def δln2(*xδs): return sum(δx**2/x**2 for x, δx in xδs)

def vδ(dδ, tδ):
    d, t = (xδ[0] for xδ in (dδ, tδ))

    v = d/t

    δv = abs(v)*sqrt(δln2(dδ, tδ))

    return (v, δv)

def ξδ(bδ, pδ):
    b, _ = bδ
    p, _ = pδ

    ξ = b/(2*p)

    δξ = ξ*sqrt(δln2(bδ, pδ))

    return (ξ, δξ)

def a_0δ(ηδ, v_fδ, ρδ, gδ):
    η, _ = ηδ
    v_f, _ = v_fδ
    ρ, _ = ρδ
    g, _ = gδ

    a_0 = sqrt((9*η*v_f)/(2*ρ*g))
    
    δa_0 = (a_0/2)*sqrt(δln2(ηδ, v_fδ, ρδ, gδ))

    return (a_0, δa_0)

def aδ(ξδ, a_0δ):
    ξ, δξ = ξδ
    a_0, δa_0 = a_0δ

    a = sqrt(ξ**2 + a_0**2) - ξ

    δa = sqrt((a_0**2)*(δa_0**2) + (a**2)*(δξ**2))/abs(a + ξ)

    return (a, δa)

def qδ(aδ, ρδ, gδ, dδ, v_fδ, v_rδ, Vδ):
    a, ρ, g, d, V = (xδ[0] for xδ in (aδ, ρδ, gδ, dδ, Vδ))
    v_f, δv_f = v_fδ
    v_r, δv_r = v_rδ

    q = (4*pi/3) * (a**3) * (ρ*g*d*(v_f + v_r))/(V*v_f)

    δq = abs(q) * sqrt(9*δln2(aδ) + δln2(ρδ, gδ, dδ) + (δv_f**2 + δv_r**2)/(v_f**2 + v_r**2) + δln2(Vδ, v_fδ))

    return (q, δq)

if __name__ == "__main__":

    Tδs = [RT_table_lintearpolate(*Rδ) for Rδ in Rδs]
    ηδs = [Tη_table_lintearpolate(*Tδ) for Tδ in Tδs]

    # for ηδ in ηδs: print(ηδ)

    v_fδs = [vδ(*vals) for vals in zip(d_fδs, t_fδs)]
    v_rδs = [vδ(*vals) for vals in zip(d_rδs, t_rδs)]

    ξδs = [ξδ(*vals) for vals in zip(bδs, pδs)]

    a_0δs = [a_0δ(*vals) for vals in zip(ηδs, v_fδs, ρδs, gδs)]

    aδs = [aδ(*vals) for vals in zip(ξδs, a_0δs)]

    qδs = [qδ(*vals) for vals in zip(aδs, ρδs, gδs, dδs, v_fδs, v_rδs, Vδs)]

    for q, δq in qδs:
        print(f"{q} ± {δq}")

    # plt.title("Thermistor Resistance Table Data in Plot Format")
    # plt.xlabel("Thermistor Resistance (MΩ)")
    # plt.ylabel("Air Temperature (°C)")
    # plt.plot(RT_table_Rs, RT_table_Ts, 'b-', label="Linear Interpolation")
    # plt.plot(RT_table_Rs, RT_table_Ts, 'k.', label="Table Entries")
    # plt.legend()
    # plt.show()