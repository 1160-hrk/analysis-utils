#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  5 10:02:50 2023

@author: hirokitsusaka
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 10:51:49 2022

@author: hirokitsusaka
"""
import math as mt
import numpy as np

    
def find_transition(f,J_max,v_max): # PPスペクトル上のピーク信号がどの準位からどの準位への遷移かを判定
    V = np.zeros(12)
    A = [0,0,0,0,1]
    V[0]=2349.1
    V[1]=2324.2
    V[2]=2299.3
    V[3]=2274.4
    V[4]=2249.5
    V[5]=2224.7
    V[6]=2200.0
    V[7]=2175.3
    V[8]=2150.7
    V[9]=2126.1
    V[10]=2101.7
    V[11]=2077.4
    # cv = ['maroon','orangered','peru','goldenrod','yellow',
    #          'greenyellow','mediumseagreen','mediumturquoise',
    #          'mediumblue','darkorchid','fuchsia','crimson']
    cv = ['brown','orangered','orange','mediumseagreen',
             'mediumblue','darkorchid','fuchsia','crimson',
             'brown','orangered','orange','mediumseagreen']
    B = 0.39022
    dB = 0.003076
    # J_max = 25
    for v in range(0,v_max):
        for j in range(0,J_max):
            if v%2 == 0:
                J = j*2
            else:
                J = j*2+1
            Jl = (B-dB*(v))*J*(J+1)
            Ju = (B-dB*(v+1))*(J)*(J-1)
            dJ = -1
            for i in range(0,2):
                diff = abs(V[v]+Ju-Jl-f)
                if diff < A[4]:
                    A[0] = V[v]+Ju-Jl
                    A[1] = v
                    A[2] = J
                    A[3] = dJ
                    A[4] = diff
                Ju = (B-dB*(v+1))*(J+1)*(J+2)
                dJ = 1
    return cv[A[1]],A[1],A[0],A[2],A[3]
def find_transition_w_range(f,w,J_max,v_max): # PPスペクトル上のピーク信号がどの準位からどの準位への遷移かを判定
    V = np.zeros(12)
    A = [0,0,0,0,w]
    C = []
    V[0]=2349.1
    V[1]=2324.2
    V[2]=2299.3
    V[3]=2274.4
    V[4]=2249.5
    V[5]=2224.7
    V[6]=2200.0
    V[7]=2175.3
    V[8]=2150.7
    V[9]=2126.1
    V[10]=2101.7
    V[11]=2077.4
    B = 0.39022
    dB = 0.003076
    # J_max = 25
    for v in range(0,v_max):
        for j in range(0,J_max):
            if v%2 == 0:
                J = j*2
            else:
                J = j*2+1
            Jl = (B-dB*(v))*J*(J+1)
            Ju = (B-dB*(v+1))*(J)*(J-1)
            dJ = -1
            for i in range(0,2):
                diff = abs(V[v]+Ju-Jl-f)
                if diff < A[4]:
                    A[0] = V[v]+Ju-Jl
                    A[1] = v
                    A[2] = J
                    A[3] = dJ
                    A[4] = diff
                    C.append(A)
                    A = [0,0,0,0,w]
                Ju = (B-dB*(v+1))*(J+1)*(J+2)
                dJ = 1
    return C
def transition_wavenumber(v_l,j_l,dj):
    V = np.zeros(12)
    A = [0,0,0,0,1]
    V[0]=2349.1
    V[1]=2324.2
    V[2]=2299.3
    V[3]=2274.4
    V[4]=2249.5
    V[5]=2224.7
    V[6]=2200.0
    V[7]=2175.3
    V[8]=2150.7
    V[9]=2126.1
    V[10]=2101.7
    V[11]=2077.4
    B = 0.39022
    dB = 0.003076
    A[0] = V[v_l] + (B-dB*(v_l+1))*(j_l+dj)*(j_l+dj+1) - (B-dB*v_l)*j_l*(j_l+1)
    return A[0]


def wl_correct_iHR320(wl,angle_cor = 0.027 * mt.pi/180,N=300*10**(-6),theta=21.26*mt.pi/180):
    angle_grating = np.arcsin(wl*N/2/mt.cos(theta/2))
    angle_grating += angle_cor
    wl_cor = np.sin(angle_grating)/N*2*mt.cos(theta/2)
    return wl_cor
def wl_array_board(lambda_ref, ind_pixel, f, theta = 21.56*mt.pi/180, N = 300*10**(-6), dx = 50*10**3, num_pixel = 256, rev = True):
    if rev:
        index = ind_pixel + np.arange(0,-num_pixel,-1)
    else:
        index = np.arange(0,num_pixel) - ind_pixel
    phi_array = np.arctan(index*dx/f)
    phi = np.arcsin(N*lambda_ref/(2*np.cos(theta/2)))
    lambda_array = 2/N*np.cos(theta/2+phi_array/2)*np.sin(phi+phi_array/2)
    if rev:
        lambda_array = np.flip(lambda_array)
    return lambda_array
