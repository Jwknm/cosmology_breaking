import random
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.stats import poisson
import sympy as sym
from scipy import stats
import scipy


# given fagn
fagn = (0.8, 0.9 ,1.0)  # 임의지정

# set of Ngw
Ngw = np.arange(0, 140)

# Ngw i당 Vi 배당
V = 10 ** 8  # Mpc^3 #레퍼런스 참고해 임의지정
rho_agn = 10 ** (-4.75)  # Mpc^-3 #fiducial AGN density

# i당 source에 맞게 해당 푸아송 분포 따르는 Nagn 배당
lambi = rho_agn * V  # V 고정이므로 lamb도 고정


for t in range(50): #monte carlo 시행 횟수
 Ngw3sigx = []
 Ngw3sigy = []
 for f in fagn: #fagn 변경
     for N in Ngw: #Ngw 변경
         L_Nagn_pd = []
         L0_Nagn_pd = []

         for i in range(N):
             # Si
             random_Nagn_from_sig = np.random.poisson(lambi, 1)
             pd1 = poisson(lambi).pmf(random_Nagn_from_sig - 1)
             # Bi
             random_Nagn_from_bg = np.random.poisson(lambi, 1)
             pd2 = poisson(lambi).pmf(random_Nagn_from_bg)
             # Likelihood
             L_Nagn_pd.append(0.9 * f * pd1 + (1 - 0.9 * f) * pd2)

         for i in range(N):
             random_Nagn_from_L0 = np.random.poisson(lambi, 1)
             pd3 = poisson(lambi).pmf(random_Nagn_from_L0)
             L0_Nagn_pd.append(pd3)

         # lambda 구하기
         L = np.prod(L_Nagn_pd)
         L0 = np.prod(L0_Nagn_pd)

         lamb = 2 * np.log10(L / L0)

         # p-value 구하기
         p_val = scipy.stats.chi2.sf(lamb, 2) #자유도는 뭘까?????
         if p_val < 0.00135: #Ngw의 median의 minimum을 구해야 되는 거????(이건 지금 fagn마다 최소 Ngw만 뽑는 중...)
             Ngw3sigx.append(f)
             Ngw3sigy.append(N)
             break

 print(Ngw3sigx, Ngw3sigy)
 plt.scatter(Ngw3sigx, Ngw3sigy)

plt.show()
