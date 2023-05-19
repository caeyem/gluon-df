import os
import sys
import time

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta as spbeta

def red(s):
    return '\033[1;31m%s\033[m' % s

def log(*m):
    print(" ".join(map(str, m)))

def log_exit(*m):
    log(red("ERROR:"), *m)
    exit(1)

def generatePI():
    # charts
    yList = np.array([0.6275]) #, 0.7 , 0.707, 0.9, 0.1])
    NList = np.array([1e7]) #, 3.25e6, 1e7, 1e7, 1e7])
    binsList = np.array([1e2]) #, 1e2, 1e2, 1e2, 1e2])

    # %
    # mean field theory of spin glasses
    # ---------------------------------------------------------------------------------
    # approach 1: using beta distribution

    # distribution check

    # try using beta distribution
    xx = np.arange(0, 1.01, 0.01)
    for idx in range(0, yList.size):
        y = yList[idx]
        for nn in range(1, 10, 1):
            fig1, ax = plt.subplots()
            fig1.suptitle(" ".join(['Sample Beta distributions: y: ', str(y), ' n: ', str(nn)]))
            x = np.linspace(spbeta.ppf(0.01, y, (1-y)*nn),
                            spbeta.ppf(0.99, y, (1-y)*nn), 100)
            ax.plot(x, spbeta.pdf(x, y, (1-y)*nn),
                            'r-', lw=5, alpha=0.6, label='beta pdf')
            fig1.show()

    for idx in range(0, yList.size):
        nn = 1
        y = float(yList[idx])
        N = int(NList[idx])
        bins = int(binsList[idx])

        xn = np.arange(0, 1 + 1/N, 1/N)
        na = np.arange(1, xn.size+1, 1)

        alpha = y
        beta1 = na * (1 - y)
        beta2 = nn * (1 - y)
        rhonx1 = np.random.beta(alpha, beta1)
        rhonx2 = np.random.beta(alpha, beta2, (1, N))

        # probabilities
        W = rhonx2
        Wmax = np.zeros((1, N + 1))
        Wcmax = np.zeros((1, N + 1))
        Y = np.zeros((1, N + 1))
        Wfac = 1
        for n in range(0, N):
            Wmax[0,n+1] = np.maximum((1 - W[0,n]) * Wmax[0,n], W[0,n])
            Wcmax[0,n+1] = np.maximum(np.minimum((1 - W[0,n]) * Wmax[0,n], W[0,n]), (1 - W[0,n]) * Wcmax[0,n])
            Y[0,n+1] = np.power(W[0,n], 2) + (np.power((1 - W[0,n]), 2)) * Y[0,n]

        fig2, (ax1, ax2, ax3) = plt.subplots(1,3)
        (cnt1, ctr1) = np.histogram(Wmax, bins)
        ax1.plot(ctr1[0:ctr1.size-1], bins * cnt1 / N)
        ax1.set_xlim(0, 1)
        ax1.set_ylim(0, 4)
        ax1.set_xlabel('W: W_{max}')
        ax1.set_ylabel('P_1(W)')
        ax1.set_title('Max Valley Weight Distribution')
        (cnt2, ctr2) = np.histogram(Wcmax, bins)
        ax2.plot(ctr2[0:ctr2.size - 1], bins * cnt2 / N)
        ax2.set_xlim(0, 1)
        ax2.set_ylim(0, 4)
        ax2.set_xlabel('W: W_{max}^c')
        ax2.set_ylabel('P_2(W)')
        ax2.set_title('Second Max Valley Weight Distribution')
        (cnt3, ctr3) = np.histogram(Y, bins)
        ax3.plot(ctr3[0:ctr3.size - 1], bins * cnt3 / N)
        ax3.set_xlim(0, 1)
        ax3.set_ylim(0, 4)
        ax3.set_xlabel('Y')
        ax3.set_ylabel('\Pi(Y)')
        ax3.set_title('Overlap Distribution')
        fig2.suptitle(''.join(['Distributions for the SG/GluonTM model ', ' y:', str(y), ' N:', str(N), ' bins:', str(bins)]))
        fig2.show()

    return


def main():
    generatePI()
    try:
        generatePI()
    except Exception:
        log_exit("Couldn't run")
    a = input('Press to continue')

if __name__ == "__main__":
    main()