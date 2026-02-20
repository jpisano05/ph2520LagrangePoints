import numpy as np

m1 = 0
m2 = 0
R = 0

def calculatePoints(m1, m2, R):
    mu = m2 / (m1 + m2)

    L1 = R * ((mu/3) ** (1/3))
    L2 = -1 * L1

    L3 = R - (R * (7/12 * mu))

    
    return L1, L2, L3
    
if __name__ == "__main__":
    #test for earth/sun system
    m1 = 1.989e30
    m2 = 5.972e24
    R = 1.496e11
    
    L1, L2, L3 = calculatePoints(m1, m2, R)
    
    print(L1)
    print(L2)
    print(L3)