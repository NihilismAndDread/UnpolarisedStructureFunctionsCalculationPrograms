import math
from scipy.integrate import quad
from scipy.optimize import brentq
def cutoff(cqm):
    mpi=135
    fpi=93
    Nc=3
    def f(A):
        pr = -1/(16*math.pi**2)
        intnum=3000
        y=[]
        dy=[]
        def y(x):
            return A**2/(1 - x*(1-x)*(mpi/cqm)**2 + A**2) - math.log((1 - x*(1-x)*(mpi/cqm)**2 + A**2)/(1 - x*(1-x)*(mpi/cqm)**2))
        def dy(x):
            return x*(1-x)*A**2/(1-x*(1-x)*(mpi/cqm)**2+A**2)**2 + x*(1-x)/(1-x*(1-x)*(mpi/cqm)**2+A**2) - x*(1-x)/(1-x*(1-x)*(mpi/cqm)**2)
        PI=pr*quad(y, 0, 1)[0]
        dPI=pr*quad(dy,0,1)[0]
        #print(PI, (mpi/cqm)**2*dPI)
        return (fpi/cqm)**2-4*Nc*PI**2/(PI+(mpi/cqm)**2*dPI)
    eps=0.1e-10
    A1=0.00001
    A2=5.0
    i=0

    while (i<50):
        i+=1
        a=f(A1)
        b=f(A2)
        #print("{:e} {:e}".format(a, b))
        if (a*b<0):
            A3=(A1+A2)/2.0
            c = f(A3)
        else:
            print("Should have opposite signs")
            exit()
        if (a*c<0):
            A2=A3
        elif(b*c<0):
            A1=A3
        if (abs(c)<eps):
            #print(A3)
            return A3
            break
    print("End of loop reached, may not have attained descired accuracy.")
    return A3
print(cutoff(400))
for i in [250, 300, 350, 400, 450, 500, 550, 600, 650, 700]:
    print([i,cutoff(i)])