import numpy as np
import matplotlib.pyplot as plt

# grid set up
n           = 100 # number of nodes
L           = 1 # side length of square block
dx          = L/n # width of each node
Xx           = np.linspace(dx/2, L-dx/2,n) # setup of grid X vector
X            = np.ones((1, n)) # setup of grid Y vector
X[0][:]      = Xx*X
#print(Xx,X)
#XX,YY        = np.meshgrid(X, Y) # setup of complete mesh grid
#print(XX,YY)

# time steps
totaltime   = 200 # total time 
timestep    = .5 # time increment for each step of the solver
# note: time steps larger than 0.5 cause the solution to diverge

# physical constants
alpha       = .0001 # thermal diffusivity

# boundary conditions
Tleft       = 100 # left BC
Tright      = 100 # right BC
Tinit       = 0 # initial temperature of block
T           = Tinit*np.ones(np.shape(X)) # the temperature of each node, initially

# plot style
style       = "pcolor" # "pcolor" or "line"

def compute_dTdt(i):
    """ 
    calculates the time derivative of the temperature at the i'th node on the grid X
    i: the i'th node of X

    """
    return alpha*(-(T[0][i]-T[0][i-1])/dx**2 + (T[0][i+1]-T[0][i])/dx**2)

def plot_solution(style="pcolor"):
    # clear the plots 
    plt.ion()
    plt.clf()
    plt.figure(1)

    if style=="pcolor":
        # color grid plot animation
        plt.pcolor(T)

    elif style=="line":
        # regular line plot animation
        plt.plot(X[0][:],T[0][:])
        plt.axis([0, L, Tleft, Tright])

    plt.show()
    plt.pause(0.0001)
    

# forward euler method for T with time step dt

# define temperature gradient at each node
dTdt = np.empty(np.shape(X))

# at each time step, we compute the temperature gradient
for ti in np.arange(0, totaltime, timestep):

    # compute temperature gradient at left boundary node
    dTdt[0][0] = alpha*(-(T[0][0]-Tleft)/dx**2 + (T[0][1]-T[0][0])/dx**2)  
    # compute temperature gradient at inner nodes
    for i in range(1,len(X[0][:])-1):
        dTdt[0][i] = compute_dTdt(i)
    # compute temperature gradient at right boundary node
    dTdt[0][n-1] = alpha*(-(T[0][n-1]-T[0][n-2])/dx**2 + (Tright-T[0][n-1])/dx**2)  

    # update the temperature vector with the time-stepped value
    T = T + dTdt*timestep

    # plot the updated temperature of the block
    plot_solution(style)
