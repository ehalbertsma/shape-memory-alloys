import numpy as np


class Model(object):

    def __init__(self):
        self.yee = True
        # idk

    def integrate(self, f, *X):
        """
        Integral of the function f on the interval X, which is a 1-D array.
        Left endpoint Riemann sum.

        """
        I = 0
        Y = f(*X)
        for yi in Y:
            I+=yi

        return I

    def sigma(t):
        """
        Returns the stress sigma at time t
        We're using sinusoidal, but we could also use other funtions

        """

        rate = 4/3 # stress rate [MPa/s]
        sigmax = 200 # amplitude of stress [MPa]
        T = 2*sigmax/rate # cycle period [s]

        return sigmax*np.sin((np.pi/T)*t)

    def temperature(t):
        """
        Returns the temperature at time t.
        TODO - requires sawtooth function

        """"
        return 270

    def sawtooth(peak, freq, t):
        """
        IDK

        """
        
        return 


math = Model()

def f(x):
    return 2

X = np.linspace(0,10,100)

I = math.integrate(f,*X)
print(I)