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

math = Model()

def f(x):
    return 2

X = np.linspace(0,10,100)

I = math.integrate(f,*X)
print(I)