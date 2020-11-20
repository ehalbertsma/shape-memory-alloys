import matplotlib.pyplot as plt

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