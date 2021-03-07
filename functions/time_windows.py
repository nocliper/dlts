
def time_windows(n_windows, X, Time, C, ax):
    """
    n_windows – number of time windows
    X – fixed relation t2/t1
    Time – time-domain points array
    C – capacitance data
    """

    import matplotlib.pyplot as plt
    from matplotlib import cm
    import numpy as np

    ax.set_ylabel('Capacitance $C, nF$')
    ax.set_xlabel('Time $t, s$')
    ax.grid(True)

    c = cm.gnuplot(np.linspace(0,1,len(C)))

    for i in range(len(C)):
        ax.plot(Time,C[i], c=c[i])

    T1 = []
    T2 = []

    for j in range(n_windows):
        ax.axvline(x = Time[j*1], color=str(j/n_windows), linestyle =':')
        T1.append(j)
        for i in range(len(Time)):
            if Time[i] >= X*Time[j*1]:
                ax.axvline(x = Time[i], color=str(j/n_windows), linestyle =':')
                T2.append(i)
                break
    #print(np.array(Time[T2])/np.array(Time[T1]))
    if len(T1)!=len(T2):
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!\n'+
            '!!!!!!!!!!!!!!!!!!!!!!!!!!\n'+
            '!!!! Too many windows !!!!\n'+
            '!!!!!!!!!!!!!!!!!!!!!!!!!!\n'+
            '!!!!!!!!!!!!!!!!!!!!!!!!!!')

    return T1, T2
