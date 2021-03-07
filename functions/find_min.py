def find_min(T_start,T_stop, DLTS, T, Time, T1, T2, X, n_windows, Doping, A_h, ax1, ax2):
    '''Returns Ea and crossection of trap
    in temperature interval T_start -> T_stop'''

    import matplotlib.pyplot as plt
    from matplotlib import cm
    import numpy as np
    from scipy.optimize import curve_fit

    c = cm.gnuplot(np.linspace(0, 0.9, n_windows))

    start = 0
    stop  = 0

    for i in range(len(T)):

        if T_start >= T[i]:
            start = i
        elif T_stop >= T[i]:
            stop = i

    from scipy.signal import argrelextrema

    Temperature = []
    Sx          = []

    for i in range(n_windows):
        loc_min = argrelextrema(DLTS[i][start:stop], np.less, order = 25)[0]
        if not loc_min:
            loc_min = 0
        else:
            Temperature.append(T[start:stop][loc_min[0]])
            temp = (Time[T2[i]]-Time[T1[i]])/(np.log(Time[T2[i]]/Time[T1[i]]))
            Sx.append(temp)

    Sx  = np.array(Sx)
    Tx  = np.array(Temperature)
    Sx  = Sx**(-1)*Tx**(-2)


    for i in range(len(DLTS)):
        loc_min = argrelextrema(DLTS[i][start:stop], np.less, order = 20)[0]
        if not loc_min:
            loc_min = 0
        else:
            ax1.plot(T[start:stop][loc_min[0]],DLTS[i][start:stop][loc_min[0]], 'x', c = c[i])

    imax = 0
    for i in range((len(DLTS))):
        a = np.amax(np.abs(DLTS[i][start:stop]))
        if imax <= a:
            imax = a

    N_t = 2*Doping*a*X**(X/(X-1))/(X-1)


    #plt.savefig('DLTS_spectra_scanned'+'.png', format='png', dpi=500)


    ##Fitting curve

    def func(x, a, b):
        return -a*x + b

    popt, pcov = curve_fit(func, 1/Tx, np.log(Sx)) # your data x, y to fit
    perr = np.sqrt(np.diag(pcov)) # a and b error

    T_average = (T_start + T_stop)/2

    ax2.set_ylabel(r'$\ln{(\tau^{-1}\cdot T^{-2})}, s^{-1}\cdot K^{-2}$')
    ax2.set_xlabel(r'Temperature $1/T, K^{-1}$')
    ax2.scatter(1/Tx, np.log(Sx), c = c, label='Original Data')
    ax2.plot(1/Tx, popt[1]-popt[0]/Tx, 'k-', label='Fitted line')
    ax2.legend()
    ax2.grid()

    ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.tight_layout()

    from tabulate import tabulate
    Table = np.array([1000/Tx, Sx])
    t = tabulate(Table.T, headers=['1000/T, K-1', 'e/T^2, s-1*K-2'])
    print(t)

    def sci_notation(a, b):
        a_pow = np.ceil(np.log10(a))-1
        b_pow = np.ceil(np.log10(b))-1
        return a/10**a_pow, b/10**a_pow, a_pow

    sigm, serr, spow = sci_notation(np.exp(popt[1])/A_h, np.exp(perr[1])/A_h)

    from IPython.display import display, Math
    display(Math(r'$$E_t - E_V = ({} \pm {}) eV$$'.format(round(popt[0]*(8.617*10**-5),3), round(perr[0]*(8.617*10**-5),3))))
    display(Math(r'$$N_t = %.2E cm^2$$'%N_t))
    display(Math(r'$$\sigma_p = (%.3f \pm %.3f)\cdot {10}^{%.0f} {cm}^2$$'%(sigm, serr, spow)))
