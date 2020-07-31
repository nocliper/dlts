def find_max(T_start,T_stop, DLTS, T, Time, T1, T2, X, n_windows, Doping, A_h):
    '''Returns Ea and crossection of trap
    in temperature interval T_start -> T_stop'''

    import matplotlib.pyplot as plt
    from matplotlib import cm
    import numpy as np
    from scipy.optimize import curve_fit

    start = 0
    stop  = 0

    for i in range(len(T)):

        if T_start >= T[i]:
            start = i
        elif T_stop >= T[i]:
            stop = i

    #print('start=',T[start])
    #print('stop=',T[stop])
    from scipy.signal import argrelextrema

    Temperature = []
    Sx          = []


    for i in range(n_windows-1):
        loc_max = argrelextrema(DLTS[i][start:stop], np.greater, order = 25)[0]
        if not loc_max:
            loc_max = 0
        else:
            Temperature.append(T[start:stop][loc_max[0]])
            temp = Time[T2[i]]-Time[T1[i]]
            Sx.append(temp)

    Sx  = np.array(Sx)
    Tx  = np.array(Temperature)
    Sx  = Sx/np.log(X) #X = tau2/tau1
    Sx  = Sx**(-1)*Tx**(-2)

    fig, axs = plt.subplots(1, 2)
    fig.set_size_inches(10,5)

    axs[0].set_ylabel('DLTS $rel.$')
    axs[0].set_xlabel('Temperature $T, K$')
    axs[0].grid()
    #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    for i in range(len(DLTS)):
        axs[0].plot(T[start:stop],DLTS[i][start:stop])

    for i in range(len(DLTS)):
        loc_max = argrelextrema(DLTS[i][start:stop], np.greater, order = 25)[0]
        if not loc_max:
            loc_max = 0
        else:
            axs[0].plot(T[start:stop][loc_max[0]],DLTS[i][start:stop][loc_max[0]],'x')


    imax = 0
    for i in range((len(DLTS))):
        a = np.amax(np.abs(DLTS[i][start:stop]))
        if imax <= a:
            imax = a

    N_t = 2*Doping*a*X**(X/(X-1))/(X-1)


    ##Fitting curve

    def func(x, a, b):
        return -a*x + b

    popt, pcov = curve_fit(func, 1/Tx, np.log(Sx)) # your data x, y to fit
    #print(popt)

    T_average = (T_start + T_stop)/2

    Error = popt[0]*(8.617*10**-5)*(-2)*(T_stop-T_start)/(1.15*T_average*(np.log10(T_start**2/(6*T_stop**2))))##Bridman p.22 eq. 1.15

    axs[1].set_ylabel(r'$\ln{(\tau^{-1}\cdot T^{-2})}, s^{-1}\cdot K^{-2}$')
    axs[1].set_xlabel(r'Temperature $1/T, K^{-1}$')
    axs[1].plot(1/Tx,np.log(Sx),'o', label='Original Data')
    axs[1].plot(1/Tx,popt[1]-popt[0]/Tx,'k-',label='Fitted line')
    axs[1].legend()
    axs[1].grid()

    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.tight_layout()

    from tabulate import tabulate
    Table = np.array([Tx, Sx])
    t = tabulate(Table.T, headers=['1000/T, K', 'e/T^2, s-1*K-2'], tablefmt='orgtbl')
    print(t)

    from IPython.display import display, Math
    display(Math(r'$$E_C - E_t = ({} \pm {})eV$$'.format(round(popt[0]*(8.617*10**-5),3),round(Error,3))))
    display(Math(r'$$N_t = %.2E cm^2$$'%N_t))
    display(Math(r'$$\sigma_p = %.2E cm^2$$'%(np.exp(popt[1])/A_h)))
