def dlts_plot(T, Time, C, T1, T2, n_windows, ax, Smooth):

    import matplotlib.pyplot as plt
    from matplotlib import cm
    import numpy as np
    from scipy.signal import savgol_filter

    ax.set_ylabel('DLTS $arb. units$')
    ax.set_xlabel('Temperature $T, K$')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    DLTS = []

    for j in range(0,n_windows-1):
        DLTSx = []
        for i in range(len(T)):
            x = (C[i][T2[j]] - C[i][T1[j]])/C[i][-1]
            DLTSx.append(x)
        if Smooth == 1:
            DLTS.append(DLTSx)
        else:
            DLTS.append(savgol_filter(DLTSx, Smooth, 1))

    DLTS = np.asarray(DLTS)

    for i in range(0,n_windows-1):
        tau = (Time[T2[i]]-Time[T1[i]])/np.log(Time[T2[i]]/Time[T1[i]])
        ax.plot(T,DLTS[i], label = r'$\tau = %.3f$ s'%(tau))

    ax.legend()
    ax.grid(True,ls="-")
    return DLTS
