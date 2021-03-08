def dlts_plot(T, Time, C, T1, T2, n_windows, ax, Smooth, Bounds):

    import matplotlib.pyplot as plt
    from matplotlib import cm
    import numpy as np
    from scipy.signal import savgol_filter

    ax.set_ylabel(r'$\Delta C/C_0$ $arb. units$')
    ax.set_xlabel('Temperature $T, K$')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    ax.axvspan(Bounds[0], Bounds[1], alpha=0.1, color='red')

    DLTS = []

    for j in range(n_windows):
        DLTSx = []
        for i in range(0, len(T)):
            x = (C[i][T2[j]] - C[i][T1[j]])/C[i][-1]
            DLTSx.append(x)

        if Smooth == 1:
            DLTS.append(DLTSx)
        else:
            DLTS.append(savgol_filter(DLTSx, Smooth, 1))

    DLTS = np.asarray(DLTS)

    c = cm.gnuplot(np.linspace(0, 0.9, n_windows))

    for i in range(n_windows):
        tau = (Time[T2[i]]-Time[T1[i]])/np.log(Time[T2[i]]/Time[T1[i]])
        ax.plot(T, DLTS[i], c = c[i], label = r'$\tau = %.3f$ s'%(tau))

    if n_windows <= 11:
        ax.legend()

    ax.grid(True,ls="-")
    plt.tight_layout()
    return DLTS
