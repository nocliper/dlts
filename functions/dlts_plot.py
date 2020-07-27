def dlts_plot(T, Time, C, T1, T2, n_windows, ax):

    import matplotlib.pyplot as plt
    from matplotlib import cm
    import numpy as np

    ax.set_ylabel('DLTS $arb. units$')
    ax.set_xlabel('Temperature $T, K$')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    DLTS = []

    for j in range(0,n_windows-1):
        DLTSx = []
        for i in range(len(T)):
            x = (C[i][T2[j]] - C[i][T1[j]])/C[i][-1]
            DLTSx.append(x)
        DLTS.append(DLTSx)

    DLTS = np.asarray(DLTS)

    for i in range(0,n_windows-1):

        ax.plot(T,DLTS[i], label = r'$\tau = %.1f s$'%(Time[T1[i]]-Time[T2[i]]))

    ax.legend()
    ax.grid(True,ls="-")
    return DLTS
