
def band_plot(Eg, Ts, Ns, Es, size = 12, thick = 0.05, border = 0.085, c = 'blue'):

    import matplotlib.pyplot as plt
    import numpy as np


    x  = np.linspace(0, 1, 20)
    Ec = -x**4/4

    fig = plt.figure(figsize=(7,5))
    ax  = fig.add_subplot(111)
    ax.set_ylim(-0.11*max(Ns), 1.1*max(Ns))
    ax.set_xlim(0-Eg*0.1, Eg*1.1)
    ax.axhline(y = 0, xmin= border, xmax= 1-border, ls = '-', c = 'k', lw = 2)

    # Bands
    ax.plot(Ec, x*max(Ns),         lw = 2, c = 'k')
    ax.plot(-Ec + Eg, (x)*max(Ns), lw = 2, c = 'k')

    # Annotations
    for i, v in enumerate(Ts):
        ax.annotate('%.4s'%v , xy = (Es[i], Ns[i]+max(Ns)*0.05), rotation = 90, size = size, ha = 'center')

    ax.annotate(r'E$_C$', xy = (Ec[-1],       max(Ns)*1.05),     rotation = 90, size = size+4,
                va = 'center', ha = 'center')
    ax.annotate(r'E$_V$', xy = (-Ec[-1] + Eg, max(Ns)*1.05), rotation = 90, size = size+4,
                va = 'center', ha = 'center')

    ax.annotate(s='', xy=(0, -max(Ns)*0.05), xytext=(Eg, -max(Ns)*0.05), arrowprops=dict(arrowstyle='<|-|>'), size=size+6, color = 'k')
    ax.annotate(r'E$_G$ = '+str(Eg)+' eV', xy=(Eg/2, -max(Ns)*0.1), rotation = -180, size = size+4, ha = 'center')

    edgecolor = ''
    ax.bar(Es, Ns, thick, edgecolor=edgecolor, color=c)

    # Scale
    scale = np.log10(Ns)+0.2
    scale = np.average(scale)//1
    print(scale)

    ax.bar(Eg*1.05, 10**scale, thick, edgecolor=edgecolor,  color = 'red')
    ax.annotate('~ 10$^{%.0f}$ ${cm}^{-3}$'%(scale) , xy = (Eg*1.05, 10**scale+max(Ns)*0.05), rotation = 90, size = size, ha = 'center')
    ###

    plt.axis('off')
    plt.tight_layout()
    plt.savefig('traps.svg')
