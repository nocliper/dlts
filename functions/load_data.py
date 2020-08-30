def load_data(path, dt = 150):
    '''awdawd'''
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import numpy as np

    csv   = np.genfromtxt(path, delimiter="\t")

    Time   = [] #s
    T      = csv[0:,0] #K
    C      = []

    for i in range(0,len(T)):
        C.append(csv[i,3:-3]) # nF


    for i in range(len(C[0])):
        Time.append(dt*(i+1)/1000) # s
    Time = np.asarray(Time)

    C = np.asarray(C)
    ## Sorting

    T = T.tolist()
    C = C.tolist()

    L = sorted(zip(T,C))

    new_T, new_C = zip(*L)

    return np.array(new_T), np.array(Time), np.array(new_C)
