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
        C.append(csv[i,2:-2]) # nF

    for i in range(len(C[0])):
        Time.append(dt*(i+2)/1000) # s
    Time = np.asarray(Time)

    return np.array(T), np.array(Time), np.array(C)
