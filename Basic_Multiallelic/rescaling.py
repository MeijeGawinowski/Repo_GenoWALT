def rescaling(init,real_min,real_max):
    init_min = min(init)
    if init_min < 0:
        init = list(np.add(-init_min,init))
    y = [(real_max-real_min)/(max(init)-min(init))*init[i]+real_min for i in range(len(init))]
    return(y)