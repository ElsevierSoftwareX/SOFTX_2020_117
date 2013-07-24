import numpy as np
import pylab as pl

def plotReducedPerformanceData(perfdata, ordered=False):
    labels = []
    times = []

    for pd in perfdata:
        if pd[0] in labels:
            times[labels.index(pd[0])] += pd[3]
        else:
            labels.append(pd[0])
            times.append(pd[3])

    if ordered:
        times,labels = (list(t) for t in zip(*sorted(zip(times,labels))))
            
    ind = np.arange(len(labels))

    fig = pl.figure()
    plt = fig.add_subplot(111)
    plt.barh(ind, times, height=0.5, log=True, align='center')
    pl.yticks(ind, labels)
    pl.xlabel("Time [s]")
    pl.title("Timing data for FINESSE")
    
    return labels, times, fig