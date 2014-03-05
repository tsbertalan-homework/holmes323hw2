
import logging
logging.basicConfig(format='%(levelname)s: %(message)s',
                    level=logging.INFO, # Turn on/off debug by switching INFO/ERROR.
                    name='log')

import numpy as np

def plotNullclines(ax, g=4, b=0.2, s1=.5, s2=.5, lims = (-1, 1), showFlipped=True):
    
    f = lambda x, g, b: 1 / (1 + np.exp(-4*g*(x-b)))

    # nullcline 1
    color = 'b'
    N =200
    x2 = np.linspace(lims[0], lims[1], N)
    x1 = -f(x2, g, b) + s1
    ax.plot(x1, x2, color+'-', label=r'$\nu_1(x_2)$')
    if showFlipped:
        ax.plot(x2, x1, color+'-.', label=r'$\nu_1^{-1}(x_1)$')

    # nullcline 2
    color = 'r'
    x1 = np.linspace(lims[0], lims[1], N)
    x2 = -f(x1, g, b) + s2
    ax.plot(x1, x2, color+'-', label=r'$\nu_2(x_1)$',zorder=-10)
    if showFlipped:
        ax.plot(x2, x1, color+'-.', label=r'$\nu_2^{-1}(x_2)$')


    ax.plot(x1, x1, 'k', lw=.1)
    for ticker in ax.set_xticks, ax.set_yticks:
        ticker(np.arange(lims[0], lims[1], .5))
    
    ax.set_title(r'$g=%.0f$, $\beta=%.1f$, $s_1=%.2f$, $s_2=%.2f$' 
                        % (g, b, s1, s2))
    ax.set_xlabel('$x_1$')
    ax.set_ylabel('$x_2$')
    
def findPeaks(traces):
    '''
    Find all the peaks in a (set of?) voltage trace(s).
    Returns their indices.
    
    >>> T = np.arange(0, 10, 100)
    >>> X1 = np.sin(T)
    >>> X2 = np.sin(T + .1)
    >>> traces = np.hstack((X1, X2))
    >>> traces.shape
    >>> findPeaks(traces)
    '''
    from scipy.signal import argrelmax
    if len(traces.shape) == 1:
        traces = traces.reshape((traces.size, 1))
    return [argrelmax(traces[:,i], axis=0, order=2)[0] for i in range(traces.shape[1])]
