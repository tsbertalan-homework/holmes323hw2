
import numpy as np
import kevrekidis as kv  # In-lab package. Lots of convenience functions that
                         #  aren't really needed here, but that I'm used to using. 
from utilities import plotNullclines, logging
if __name__ != "__main__":
    logging.disable(logging.INFO)  # Another way of turning of debug. Unnecessary.


def problem1(fname='hw2-problem1.pdf'):
    #make some figure objects
    fig = kv.plotting.plt.figure(figsize=(4,13))
    A = [fig.add_subplot(4,1,1),
            fig.add_subplot(4,1,2),
            fig.add_subplot(4,1,3),
            fig.add_subplot(4,1,4),
    ]           

    # add nullclines
    plotNullclines(A[0], showFlipped=False)
    plotNullclines(A[1],
                   s1=.55, s2=.45)
    plotNullclines(A[2],
                   s1=.75, s2=.25)
    plotNullclines(A[3],
                   s1=.99, s2=.01)
    
    # further manipulate the figures
    for i in range(len(A)):
        if i != len(A)-1:
            A[i].set_xticklabels([])
            A[i].set_xlabel('')
        if i == 0:
            A[i].legend(loc='lower left', prop={'size':10})
    fig.subplots_adjust(left=.2, bottom=.05)
    fig.savefig(fname)


def problem1Bifurc(fname='hw2-problem1Bifurc.pdf'):
    fig, ax = kv.fa(figsize=(4,2))  # figures

    # a cubic equation
    y = lambda x: -(x - 1) * (x + 1) * x
    ax.set_xlim(-1,1)
    s3 = np.sqrt(1./3)  # regions for changing plotting style
    offset = 0
    stablex1 = np.linspace(-2., -s3, 400) + offset
    unstablex = np.linspace(-s3, s3, 400) + offset
    stablex2 = np.linspace(s3, 2., 400) + offset
    for (x, style) in zip(
                        (stablex1, unstablex, stablex2),
                        ('k', 'r--', 'k')
                    ):
        ax.plot(y(x-offset), x, style)
    
    # make everything wonderful
    ax.set_xlim(0,.5)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(r'$\Delta s$'); ax.set_ylabel(r'$\tilde x_2$')
    ax.set_xticklabels
    fig.savefig(fname)


def problem2(fname='hw2-problem2.pdf',
             g=1.7, b=0.2, s1=0.0, s2=0.0, N=32):
    f = lambda x: 1 / (1 + np.exp(-4*g*(x-b)))
    
    # domain
    x1lim = x2lim = -2, 2
    x1, x2 = np.meshgrid(
        np.linspace(x1lim[0], x1lim[1], N),
        np.linspace(x2lim[0], x2lim[1], N),
    )
    
    # construct V
    integ = lambda x: -(
                        4*(-1 + 1/
                                  (1+np.exp(4*g*(x - b)))
                           )*g*x
                        + 
                        np.log(1+np.exp(4*g* (x-b)))
                        ) / (4*g)
    V = integ(x1) + integ(x2) + f(x1) * f(x2)
    
    # construct vdot
    Dfdx = lambda x: (4* np.exp(-4* g *(x-b))* g)/\
                     (1+np.exp(-4* g* (x-b)))**2
    xdot = lambda xi, xj, si: si - f(xj) - xi
    Vdot = Dfdx(x1) * xdot(x1, x2, s1) * (x1 + f(x1)) +\
           Dfdx(x2) * xdot(x2, x1, s2) * (x2 + f(x2))
           
    
    # plot V and Vdot 
    fig = kv.plotting.plt.figure(figsize=(4,10))
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
    for ax, mat in zip((ax1, ax2),(Vdot, V)):
        
    #fig, ax = kv.showMat(V, extent=x1lim+x2lim, zorder=-30)
        kv.showMat(mat, extent=x1lim+x2lim, zorder=-30, figax=(fig, ax), cmap='jet',
                   alpha=.9, interpolation='bicubic')
        # Add the nullclines to both plots
        plotNullclines(ax, g=g, b=b, showFlipped=False, s1=s1, s2=s2, lims=x1lim)
    
        # Set titles and labels; save
        ax.set_xlim(x1lim)
        ax.set_ylim(x2lim)
        ax.set_yticks(np.linspace(x2lim[0], x2lim[1], 2))
        
    ax1.set_xticks([])
    ax2.set_xticks(np.linspace(x1lim[0], x1lim[1], 2))
    
    ax1.set_title(r'$\frac{d V(x_1, x_2)}{dt}$')
    ax2.set_title('$V(x_1,x_2)$')
    fig.subplots_adjust(left=.2, top=.95, bottom=.05)
    fig.savefig(fname)



def treatment(fname, fig, axes):
    '''some tasks for both of problem 3's figures'''
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
    axes[-1].set_ylabel('$r^e$')
    axes[-1].set_xlabel(r'$\mu$')
    fig.savefig(fname)


def problem3o3(fname='hw2-problem3-o3.pdf'):
    fig = kv.plotting.plt.figure()
    ax = fig.add_subplot(2, 1, 1)
    axes = [fig.add_subplot(2,2,i) for i in 3,4]; axes.append(ax)  # make a grid of axes
    
    f = lambda x: -x**2  # a quadratic pair of branches (QUALITATIVE! probably)
    #Plot the branches
    x = np.linspace(-1,1, 100) 
    ax.plot(f(x), x, 'k--')
    # plot the horizontal axis
    x = (-2,f(0))
    ax.plot(x, (0,0), 'k-')
    x = (f(0),1)
    ax.plot(x, (0,0), 'k--')
    ax.set_xlim(-1,1)
    
    treatment(fname, fig, axes)


def problem3o5(fname='hw2-problem3-o5.pdf'):
    fig = kv.plotting.plt.figure()
    ax = fig.add_subplot(2, 1, 1)
    axes = [fig.add_subplot(2,3,i) for i in 4,5,6]; axes.append(ax)
    borders = -3, -np.sqrt(5./2), np.sqrt(5./2), 3  # regions for switching plotting style
    f = lambda x: 4 - 5*x**2 + x**4
    styles = 'k-', 'k--'
    
    # Plot the branches
    for i in range(len(borders)-1):
        x = np.linspace(borders[i], borders[i+1], 100) 
        ax.plot(f(x), x, styles[i%2])
    
    # Plot the horizontal axis
    x = (-5,f(0))
    ax.plot(x, (0,0), 'k-')
    x = (f(0),6)
    ax.plot(x, (0,0), 'k--')
    ax.set_xlim(-5,6)
    treatment(fname, fig, axes)


def problem5fig1(fig, fname='hw2-problem5-Vtraces.pdf'):
    fig.savefig(fname)  # this is very ugly, I know. The point of doing it this
                        #  way is that I can call this code from LaTeX
                        #  (via ctan.org/tex-archive/macros/latex/contrib/python),
                        #  and, for organization's sake, I
                        #  am trying to have one function per figure, so I can
                        #  call generate it directly there, with the option
                        #  to pass in parameters if desired.
    
def problem5fig2(fig, fname='hw2-problem5-DvsL.pdf'):
    fig.savefig(fname)
    
def problem5():
    from scipy.optimize import newton_krylov  # Simple Newton fails horribly.
    from utilities import findPeaks
    def Dxdt(lbase, D, tD):
        lplus = lbase + D
        def dxdt(X, t):
            if t > tD:
                l = lplus
            else:
                l = lbase
            c, h, b, a, p, g = tuple(X)
            return np.array([
                             (-c-p*h+l) / 10.,
                             (-h+c)/100.,
                             (-b+(6.*c-5.*h)/(1.+9.*a)),
                             (-a+b)/80.,
                             (-p+a/10.)/4000.,
                             (-g+50.*b/(13. + b))/10.
                             ])
        return dxdt
    
    ## We will plot...
    #   ... the G(t) traces ...
    fig1, ax1 = kv.fa(figsize=(8,2.5))
    #   ... and the D vs L dependence.
    fig2, ax2 = kv.fa(figsize=(8,4))
    
    
    lvals = 1, 2, 10, 20, 100, 200, 1000, 2000, 10000
    colors = kv.plotting.colors(numColors=len(lvals))
    Dvals = []
    for i in range(len(lvals)):
        l = lvals[i]
        logging.info("l = %.0f" % l)
        color = colors[i]
        ## FIND A FP
        # First get b as a zero of a cubic, then find the other SS values
        # Using 2.5.5 and adjacent text, for this value of l. 
        b = newton_krylov(lambda b: 9*b**3+91*b**2+10*b-10*l, 10) # nk is overkill here, obviously
        c = 10*l/(10+b); h=c; a=b; p=a/10; g=10*b/(13+b) 
        # Clean this up a bit (perhaps unnecessary) by ensuring that it's
        #  actually a fixed point of the dynamics from Dxdt:
        Dguess = l
        dxdtwt = Dxdt(l, Dguess, 1)
        fpF = lambda X: dxdtwt(X, 0)
        initial = newton_krylov(fpF, (c,h,b,a,p,g))
        #initial = (c,h,b,a,p,g)  # <- This also sortof works, but not as well.
        def doTest(D):
            ## Integrate for this value of D. This uses scipy.odeint
            X, T = kv.integrate(initial, Dxdt(l, D, 10),
                                tmin=0, tmax=200, giveTime=True)
            ## Find the peak value of g.
            peaks = findPeaks(X[:,-1])[0]
            peakInd = peaks[np.argmax([X[peak,-1] for peak in peaks])]
            peakg = X[peakInd,-1]
            return X, T, peakInd, peakg
        def zeroMe(D):
            D = abs(D)  # Counteract excessive cleverness from the solver.
            X, T, peakInd, peakg = doTest(D)
            spikeMag = peakg - initial[-1]
            # spikeMag = abs(peakg - initial[-1])  # this gives the same results
            val = spikeMag - 1.0  # Optimize D s.t. this is zero.
            return val
        D = abs(  # We only allowed semipositive values up there.
                newton_krylov(zeroMe, Dguess)
                )
        Dvals.append(D)
        X, T, peakInd, peakg = doTest(D)
        label = None
        if l==lvals[0] or l == lvals[-1]:
            label = "$l=%.0f$" % l
        ax1.plot(T, X[:,-1] - initial[-1], label=label, lw=4, color=color)
        ax1.axvline(x=T[peakInd], color=color, lw=1)
        ax1.set_xlabel("$t$ [ms]")
        ax1.set_ylabel("$G(t) - G(0)$ [Hz?]")
        
    # Do some more figure manipulations
    fig1.subplots_adjust(right=.95, bottom=.2, left=.1)
    ax1.legend()
    
    ax2.scatter(lvals, Dvals, color='black', label='simulated')
    ax2.set_xlabel('$L$')
    ax2.set_ylabel('$\Delta$')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    fig2.subplots_adjust(right=.95, bottom=.1, left=.1)
    
    L = np.logspace(np.log10(min(lvals)), np.log10(max(lvals)), num=50)
    
    # Let's add some polynomial "trendlines", then wish we hadn't:
    def addPoly(x0, y0, ord, labelord):
        x = L
        a = np.log10(x)
        a0 = np.log10(x0)
        b0 = np.log10(y0)
        b = ord * (a - a0) + b0
        y = 10 ** b
        ax2.plot(x, y, lw=.5, label=r'$\log_{10}(\Delta)=%s\log_{10}(L)$' % labelord)
    addPoly(min(lvals), min(Dvals), 0.33333, r'\frac{1}{3}')
    addPoly(max(lvals), max(Dvals), 1.0, '1')
    
    # Compare to the approxmation in Eqn. (2.5.7).
    Dapprox = (13.+L)**2/(50.*13.)
    ax2.plot(L, Dapprox, 'k--', label='eqn (2.5.7)')
    ax2.set_ylim(min((
                     min(Dvals),
                     min(Dapprox),
                     ))/2,
                     max(Dvals)*2,
                 )
    ax2.legend(loc='upper left')
    return fig1, fig2
    
    
if __name__ == "__main__":
    # If this is being run as a script, do and show everything.
    problem1()
#     problem1Bifurc()
#     problem2()
#     problem3o3()
#     problem3o5()
#     fig1, fig2 = problem5()
#     problem5fig1(fig1)
#     problem5fig2(fig2)
    
    kv.plotting.show()