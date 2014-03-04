
import numpy as np
import kevrekidis as kv
from utilities import plotNullclines, logging
if __name__ != "__main__":
    logging.disable(logging.INFO)


def problem1(fname='hw2-problem1.pdf'):
    #fig, A = kv.fa(numAxes=3)
    fig = kv.plotting.plt.figure(figsize=(4,13))
    A = [fig.add_subplot(4,1,1),
            fig.add_subplot(4,1,2),
            fig.add_subplot(4,1,3),
            fig.add_subplot(4,1,4),
    ]           

    plotNullclines(A[0])
    #plotNullclines(A[1], b=.75)
    #plotNullclines(A[2], b=.75, ds=0.05)
    plotNullclines(A[1],
                   s1 = .5 + .1 / 2,
                   s2 = 1 - .5 + .1 / 2)
    plotNullclines(A[2],
                   s1 = .5 + .5 / 2,
                   s2 = 1 - .5 + .5 / 2)
    plotNullclines(A[3],
                   s1 = .5 + 0.98 / 2,
                   s2 = 1 - .5 + 0.98 / 2)
    for i in range(len(A)):
        if i != len(A)-1:
            A[i].set_xticklabels([])
            A[i].set_xlabel('')
        if i == 0:
            A[i].legend(loc='lower left', prop={'size':10})
    fig.subplots_adjust(left=.2)
    fig.savefig(fname)


def problem1Bifurc(fname='hw2-problem1Bifurc.pdf'):
    fig, ax = kv.fa(figsize=(4,2))

    y = lambda x: -(x - 1) * (x + 1) * x
    ax.set_xlim(-1,1)
    s3 = np.sqrt(1./3)
    offset = 0
    stablex1 = np.linspace(-2., -s3, 400) + offset
    unstablex = np.linspace(-s3, s3, 400) + offset
    stablex2 = np.linspace(s3, 2., 400) + offset
    pairs = zip(
                        (stablex1, unstablex, stablex2),
                        ('k', 'r--', 'k')
                    )
    for (x, style) in pairs:
        ax.plot(y(x-offset), x, style)
    
    ax.set_xlim(0,.5)
    #ax.set_ylim(-1,1)
    #ax.axvline(0, color='gray', lw=.01)
    #ax.axhline(0, color='gray', lw=.01)
    ax.set_xticks([]); ax.set_yticks([])
    #ax.set_xticklabels(['0', '1'])
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
           
    
    # plot ... 
    fig = kv.plotting.plt.figure(figsize=(4,10))
    #  ... Vdot
    ax1 = fig.add_subplot(2,1,1)
    #fig, ax = kv.showMat(V, extent=x1lim+x2lim, zorder=-30)
    kv.showMat(Vdot, extent=x1lim+x2lim, zorder=-30, figax=(fig, ax1), cmap='RdYlBu',
               alpha=.9, interpolation='bicubic')
    
    # show contour lines
#     CS = ax1.contour(x1, x2, VdotSimp)
#     kv.plotting.plt.clabel(CS, fontsize=9, inline=1)
#     kv.plotting.plt.clabel(CS, fontsize=9, inline=1)

    #  ... V
    from mpl_toolkits.mplot3d import axes3d
    ax2 = fig.add_subplot(2,1,2,projection='3d')
    ax2.plot_surface(x1, x2, V, rstride=1, cstride=1,
                 linewidth=0, antialiased=True, cmap=kv.plotting.getCmap('jet'),
                 alpha=.7)
    
#     ax2.contour(x1, x2, V)
    
    # Add the nullclines to both plots
    plotNullclines(ax1, g=g, b=b, showFlipped=False, s1=s1, s2=s2, lims=x1lim)
    plotNullclines(ax2, g=g, b=b, showFlipped=False, s1=s1, s2=s2, lims=x1lim)
    
    # Set titles and labels; save
    for ax in ax1, ax2:
        ax.set_xticks(np.linspace(x1lim[0], x1lim[1], 2))
        ax.set_yticks(np.linspace(x2lim[0], x2lim[1], 2))
    
    ax2.set_zticks(np.linspace(V.min(), V.max(), 2))
    ax2.view_init(elev=20, azim=250)
    
    ax1.set_title(r'$\frac{d V(x_1, x_2)}{dt}$')
    ax2.set_title('$V(x_1,x_2)$')
    fig.subplots_adjust(left=.2, top=.95, bottom=.02)
    ax1.set_xlim(x1lim)
    ax1.set_ylim(x2lim)
    ax2.set_xlim(x1lim)
    ax2.set_ylim(x2lim)
    fig.savefig(fname)


def problem3(fname='hw2-problem3.pdf'):
    def Dx(mu, om, al):
        def dx(X):
            x = X[0]
            y = X[1]
            out = np.array([
             mu*x-om*y+al*(x**2 + y**2)*x,
             mu*y+om*x+al*(x**2 + y**2)*y,
                             ])
            return out
        return dx
    fig = kv.plotting.plt.figure()
    bifurcAx = fig.add_subplot(2, 1, 1)
    phaseAx1 = fig.add_subplot(2, 2, 3)
    phaseAx2 = fig.add_subplot(2, 2, 4)
    om = 1
    al = .0001
    for mu, ax, tmax in zip((-1, 1), (phaseAx1, phaseAx2), (64, 1)):
        dx = Dx(mu, om, al)
        NX = NY = 8
        width = 4
        for x0 in np.linspace(-width/2, width/2, NX):
            for y0 in np.linspace(-width/2, width/2, NY):
                X, T = kv.integrate((x0, y0), dx, tmin=0, tmax=tmax, nstep=1e3)
                ax.plot(X[:,0], X[:,1], 'k')
                ax.scatter(x0, y0, color='green')
        ax.set_xlim(-width/2, width/2)
        ax.set_ylim(-width/2, width/2)
    fig.suptitle(r'$\omega=%.1f$, $\alpha=%.1f$' % (om, al))
    fig.savefig(fname)
        

def problem5fig1(fig, fname='hw2-problem5-Vtraces.pdf'):
    fig.savefig(fname)
    
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
            ## Integrate for this value of D. 
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
    fig1.subplots_adjust(right=.95, bottom=.2, left=.1)
    ax1.legend()
    
    ax2.scatter(lvals, Dvals, color='black', label='simulated')
    ax2.set_xlabel('$L$')
    ax2.set_ylabel('$\Delta$')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    fig2.subplots_adjust(right=.95, bottom=.1, left=.1)
    
    L = np.logspace(np.log10(min(lvals)), np.log10(max(lvals)), num=50)
    
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
        
    Dapprox = (13.+L)**2/(50.*13.)
    ax2.plot(L, Dapprox, 'k--', label='eqn (2.5.7)')
    ylims = ax2.get_ylim()
    ax2.set_ylim(min((
                     min(Dvals),
                     min(Dapprox),
                     ))/2,
                     max(Dvals)*2,
                 )
    ax2.legend(loc='upper left')
    #fig2.subplots_adjust()
    return fig1, fig2
    
if __name__ == "__main__":
#     problem1()
#     problem1Bifurc()
    problem2()
#     problem3()
#     fig1, fig2 = problem5()
#     problem5fig1(fig1)
#     problem5fig2(fig2)
    kv.plotting.show()