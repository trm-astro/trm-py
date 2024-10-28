"""
Class to help with plots
"""

import trm.subs as subs
from ppgplot import *

def defcol():
    """
    Sets up default pgplot colours
    """
    pgscr(0,1,1,1)
    pgscr(1,0,0,0)
    pgscr(2,0.5,0,0)
    pgscr(3,0,0.5,0)
    pgscr(4,0,0,0.5)
    
class MultiPanel(object):
    """
    This class is to help when setting the world coordinates for
    multi-panel plots. You create a MultiPanel giving it the outer
    limits of the panels, the number of panels and the amount of
    gap between them in x and y as a fraction of the total, then you
    can call 'xylims' to return the bounds for any individual panel.
    """
    
    def __init__(self, xv1, xv2, yv1, yv2, nx, ny, xgap=0., ygap=0.):
        self.xv1  = xv1
        self.xv2  = xv2
        self.yv1  = yv1
        self.yv2  = yv2
        if nx < 1:
            raise MultiPanelError('MultiPanel.__init__: cannot have nx < 1')
        self.nx   = nx
        if ny < 1:
            raise MultiPanelError('MultiPanel.__init__: cannot have ny < 1')
        self.ny   = ny
        self.xgap = xgap
        self.ygap = ygap

    def xylims(self, ix, iy):
        """
        (x1,x2,y1,y2) = self.xylims(ix,iy) returns panel limits of panel
        ix, iy, starting from 1,1 in the bottom-left.
        """

        xr     = self.xv2-self.xv1
        xpanel = xr*(1.-self.xgap)/self.nx
        if self.nx > 1:
            xg = xr*self.xgap/(self.nx-1)
        else:
            xg = 0.
            
        x1 = self.xv1 + (xpanel+xg)*(ix-1)
        x2 = x1 + xpanel

        yr     = self.yv2-self.yv1
        ypanel = yr*(1.-self.ygap)/self.ny
        if self.ny > 1:
            yg = yr*self.ygap/(self.ny-1)
        else:
            yg = 0.
            
        y1 = self.yv1 + (ypanel+yg)*(iy-1)
        y2 = y1 + ypanel
        
        return (x1,x2,y1,y2)

def gpanel(xv1, xv2, yv1, yv2, nx, ny, xgap=0., ygap=0.):
    """
    generator that returns normalised world coordinates
    of successive panels of a multi-panel plot. full return
    is (x1,x2,y1,y2,xlab,ylab,tlab). x1,x2,y1,y2 are the normalised
    device coordinates as used in a call to pgsvp, xlab, ylab, tlab
    are bools to indicate whether the x and y axes and the top of the
    plot should be labelled given potential overlap by other panels.
    """
    mpanel = MultiPanel(xv1, xv2, yv1, yv2, nx, ny, xgap, ygap)
    iy = 0
    while iy < ny:
        iy += 1
        ix = 0
        while ix < nx:
            ix += 1
            (x1,x2,y1,y2) = mpanel.xylims(ix,iy)
            yield (x1, x2, y1, y2, iy==1, ix==1, iy==ny)

def axes(xlabel,ylabel,tlabel,xlab=True,ylab=True,tlab=True,axcol=4,lcol=2):
    """
    Uses pgplot commands to draw and label axes.
    
    xlabel, ylabel, tlabel -- the labels to use
    xlab, ylab, tlab       -- whether to use them (to interface with gpanel)
    axcol, lcol            -- axis and label colours
    """
    
    if xlab:
        xaxis = 'bcnst'
    else:
        xaxis = 'bcst'
    if ylab:
        yaxis = 'bcnst'
    else:
        yaxis = 'bcst'
    pgsci(axcol)
    pgbox(xaxis,0,0,yaxis,0,0)
    pgsci(lcol)
    pglab(xlabel if xlab else '', ylabel if ylab else '', tlabel if tlab else '')
    
class MultiPanelError(subs.SubsError):
    pass


        
