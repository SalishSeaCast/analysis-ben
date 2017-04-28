import numpy                as np
from graphviz             import Source


def calc_rho(Sal, TempC, P):
    """ Calculate rho: Based on SOG code
    """
    
    # Calculate the square root of the salinities
    sqrSal = np.sqrt(Sal)

    # Calculate the density profile at the grid point depths
    # Pure water density at atmospheric pressure
    # (Bigg P.H., (1967) Br. J. Applied Physics 8 pp 521-537)
    R1 = ((((6.536332e-9 * TempC - 1.120083e-6) * TempC + 1.001685e-4)
           * TempC - 9.095290e-3) * TempC + 6.793952e-2) * TempC - 28.263737
    R2 = (((5.3875e-9 * TempC - 8.2467e-7) * TempC + 7.6438e-5)
          * TempC - 4.0899e-3) * TempC + 8.24493e-1
    R3 = (-1.6546e-6 * TempC + 1.0227e-4) * TempC - 5.72466e-3

    # International one-atmosphere equation of state of seawater
    SIG = (4.8314e-4 * Sal + R3 * sqrSal + R2) * Sal + R1

    # Specific volume at atmospheric pressure
    V350P = 1.0 / 1028.1063
    SVA   = -SIG * V350P / (1028.1063 + SIG)

    # Density anomoly at atmospheric pressure
    rho = 28.106331 - SVA / (V350P * (V350P + SVA)) + 1000
    
    return rho


def flow_diagram():
    """
    """
    
    diagram = Source('digraph G { ' +
            'rankdir=TB ' +
            'node[shape=circle]; a1 b1 c1 d1 ' +
            'node[shape=box] ' +
            'a1 -> b1 -> c1 -> d1 ' +
            'b1 -> b2 -> b1 ' +
            'c1 -> c2 -> c3 -> c2 -> c1 ' +
            'a1 [label="Load Data"] ' +
            'b1 [label="Make Figure"] ' +
            'c1 [label="Animate"] ' +
            'd1 [label="Save"] ' +
            'b2 [label="def make_figure( ):"] ' +
            'c2 [label="def next_frame( t ):"] ' +
            'c3 [label="def update_figure( ):"] ' +
            '{rank=same; b1 b2} ' +
            '{rank=same; c1 c2 c3}}')
    
    return diagram
