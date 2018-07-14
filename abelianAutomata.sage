load("tim/graph_loop.sage")

#{{{ matrix definitions
"""
All the 2d matrices for ease of use (order taken from Okano's paper)

lists are poly coefficients starting with constant term
"""
A1 = companion_matrix([-1/2,    0, 1], format='left')
A2 = companion_matrix([ 1/2,    0, 1], format='left')
A3 = companion_matrix([ 1/2, -1/2, 1], format='left')
A4 = companion_matrix([ 1/2,  1/2, 1], format='left')
A5 = companion_matrix([ 1/2,   -1, 1], format='left')
A6 = companion_matrix([ 1/2,    1, 1], format='left')

A  = A6  # default A to be the CCC^3_2 matrix

"""
Might as well do the 3d matrices while we're here...
"""
B1  = companion_matrix([-1/2, -1/2,  1/2, 1], format='left')
B2  = companion_matrix([-1/2,    0,    0, 1], format='left')
B3  = companion_matrix([-1/2,    0,  1/2, 1], format='left')
B4  = companion_matrix([-1/2,    0,    1, 1], format='left')
B5  = companion_matrix([-1/2,  1/2, -1/2, 1], format='left')
B6  = companion_matrix([-1/2,  1/2,    0, 1], format='left')
B7  = companion_matrix([-1/2,    1,   -1, 1], format='left')
B8  = companion_matrix([-1/2, -1/2,  1/2, 1], format='left')
B9  = companion_matrix([ 1/2, -1/2, -1/2, 1], format='left')
B10 = companion_matrix([ 1/2,    0, -1/2, 1], format='left')
B11 = companion_matrix([ 1/2,    0,    0, 1], format='left')
B12 = companion_matrix([ 1/2,    0,    0, 1], format='left')
B13 = companion_matrix([ 1/2,  1/2,  1/2, 1], format='left')
B14 = companion_matrix([ 1/2,    1,    1, 1], format='left')
# }}}

#{{{ Stuff for testing

def binLen(n):
    """
    All binary strings of length n
    """
    if n == 0:
        return [""]
    else:
        recur = binLen(n-1)
        return map(lambda xs: "1"+xs, recur) + map(lambda xs: "0"+xs, recur)

#}}}

# plot settings
PLOT_DPI = 200
PLOT_ITERS = 1000

MAX_DEPTH=10000
sys.setrecursionlimit(MAX_DEPTH)

# Define polynomial rings
RZ.<z> = PolynomialRing(ZZ)
RQ.<x> = PolynomialRing(QQ)

class ADiGraph(DiGraph):
    def plot2(self, **kwargs):
        """
        Pretty print a digraph

        Edges have the following coloring:
            0 | 0 -> black
            1 | 1 -> grey
            0 | 1 -> green
            1 | 0 -> blue
        """
        edge_colormap = {"0|0": "black"
                        ,"1|1": "grey"
                        ,"0|1": "green"
                        ,"1|0": "blue"
                        }

        return self.graphplot(layout         = 'spring'
                             ,color_by_label = edge_colormap
                             ,vertex_color   = 'white'
                             ,iterations     = PLOT_ITERS
                             ,dpi            = PLOT_DPI
                             ,**kwargs
                             ).plot()


class CompleteAutomaton(object):
    def __init__(self, A, r):
        """
        Construct the Complete Automaton for a given (A,r) pair.
        """
        assert A.dimensions()[0] == len(r)
        self.A   = A
        self.r   = vector(r)

        self.m   = len(self.r)
        self.chi = self.A.charpoly()

        self.Ai   = self.A.inverse()
        self.chii = RZ(self.Ai.charpoly())

        self.endo = RQ.quo(self.chii)
        self.x    = self.endo.gens()[0]  # the variable for self.R


    def scaleByPoly(self, p):
        """
        Return a new CompleteAutomaton which is this one scaled by p
        """
        p = RZ(list(p))
        q = RZ(list(self.r))

        return CompleteAutomaton(self.A, vector((p * q) % RZ(self.chii)))

    def wreath(self, v):
        """
        Get the wreath representation of v

        (v0,v1), toggle
        """
        v = vector(v)

        if v[0] % 2 == 0:
            return (self.A * v, self.A*v), False
        else:
            return (self.A * (v - self.r), self.A * (v + self.r)), True

    def run(self,v,w):
        """
        Treat v as a function acting on w
        """
        v = vector(v)

        def flip(x):
            if   x == "0": 
                return "1"
            elif x == "1": 
                return "0"
            else:
                raise ValueError("Can't flip {0}".format(x))

        try:
            hd, tl = w[0], w[1:]
            (v0,v1), t = self.wreath(v)

            if t:
                return flip(hd) + self.run((v0 if hd == "0" else v1), tl)
            else:
                return hd + self.run(v0, tl)
        except IndexError: # Happens if w is empty
            return "" 

    def wordCoord(self,w):
        """
        Return v such that v(0) = w
        """
        return sum([int(w[i]) * self.Ai^i * self.r for i in range(len(w))])

    def iterorbit(self,v,w):
        """
        Return an iterator for the orbit of v at w
        """
        y = self.run(v,w)
        while y != w:
            yield y
            y = self.run(v,y)

    def completeGraph(self):
        """
        Return a graph of all the vectors within the spectral radius of the matrix
        """
        spectralRadius = abs(max(self.chi.complex_roots()))
        searchRadius   = floor((spectralRadius^self.m) / (1 - spectralRadius))

        edges = {}

        def addEdges(v):
            (v0,v1), t = self.wreath(v)
            if tuple(v) in edges:
                return
            if t:
                edges[tuple(v)] = {tuple(v0): '0|1', tuple(v1): '1|0'}
            else:
                edges[tuple(v)] = {tuple(v0): '0|0', tuple(v1): '1|1'}
            addEdges(v0)
            addEdges(v1)
            
        def increment(v):
            try:
                hd, tl = v[0], v[1:]
                if hd != searchRadius:
                    return True, [(hd+1)] + tl
                else:
                    b, xs = increment(tl)
                    return b, [-searchRadius] + xs

            except IndexError:
                return False, []

        v = [-searchRadius] * self.m
        going = True
        while going:
            addEdges(v)
            going, v = increment(v)

        return ADiGraph(edges)

    def SCCs(self):
        """
        Return each terminal strongly connected component
        """

        complete = self.completeGraph()
        sccs = complete.strongly_connected_components()
        return sccs

    def anchorAt(self, v):
        """
        Return a graph representing the automaton anchored at v
        """
        edges = {}
        def getClosure(v):
            if tuple(v) in edges:
                return
            else:
                (v0,v1), t = self.wreath(v)
                if t:
                    edges[tuple(v)] = {tuple(v0): '0|1', tuple(v1): '1|0'}
                else:
                    edges[tuple(v)] = {tuple(v0): '0|0', tuple(v1): '1|1'}
                getClosure(v0)
                getClosure(v1)

        getClosure(v)
        return ADiGraph(edges)
