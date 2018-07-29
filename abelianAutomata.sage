load("tim/graph_loop.sage")

#{{{ matrix definitions
matrices2 = 
[ companion_matrix([-(1/2) , 0      , 1], align='left')
, companion_matrix([1/2    , -1     , 1], align='left')
, companion_matrix([1/2    , -(1/2) , 1], align='left')
, companion_matrix([1/2    , 0      , 1], align='left')
, companion_matrix([1/2    , 1/2    , 1], align='left')
, companion_matrix([1/2    , 1      , 1], align='left')
]

matrices3 = 
[ companion_matrix([-(1/2) , -(1/2) , 1/2    , 1], align='left')
, companion_matrix([-(1/2) , 0      , 0      , 1], align='left')
, companion_matrix([-(1/2) , 0      , 1/2    , 1], align='left')
, companion_matrix([-(1/2) , 0      , 1      , 1], align='left')
, companion_matrix([-(1/2) , 1/2    , -(1/2) , 1], align='left')
, companion_matrix([-(1/2) , 1/2    , 0      , 1], align='left')
, companion_matrix([-(1/2) , 1      , -1     , 1], align='left')
, companion_matrix([1/2    , -(1/2) , -(1/2) , 1], align='left')
, companion_matrix([1/2    , 0      , -1     , 1], align='left')
, companion_matrix([1/2    , 0      , -(1/2) , 1], align='left')
, companion_matrix([1/2    , 0      , 0      , 1], align='left')
, companion_matrix([1/2    , 1/2    , 0      , 1], align='left')
, companion_matrix([1/2    , 1/ 2   , 1/2    , 1], align='left')
, companion_matrix([1/2    , 1      , 1      , 1], align='left')
]

matrices4 = 
[ companion_matrix([-(1/2) , -1     , 0      , 1      , 1], align='left')
, companion_matrix([-(1/2) , -(1/2) , 0      , 1/2    , 1], align='left')
, companion_matrix([-(1/2) , 0      , 0      , 0      , 1], align='left')
, companion_matrix([-(1/2) , 0      , 1/2    , -(1/2) , 1], align='left')
, companion_matrix([-(1/2) , 0      , 1/2    , 1/2    , 1], align='left')
, companion_matrix([-(1/2) , 1/2    , 0      , -(1/2) , 1], align='left')
, companion_matrix([-(1/2) , 1      , 0      , -1     , 1], align='left')
, companion_matrix([1/2    , -1     , 1      , -1     , 1], align='left')
, companion_matrix([1/ 2   , -1     , 3/2    , -(3/2) , 1], align='left')
, companion_matrix([1/2    , -(1/2) , -(1/2) , 0      , 1], align='left')
, companion_matrix([1/2    , -(1/2) , 0      , -(1/2) , 1], align='left')
, companion_matrix([1/2    , -(1/2) , 0      , 0      , 1], align='left')
, companion_matrix([1/2    , -(1/2) , 1/2    , -1     , 1], align='left')
, companion_matrix([1/2    , -(1/2) , 1/2    , -(1/2) , 1], align='left')
, companion_matrix([1/2    , -(1/2) , 1      , -1     , 1], align='left')
, companion_matrix([1/ 2   , -(1/2) , 1      , -(1/2) , 1], align='left')
, companion_matrix([1/2    , 0      , -1     , 0      , 1], align='left')
, companion_matrix([1/2    , 0      , -(1/2) , -(1/2) , 1], align='left')
, companion_matrix([1/2    , 0      , -(1/2) , 0      , 1], align='left')
, companion_matrix([1/2    , 0      , -(1/2) , 1/2    , 1], align='left')
, companion_matrix([1/2    , 0      , 0      , -(1/2) , 1], align='left')
, companion_matrix([1/2    , 0      , 0      , 0      , 1], align='left')
, companion_matrix([1/2    , 0      , 0      , 1/2    , 1], align='left')
, companion_matrix([1/2    , 0      , 1/2    , -(1/2) , 1], align='left')
, companion_matrix([1/2    , 0      , 1/2    , 0      , 1], align='left')
, companion_matrix([1/2    , 0      , 1/2    , 1/ 2   , 1], align='left')
, companion_matrix([1/2    , 0      , 1      , 0      , 1], align='left')
, companion_matrix([1/2    , 1/2    , -(1/2) , 0      , 1], align='left')
, companion_matrix([1/2    , 1/2    , 0      , 0      , 1], align='left')
, companion_matrix([1/2    , 1/2    , 0      , 1/2    , 1], align='left')
, companion_matrix([1/2    , 1/2    , 1/2    , 1/2    , 1], align='left')
, companion_matrix([1/2    , 1/2    , 1/2    , 1      , 1], align='left')
, companion_matrix([1/2    , 1/2    , 1      , 1/2    , 1], align='left')
, companion_matrix([1/2    , 1/2    , 1      , 1      , 1], align='left')
, companion_matrix([1/2    , 1      , 1      , 1      , 1], align='left')
, companion_matrix([1/2    , 1      , 3/2    , 3/2    , 1], align='left')
]

matrices5 = 
[ companion_matrix([-(1/2) , -1     , -(1/2) , 1/2    , 1      , 1], align='left')
, companion_matrix([-(1/2) , -(1/2) , -(1/2) , 0      , 1      , 1], align='left')
, companion_matrix([-(1/2) , -(1/2) , -(1/2) , 1/2    , 1/2    , 1], align='left')
, companion_matrix([-(1/2) , -(1/2) , -(1/2) , 1/2    , 1      , 1], align='left')
, companion_matrix([-(1/2) , -(1/2) , 0      , 0      , 1/2    , 1], align='left')
, companion_matrix([-(1/2) , -(1/2) , 0      , 1/2    , 1      , 1], align='left')
, companion_matrix([-(1/2) , -(1/2) , 0      , 1      , 3/2    , 1], align='left')
, companion_matrix([-(1/2) , 0      , -(1/2) , 1/2    , 0      , 1], align='left')
, companion_matrix([-(1/2) , 0      , 0      , -(1/2) , 1/2    , 1], align='left')
, companion_matrix([-(1/2) , 0      , 0      , 0      , 0      , 1], align='left')
, companion_matrix([-(1/2) , 0      , 0      , 0      , 1/2    , 1], align='left')
, companion_matrix([-(1/2) , 0      , 0      , 1/2    , 0      , 1], align='left')
, companion_matrix([-(1/2) , 0      , 0      , 1/2    , 1/2    , 1], align='left')
, companion_matrix([-(1/2) , 0      , 1/2    , -(1/2) , 0      , 1], align='left')
, companion_matrix([-(1/2) , 0      , 1/2    , 0      , 0      , 1], align='left')
, companion_matrix([-(1/2) , 0      , 1/2    , 0      , 1/2    , 1], align='left')
, companion_matrix([-(1/2) , 0      , 1/2    , 1/2    , 1      , 1], align='left')
, companion_matrix([-(1/2) , 1/2    , -1     , 1      , -(1/2) , 1], align='left')
, companion_matrix([-(1/2) , 1/2    , -(1/2) , 1/ 2   , -(1/2) , 1], align='left')
, companion_matrix([-(1/2) , 1/2    , -(1/2) , 1      , -(1/2) , 1], align='left')
, companion_matrix([-(1/2) , 1/2    , 0      , 0      , -(1/2) , 1], align='left')
, companion_matrix([-(1/2) , 1/2    , 0      , 0      , 0      , 1], align='left')
, companion_matrix([-(1/2) , 1/2    , 0      , 1/2    , 0      , 1], align='left')
, companion_matrix([-(1/2) , 1/2    , 1/2    , -1     , 0      , 1], align='left')
, companion_matrix([-(1/2) , 1/2    , 1/2    , -(1/2) , -(1/2) , 1], align='left')
, companion_matrix([-(1/2) , 1/2    , 1/2    , -(1/2) , 0      , 1], align='left')
, companion_matrix([-(1/2) , 1/2    , 1      , -1     , -(1/2) , 1], align='left')
, companion_matrix([-(1/2) , 1      , -1     , 1      , -1     , 1], align='left')
, companion_matrix([-(1/2) , 1      , -(1/2) , 1/2    , -1     , 1], align='left')
, companion_matrix([1/ 2   , -1     , 1/2    , 1/2    , -1     , 1], align='left')
, companion_matrix([1/2    , -(1/2) , 0      , 0      , -(1/2) , 1], align='left')
, companion_matrix([1/ 2   , -(1/2) , 0      , 1/2    , -1     , 1], align='left')
, companion_matrix([1/2    , -(1/2) , 0      , 1      , -(3/2) , 1], align='left')
, companion_matrix([1/ 2   , -(1/2) , 1/2    , 0      , -1     , 1], align='left')
, companion_matrix([1/2    , -(1/2) , 1/2    , 1/2    , -1     , 1], align='left')
, companion_matrix([1/ 2   , -(1/2) , 1/2    , 1/2    , -(1/2) , 1], align='left')
, companion_matrix([1/2    , 0      , -(1/2) , -(1/2) , 0      , 1], align='left')
, companion_matrix([1/ 2   , 0      , -(1/2) , 0      , -(1/2) , 1], align='left')
, companion_matrix([1/2    , 0      , -(1/2) , 0      , 0      , 1], align='left')
, companion_matrix([1/2    , 0      , -(1/2) , 1/2    , -1     , 1], align='left')
, companion_matrix([1/2    , 0      , 0      , -(1/2) , -(1/2) , 1], align='left')
, companion_matrix([1/2    , 0      , 0      , 0      , -(1/2) , 1], align='left')
, companion_matrix([1/2    , 0      , 0      , 0      , 0      , 1], align='left')
, companion_matrix([1/2    , 0      , 0      , 1/2    , -(1/2) , 1], align='left')
, companion_matrix([1/2    , 0      , 0      , 1/2    , 0      , 1], align='left')
, companion_matrix([1/2    , 0      , 1/2    , 1/2    , 0      , 1], align='left')
, companion_matrix([1/2    , 1/ 2   , -1     , -1     , 1/2    , 1], align='left')
, companion_matrix([1/2    , 1/2    , -(1/2) , -1     , 0      , 1], align='left')
, companion_matrix([1/2    , 1/ 2   , -(1/2) , -(1/2) , 0      , 1], align='left')
, companion_matrix([1/2    , 1/2    , -(1/2) , -(1/2) , 1/2    , 1], align='left')
, companion_matrix([1/2    , 1/2    , 0      , 0      , 0      , 1], align='left')
, companion_matrix([1/2    , 1/2    , 0      , 0      , 1/2    , 1], align='left')
, companion_matrix([1/2    , 1/2    , 0      , 1/2    , 0      , 1], align='left')
, companion_matrix([1/2    , 1/2    , 1/2    , 1/2    , 1/2    , 1], align='left')
, companion_matrix([1/2    , 1/2    , 1/2    , 1      , 1/2    , 1], align='left')
, companion_matrix([1/2    , 1/2    , 1      , 1      , 1/2    , 1], align='left')
, companion_matrix([1/2    , 1      , 1/2    , 1/2    , 1      , 1], align='left')
, companion_matrix([1/2    , 1      , 1      , 1      , 1      , 1], align='left')
]

# }}}

# {{{ Miscelaneous (useful) things

# Define polynomial rings
RZ.<z> = PolynomialRing(ZZ)
RQ.<x> = PolynomialRing(QQ)

def binLen(n):
    """
    All binary strings of length n
    """
    if n == 0:
        return [""]
    else:
        recur = binLen(n-1)
        return map(lambda xs: "1"+xs, recur) + map(lambda xs: "0"+xs, recur)

def expandTransPeriod(self,u,v):
    """
    Return a generator uv* given u and v
    """
    for c in u:
        yield c

    if v == "":
        return

    i = 0
    while v != "":
        yield v[i % len(v)]
        i += 1


# plot settings
PLOT_DPI = 200
PLOT_ITERS = 1000

MAX_DEPTH=10000
sys.setrecursionlimit(MAX_DEPTH)

class ADiGraph(DiGraph):
    def plot2(self, labeled=False, **kwargs):
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
                             ,vertex_labels  = labeled
                             ,**kwargs
                             ).plot()


# }}}


class CompleteAutomaton(object):
    def __init__(self, A, r=None):
        """
        Construct the Complete Automaton for a given (A,r) pair.

        Principal (r = e1) by default
        """
        self.A    = A
        self.r    = vector(r)

        self.m    = len(self.r)
        self.chi  = self.A.charpoly()

        self.Ai   = self.A.inverse()
        self.chii = RZ(self.Ai.charpoly())

        self.endo = RQ.quo(self.chii)
        self.x    = self.endo.0  # the variable for self.R

    def scaleByPoly(self, p):
        """
        Return a new CompleteAutomaton which is this one scaled by p
        """
        print("WARNING: busted")
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
            return (self.A * v, self.A * v), False
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

    def wordCoord(self,u,v=None):
        """
        Return (f,p) such that f(0^omega) = uv* in (p.G) (G is the current group)
        """
        pu = self.endo([int(ui) for ui in u])
        pv = self.endo([int(vi) for vi in v]) if v else self.endo(0)

        quo = self.endo(1-self.x^len(v)) if v else self.endo(1)
        f = vector(list(self.endo(pu*quo + pv*(self.x^len(u)))))
        return f, quo


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

    def knfApprox(self, u, v, f=None):
        """
        Return a generator of approximations for omegaKNF 
        with coefficients given by uv* starting at f (defaults to delta)
        """
        if f == None:
            f = self.r
        else:
            f = vector(f)

        fn = vector([0]*self.m)  # start at identity
        for c in expandTransPeriod(u,v):
            fn += int(c) * f
            f = self.Ai * f
            yield fn

    def knfExact(self, u, v, f=None):
        """
        Return a (vector, polynomial) pair (x, p) such that
        x \in p . G (where G is the current group) is the
        function corresponding to the limit of the omegaKNF 
        uv* starting at f (delta by default)
        """
        if f == None:
            f = self.r
        else:
            f = vector(f)
