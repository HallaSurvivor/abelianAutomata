load("tim/graph_loop.sage")

#{{{ matrix definitions
matrices2 = \
[ companion_matrix([-(1/2) , 0      , 1], format='left')
, companion_matrix([1/2    , -1     , 1], format='left')
, companion_matrix([1/2    , -(1/2) , 1], format='left')
, companion_matrix([1/2    , 0      , 1], format='left')
, companion_matrix([1/2    , 1/2    , 1], format='left')
, companion_matrix([1/2    , 1      , 1], format='left')
]

matrices3 = \
[ companion_matrix([-(1/2) , -(1/2) , 1/2    , 1], format='left')
, companion_matrix([-(1/2) , 0      , 0      , 1], format='left')
, companion_matrix([-(1/2) , 0      , 1/2    , 1], format='left')
, companion_matrix([-(1/2) , 0      , 1      , 1], format='left')
, companion_matrix([-(1/2) , 1/2    , -(1/2) , 1], format='left')
, companion_matrix([-(1/2) , 1/2    , 0      , 1], format='left')
, companion_matrix([-(1/2) , 1      , -1     , 1], format='left')
, companion_matrix([1/2    , -(1/2) , -(1/2) , 1], format='left')
, companion_matrix([1/2    , 0      , -1     , 1], format='left')
, companion_matrix([1/2    , 0      , -(1/2) , 1], format='left')
, companion_matrix([1/2    , 0      , 0      , 1], format='left')
, companion_matrix([1/2    , 1/2    , 0      , 1], format='left')
, companion_matrix([1/2    , 1/ 2   , 1/2    , 1], format='left')
, companion_matrix([1/2    , 1      , 1      , 1], format='left')
]

matrices4 = \
[ companion_matrix([-(1/2) , -1     , 0      , 1      , 1], format='left')
, companion_matrix([-(1/2) , -(1/2) , 0      , 1/2    , 1], format='left')
, companion_matrix([-(1/2) , 0      , 0      , 0      , 1], format='left')
, companion_matrix([-(1/2) , 0      , 1/2    , -(1/2) , 1], format='left')
, companion_matrix([-(1/2) , 0      , 1/2    , 1/2    , 1], format='left')
, companion_matrix([-(1/2) , 1/2    , 0      , -(1/2) , 1], format='left')
, companion_matrix([-(1/2) , 1      , 0      , -1     , 1], format='left')
, companion_matrix([1/2    , -1     , 1      , -1     , 1], format='left')
, companion_matrix([1/ 2   , -1     , 3/2    , -(3/2) , 1], format='left')
, companion_matrix([1/2    , -(1/2) , -(1/2) , 0      , 1], format='left')
, companion_matrix([1/2    , -(1/2) , 0      , -(1/2) , 1], format='left')
, companion_matrix([1/2    , -(1/2) , 0      , 0      , 1], format='left')
, companion_matrix([1/2    , -(1/2) , 1/2    , -1     , 1], format='left')
, companion_matrix([1/2    , -(1/2) , 1/2    , -(1/2) , 1], format='left')
, companion_matrix([1/2    , -(1/2) , 1      , -1     , 1], format='left')
, companion_matrix([1/ 2   , -(1/2) , 1      , -(1/2) , 1], format='left')
, companion_matrix([1/2    , 0      , -1     , 0      , 1], format='left')
, companion_matrix([1/2    , 0      , -(1/2) , -(1/2) , 1], format='left')
, companion_matrix([1/2    , 0      , -(1/2) , 0      , 1], format='left')
, companion_matrix([1/2    , 0      , -(1/2) , 1/2    , 1], format='left')
, companion_matrix([1/2    , 0      , 0      , -(1/2) , 1], format='left')
, companion_matrix([1/2    , 0      , 0      , 0      , 1], format='left')
, companion_matrix([1/2    , 0      , 0      , 1/2    , 1], format='left')
, companion_matrix([1/2    , 0      , 1/2    , -(1/2) , 1], format='left')
, companion_matrix([1/2    , 0      , 1/2    , 0      , 1], format='left')
, companion_matrix([1/2    , 0      , 1/2    , 1/ 2   , 1], format='left')
, companion_matrix([1/2    , 0      , 1      , 0      , 1], format='left')
, companion_matrix([1/2    , 1/2    , -(1/2) , 0      , 1], format='left')
, companion_matrix([1/2    , 1/2    , 0      , 0      , 1], format='left')
, companion_matrix([1/2    , 1/2    , 0      , 1/2    , 1], format='left')
, companion_matrix([1/2    , 1/2    , 1/2    , 1/2    , 1], format='left')
, companion_matrix([1/2    , 1/2    , 1/2    , 1      , 1], format='left')
, companion_matrix([1/2    , 1/2    , 1      , 1/2    , 1], format='left')
, companion_matrix([1/2    , 1/2    , 1      , 1      , 1], format='left')
, companion_matrix([1/2    , 1      , 1      , 1      , 1], format='left')
, companion_matrix([1/2    , 1      , 3/2    , 3/2    , 1], format='left')
]

matrices5 = \
[ companion_matrix([-(1/2) , -1     , -(1/2) , 1/2    , 1      , 1], format='left')
, companion_matrix([-(1/2) , -(1/2) , -(1/2) , 0      , 1      , 1], format='left')
, companion_matrix([-(1/2) , -(1/2) , -(1/2) , 1/2    , 1/2    , 1], format='left')
, companion_matrix([-(1/2) , -(1/2) , -(1/2) , 1/2    , 1      , 1], format='left')
, companion_matrix([-(1/2) , -(1/2) , 0      , 0      , 1/2    , 1], format='left')
, companion_matrix([-(1/2) , -(1/2) , 0      , 1/2    , 1      , 1], format='left')
, companion_matrix([-(1/2) , -(1/2) , 0      , 1      , 3/2    , 1], format='left')
, companion_matrix([-(1/2) , 0      , -(1/2) , 1/2    , 0      , 1], format='left')
, companion_matrix([-(1/2) , 0      , 0      , -(1/2) , 1/2    , 1], format='left')
, companion_matrix([-(1/2) , 0      , 0      , 0      , 0      , 1], format='left')
, companion_matrix([-(1/2) , 0      , 0      , 0      , 1/2    , 1], format='left')
, companion_matrix([-(1/2) , 0      , 0      , 1/2    , 0      , 1], format='left')
, companion_matrix([-(1/2) , 0      , 0      , 1/2    , 1/2    , 1], format='left')
, companion_matrix([-(1/2) , 0      , 1/2    , -(1/2) , 0      , 1], format='left')
, companion_matrix([-(1/2) , 0      , 1/2    , 0      , 0      , 1], format='left')
, companion_matrix([-(1/2) , 0      , 1/2    , 0      , 1/2    , 1], format='left')
, companion_matrix([-(1/2) , 0      , 1/2    , 1/2    , 1      , 1], format='left')
, companion_matrix([-(1/2) , 1/2    , -1     , 1      , -(1/2) , 1], format='left')
, companion_matrix([-(1/2) , 1/2    , -(1/2) , 1/ 2   , -(1/2) , 1], format='left')
, companion_matrix([-(1/2) , 1/2    , -(1/2) , 1      , -(1/2) , 1], format='left')
, companion_matrix([-(1/2) , 1/2    , 0      , 0      , -(1/2) , 1], format='left')
, companion_matrix([-(1/2) , 1/2    , 0      , 0      , 0      , 1], format='left')
, companion_matrix([-(1/2) , 1/2    , 0      , 1/2    , 0      , 1], format='left')
, companion_matrix([-(1/2) , 1/2    , 1/2    , -1     , 0      , 1], format='left')
, companion_matrix([-(1/2) , 1/2    , 1/2    , -(1/2) , -(1/2) , 1], format='left')
, companion_matrix([-(1/2) , 1/2    , 1/2    , -(1/2) , 0      , 1], format='left')
, companion_matrix([-(1/2) , 1/2    , 1      , -1     , -(1/2) , 1], format='left')
, companion_matrix([-(1/2) , 1      , -1     , 1      , -1     , 1], format='left')
, companion_matrix([-(1/2) , 1      , -(1/2) , 1/2    , -1     , 1], format='left')
, companion_matrix([1/ 2   , -1     , 1/2    , 1/2    , -1     , 1], format='left')
, companion_matrix([1/2    , -(1/2) , 0      , 0      , -(1/2) , 1], format='left')
, companion_matrix([1/ 2   , -(1/2) , 0      , 1/2    , -1     , 1], format='left')
, companion_matrix([1/2    , -(1/2) , 0      , 1      , -(3/2) , 1], format='left')
, companion_matrix([1/ 2   , -(1/2) , 1/2    , 0      , -1     , 1], format='left')
, companion_matrix([1/2    , -(1/2) , 1/2    , 1/2    , -1     , 1], format='left')
, companion_matrix([1/ 2   , -(1/2) , 1/2    , 1/2    , -(1/2) , 1], format='left')
, companion_matrix([1/2    , 0      , -(1/2) , -(1/2) , 0      , 1], format='left')
, companion_matrix([1/ 2   , 0      , -(1/2) , 0      , -(1/2) , 1], format='left')
, companion_matrix([1/2    , 0      , -(1/2) , 0      , 0      , 1], format='left')
, companion_matrix([1/2    , 0      , -(1/2) , 1/2    , -1     , 1], format='left')
, companion_matrix([1/2    , 0      , 0      , -(1/2) , -(1/2) , 1], format='left')
, companion_matrix([1/2    , 0      , 0      , 0      , -(1/2) , 1], format='left')
, companion_matrix([1/2    , 0      , 0      , 0      , 0      , 1], format='left')
, companion_matrix([1/2    , 0      , 0      , 1/2    , -(1/2) , 1], format='left')
, companion_matrix([1/2    , 0      , 0      , 1/2    , 0      , 1], format='left')
, companion_matrix([1/2    , 0      , 1/2    , 1/2    , 0      , 1], format='left')
, companion_matrix([1/2    , 1/ 2   , -1     , -1     , 1/2    , 1], format='left')
, companion_matrix([1/2    , 1/2    , -(1/2) , -1     , 0      , 1], format='left')
, companion_matrix([1/2    , 1/ 2   , -(1/2) , -(1/2) , 0      , 1], format='left')
, companion_matrix([1/2    , 1/2    , -(1/2) , -(1/2) , 1/2    , 1], format='left')
, companion_matrix([1/2    , 1/2    , 0      , 0      , 0      , 1], format='left')
, companion_matrix([1/2    , 1/2    , 0      , 0      , 1/2    , 1], format='left')
, companion_matrix([1/2    , 1/2    , 0      , 1/2    , 0      , 1], format='left')
, companion_matrix([1/2    , 1/2    , 1/2    , 1/2    , 1/2    , 1], format='left')
, companion_matrix([1/2    , 1/2    , 1/2    , 1      , 1/2    , 1], format='left')
, companion_matrix([1/2    , 1/2    , 1      , 1      , 1/2    , 1], format='left')
, companion_matrix([1/2    , 1      , 1/2    , 1/2    , 1      , 1], format='left')
, companion_matrix([1/2    , 1      , 1      , 1      , 1      , 1], format='left')
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
# }}}


class CompleteAutomaton(object):
    def __init__(self, A, e=None):
        """
        Construct the Complete Automaton for a given (@A, @e) pair.

        @e defaults to e1
        (so we default to the principal group)
        """
        if e == None: e = A.columns()[1]

        self.A    = A
        self.e    = vector(e)

        self.m    = A.dimensions()[0]

        self.chi  = self.A.charpoly()

        self.Ai   = self.A.inverse()
        self.chii = RZ(self.Ai.charpoly())

        self.endo = RZ.quo(self.chii)

    def __repr__(self):
        return \
"""
===============================
CompleteAutomaton with matrix:
{}
and residuation vector:
{}
===============================
""".format(self.A,self.e)

    def scaleByPoly(self, p):
        """
        Returns a new CompleteAutomaton: self scaled by @p
        """
        p = self.endo(list(p))
        q = self.endo(list(self.e))

        return CompleteAutomaton(self.A, p * q)

    def wreath(self, f):
        """
        Get the wreath representation of a function @f

        Returns: 
            f0, f1, toggle
        """
        f = vector(f)

        if f[0] % 2 == 0:
            return self.A * f, self.A * f, False
        else:
            return self.A * (f - self.e), self.A * (f + self.e), True

    def run(self,f,u):
        """
        Run @f on a string @u

        Returns:
            f u, del_u f
        """
        f = vector(f)

        if u == "":
            return "", f
        else:
            uNew,  fNew     = self.run(f,u[:-1])
            fNew0, fNew1, t = self.wreath(fNew)
            
            if t:
                if u[-1] == "0": return uNew + "1", fNew0
                if u[-1] == "1": return uNew + "0", fNew1
            else:
                return uNew + u[-1], fNew0

    def norm(self,f):
        """
        Get the norm of a function @f

        Returns: 
            1/2^k where k is the position of the first 1 in f(0^omega)
        """
        f = vector(f)

        if f.is_zero(): 
            return 0
        else:
            f0, _, t = self.wreath(f)
            if t: return 1
            else: return 1/2 * self.norm(f0)

    def wordCoord(self,u,v=None):
        """
        Returns: 
            (f,p) such that f(0^omega) = uv* in (p.G) (G is the current group)
        """
        if v == None: v = ""

        pu = self.endo([int(ui) for ui in u])
        pv = self.endo([int(vi) for vi in v])

        x = self.endo.0

        quo = 1 - x^len(v) if v else self.endo(1)
        f = vector(pu*quo + pv*(x^len(u)))
        return f, quo


    def iterorbit(self,f,u):
        """
        Returns an iterator for the orbit of @f at @w
        """
        y = self.run(f,u)
        while y != u:
            yield y
            y = self.run(f,y)

    def anchorAt(self, f=None, plot=True):
        """
        Plots a graph representing the automaton anchored at @f

        if @plot is False, then we return the graph itself

        @f defaults to e1
        """
        if f == None:
            f = self.A.columns()[1] # The worst possible way to say e1

        edges = {}
        def getClosure(v):
            if tuple(v) in edges:
                return
            else:
                v0, v1, t = self.wreath(v)
                if t:
                    edges[tuple(v)] = {tuple(v0): '0|1', tuple(v1): '1|0'}
                else:
                    edges[tuple(v)] = {tuple(v0): '0|0', tuple(v1): '1|1'}
                getClosure(v0)
                getClosure(v1)

        getClosure(f)

        D = DiGraph(edges)
        if plot:
            vertex_colormap = {"red": [tuple(f)]}
            edge_colormap = {"0|0": "grey", "0|1":"green", "1|0":"blue"}

            size = max((D.order() + 1)/2, 10)

            return D.graphplot(layout='spring'
                              ,iterations=1000
                              ,dpi=200
                              ,vertex_labels=False
                              ,color_by_label=edge_colormap
                              ,vertex_colors=vertex_colormap
                              ,vertex_color='white'
                              ,figsize=[size, size]
                              ).plot()
        else:
            return D
