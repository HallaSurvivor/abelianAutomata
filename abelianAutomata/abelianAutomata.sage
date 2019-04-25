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

# {{{ miscelaneous (useful) things

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

def triLen(n):
    """
    All trit strings of length n
    """
    if n == 0:
        return [[]]
    else:
        recur = triLen(n-1)
        return map(lambda xs: [ 1]+xs, recur) + \
               map(lambda xs: [ 0]+xs, recur) + \
               map(lambda xs: [-1]+xs, recur)

def randVect(n=2,s=100):
    """
    A random integer vector of size n
    """
    return vector([int(s * random()) for _ in range(n)])

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

def borwein(n):
    """
    Return (ns,ps), where ns is a list of borwein polynomials of degree @n
    and ps is a list of borwein polynomials of degree < @n.

    Recall borwein polynomials are those polynomials in {-1,0,1}[x]
    """
    if n == 0:
        return ([1], [])
    else:
        (ns,ps) = borwein(n-1)
        ns2 = map(lambda p: x*p  , ns) + \
              map(lambda p: x*p+1, ns) + \
              map(lambda p: x*p-1, ns)
        return (ns2, ns+ps)
# }}}

#{{{ the automaton group class
class CompleteAutomaton(object):
    def __init__(self, A, e=None):
        """
        Construct the Complete Automaton for a given (@A, @e) pair.

        @e defaults to e1
        (so we default to the principal group)
        """
        if e == None: e = A.columns()[1]
        assert e[0] % 2 == 1

        self.A    = A
        self.e    = vector(e)

        self.m    = A.dimensions()[0]

        self.chi  = self.A.charpoly()

        self.Ai   = self.A.inverse()
        self.chii = RZ(self.Ai.charpoly())

        self.endo  = RZ.quo(self.chii)
        self.endo2 = ZZ.extension(self.chii,'a')

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

    def plot(self, f=None, plot=True):
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

        D = DiGraph(edges, loops=True)
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
                              ,vertex_size=10
                              ,figsize=[2*size,2*size]
                              ).plot()
        else:
            return D

#}}}

#{{{ testing conjectures, building intuition, etc
def leftResAlmostHom(trials=10, depth=1000, verbose=False):
    """
    Left residuation is `almost` a homomorphism in the following sense:
    (f+g)_0 != f_0 + g_0 iff f,g both odd

    Notice that when they ARE both odd, (f+g)_0 = f_0 + g_0 + gamma

    So then |(f+g)_0 - (f_0 + g_0)| <= |gamma|

    Given n, can we always find k such that
    |(f+g)_{0^k} - (f_{0^k} + g_{0^k})| <= 1/2^n

    (
        note: in principal machine this is clear since 
        EVERY function eventually residuates to I
    )

    ===========================================================

    looks true. Tried with 100 trials and depth 100,000 and every
    hom eventually had difference 0

    note though, this might be for ``dumb reasons'' because A^100,000
    is SUPER shrinking, and it's possible we just collapsed everything
    to the zero vector by the end
    """
    for m in matrices2 + matrices3 + matrices4 + matrices5:
        dim = m.dimensions()[0]
        aut = CompleteAutomaton(m,2*randVect(dim,100)+vector([1]*dim))
        print(aut)
        for i in range(trials):
            f = 2 * randVect(dim,1000) + vector([1]*dim) # make them odd
            g = 2 * randVect(dim,1000) + vector([1]*dim) # make them odd
            t = f + g

            fStart = f
            gStart = g

            k, min_ = 0, 1  # k such that the difference in del_0^k is min_
            for j in range(depth):
                f, _, _ = aut.wreath(f)
                g, _, _ = aut.wreath(g)
                t, _, _ = aut.wreath(t)
                n = aut.norm(t - (f+g))

                if n < min_:
                    k = j
                    min_ = n
                if min_ == 0:
                    break
            if verbose or min_ != 0:
                print("min: {} k: {} f: {} g: {} fk: {} gk: {}".format(min_,k,fStart,gStart,f,g))


def oddElementsOfPrincipalAreUnits():
    """
    It is true that every odd state of the principal is a unit,
    but is it true that every unit is an odd state of the principal?

    ======================================

    I ran this and started getting a BUNCH of counterexamples in dim3,
    but upon checking some by hand, it looks like they may not have been
    units after all?

    ======================================

    Nevermind, we're SOL... turns out these rings have infinitely 
    many units. Comes from dirichlet's unit thm + monogeneity conjecture
    """

    @cached_function
    def getVs(dim):
        if dim == 0:
            return [[]]
        else:
            olds = getVs(dim - 1)
            new = []
            for i in range(-50,50): # This is probably all of them...
                new += [old+[i] for old in olds]
            return new

    for m in matrices2 + matrices3: # + matrices4 + matrices5:
        aut = CompleteAutomaton(m)
        print(aut)
        units = [ tuple(v) for v in getVs(aut.m) if abs(aut.endo(v).norm()) == 1 ]
        D = aut.plot(plot=False)
        for u in units:
            if not D.has_vertex(u):
                print(u)


def tile(aut, n=None, save=None):
    """
    Returns a list representing the vectors in the nth approximation
    of the attractor of the IFS whose functions are given by 
    inverse residuation.

    @n is 10 by default (warning: gets slow quickly)
    """

    if n == None:
        n = 10 

    out = [vector([0,0,0])]
    old = out
    for i in range(n):
        print(i)
        new = []
        for v in old:
            new += [aut.A * (v + aut.e), aut.A * (v - aut.e), aut.A * v]
        out += new
        old = new

    g = Graphics()
    g += point3d(out, opacity=.6)
    g += point3d(aut.principal_vectors, color='red')
    if save:
        g.save(save)
    else:
        g.show()

def mahler(p):
    """
    Returns the mahler measure of the polynomial @p

    A polynomial with coefficients in {-1,0,1} with a root alpha
    exists if the ``mahler measure'' of alpha is < 2
    (Sufficient, not necessary).
    We take in the minimal polynomial of alpha @p.

    Interestingly, it seems all of our matrices have mahler measure = 2
    """
    return reduce(lambda x,y: x*y
                 ,map(lambda x: max(abs(x[0]), 1)
                     , p.roots(CC)
                     )
                 )

def plotBorwein(n, c=0):
    """
    Plot the solutions of b(x) = @c for all borwein polys b of degree <= @n
    """
    (deg_n, deg_lt_n) = borwein(n)
    polys = deg_n + deg_lt_n[:-1] # remove the constant 1 polynomial
    roots = [r[0] for rs in [(p-c).roots(CC) for p in polys] for r in rs]
    return points([CDF(r) for r in roots], aspect_ratio=1, transparent=True)

def plotMatrixRoots():
    """
    Plot the roots of chii for our matrices in red, except for x^m - 2 in orange
    """ 
    polys = [m.inverse().charpoly() for m in matrices2+matrices3+matrices4+matrices5]
    roots = [r[0] for rs in [p.roots(CC) for p in polys] for r in rs]
    
    badPolys = [x^n - 2 for n in range(2,6)]
    badRoots = [r[0] for rs in [p.roots(CC) for p in badPolys] for r in rs]
    return points([CDF(r) for r in roots], aspect_ratio=1, color='red', transparent=True) +\
           points([CDF(r) for r in badRoots], aspect_ratio=1, color='yellow', transparent=True)

def getRootsOfPathPoly():
    """
    Save the roots of each borwein polynomial b of degree <= 10 which are
    1 mod poly. Also save the roots of b-1 for comparison
    """
    for m in matrices2+matrices3+matrices4+matrices5:
        f = open("roots","a+")
        p = m.inverse().charpoly()
        aut = CompleteAutomaton(m)
        f.write("===============\n{0}\n===============\n\n".format(p))
        for b in B:
            if aut.endo(b) == 1 and b != 1:
                f.write("{0}\n".format(b))
                for x in b.roots(CC):
                    f.write("{0}\n".format(x[0]))
                f.write("\n{0}\n".format(b-1))
                for x in (b-1).roots(CC):
                    f.write("{0}\n".format(x[0]))
                f.write("\n\n")
        f.close()

def howCloseAreTheRootsExactly(m, n=10, plot=False):
    """
    Given a matrix @m, find the biggest distance between roots of
    b and b-1 with b = 1 mod m borwein of degree <= @n

    If @plot, then plot the roots of the maxDiff polynomial.
    Blue is b, Red is b-1
    """
    aut = CompleteAutomaton(m)
    maxDiff = 0
    maxPoly = None
    B = borwein(n)
    for b in B[0] + B[1]:
        if aut.endo(b) == 1 and b != 1:
            roots1 = [x[0] for x in b.roots(CC)]
            roots2 = [x[0] for x in (b-1).roots(CC)]
            if len(roots1) != len(roots2):
                print b
                print "Different number of roots!\n"
            else:
                diffs = [abs(r1 - r2) for r1 in roots1 for r2 in roots2]
                if maxDiff < max(diffs):
                    maxDiff = max(diffs)
                    maxPoly = b

    if plot:
        print maxDiff
        print maxPoly
        return points([CDF(r[0]) for r in maxPoly.roots(CC)], aspect_ratio=1, color='blue') +\
               points([CDF(r[0]) for r in (maxPoly-1).roots(CC)], aspect_ratio=1, color='red')

    else:
        return maxDiff, maxPoly

def howCloseCanTheRootsBe(m, n=10, plot=False):
    """
    Same as above, but look for the minimum of the maximums
    """
    aut = CompleteAutomaton(m)
    minDiff = 1000
    minPoly = None
    B = borwein(n)
    for b in B[0] + B[1]:
        if aut.endo(b) == 1 and b != 1:
            roots1 = [x[0] for x in b.roots(CC)]
            roots2 = [x[0] for x in (b-1).roots(CC)]
            if len(roots1) != len(roots2):
                print ""
                print b
                print "Different number of roots!\n"
            else:
                diffs = [abs(r1 - r2) for r1 in roots1 for r2 in roots2]
                if minDiff > max(diffs):
                    minDiff = max(diffs)
                    minPoly = b

    if plot:
        print minDiff
        print minPoly
        return points([CDF(r[0]) for r in minPoly.roots(CC)], aspect_ratio=1, color='blue') +\
               points([CDF(r[0]) for r in (minPoly-1).roots(CC)], aspect_ratio=1, color='red')

    else:
        return minDiff, minPoly
#}}}
