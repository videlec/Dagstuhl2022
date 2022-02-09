r"""
The shearing game

The main function is :func:`shearing_graph`
"""

from surface_dynamics import AbelianStratum, Origami
from surface_dynamics.misc.permutation import perm_cycles
import itertools

def perm_compose_cycle(u, c):
    r"""
    return u pre-composed with the cycle c
    """
    u0 = u[c[0]]
    for i in range(len(c)-1):
        u[c[i]] = u[c[i+1]]
    u[c[-1]] = u0

def neighbors(o):
    r"""
    This is the definition of neighbor where we do a single click in a single cylinder
    """
    r = o.r_tuple()
    u = o.u_tuple()
    for c in perm_cycles(r):
        uu = list(u)
        perm_compose_cycle(uu, c)
        oo = Origami(r, uu, as_tuple=True)
        yield oo

def shearing_graph(C, n):
    r"""
    INPUT:

    - ``C`` - component of stratum

    - ``n`` - number of squares

    EXAMPLES::

        sage: from surface_dynamics import AbelianStratum
        sage: C = AbelianStratum(2).unique_component()
        sage: print('|  n | size | diameter |')
        sage: print('|----|------|----------|')
        sage: for n in range(3, 10):
        ....:     G = shearing_graph(C, n)
        ....:     print('| {:>2} | {:>4} | {:>8} |'.format(n, G.num_verts(), G.diameter()))
    """
    # list of origamis with n squares in C
    origamis = set(C.origami_iterator(n, reduced=False, primitive=False))

    graphs = []
    while origamis:
        o = origamis.pop()
        orbit = {o: 0}
        todo = [o]
        G = DiGraph(multiedges=True, loops=True)
        while todo:
            o = todo.pop()
            oo = o.mirror()
            oo._set_standard_form()
            if oo not in orbit:
                orbit[oo] = len(orbit)
                todo.append(oo)
            G.add_edge(o, oo, 0)
            for oo in neighbors(o):
                assert oo.stratum_component() == C
                oo._set_standard_form()
                if oo not in orbit:
                    orbit[oo] = len(orbit)
                    todo.append(oo)
                G.add_edge(o, oo, 1)
        n0 = len(origamis)
        origamis.difference_update(orbit)
        assert n0 + 1 - len(origamis) == len(orbit)
        graphs.append(G)

    assert len(graphs) == 1 # conjecture
    return graphs[0]
