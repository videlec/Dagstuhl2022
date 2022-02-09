r"""
The layout problem for square-tiled surfaces

Square-tiled surfaces are encoded as pairs ``(r, u)`` where both ``r`` and
``u`` are lists of integers corresponding to the permutations ``i -> r[i]``
and ``i -> u[i]``.
"""

from sage.all import MixedIntegerLinearProgram, Graph, DiGraph
from sage.numerical.mip import MIPSolverException

def origami_digraph(r, u, loops=True, multiedges=True):
    n = len(r)
    G = Digraph(loops=loops, multiedges=multiedges)
    for i in range(n):
        if loops or r[i] != i:
            G.add_edge(i, r[i], 'r')
        if loops or u[i] != i:
            G.add_edge(i, u[i], 'u')
    return G


def staircase_layout(r, u):
    r"""
    Return a staircase layout if any
    """
    G = origami_digraph(r, u, loops=False, multiedges=False)
    H = G.hamiltonian_path(maximize=True)
    return H if (H is not None and len(H) == len(G)) else None


def soft_layout(r, u, cut_first=False, all_overlaps=False, verbose=False):
    r"""
    Return a soft layout if any

    This function uses integer linear programming and it is strongly advised to
    use a decent solver. The only one installed by default in sage (GLPK) turns
    out to be very slow.
    """
    n = len(r)
    r_inv = [None] * n
    u_inv = [None] * n
    for i in range(n):
        r_inv[r[i]] = i
        u_inv[u[i]] = i

    # Variables
    M = MixedIntegerLinearProgram()
    x = M.new_variable(integer=True)    # x coordinate of each square
    xdiff = M.new_variable(binary=True) # sign of each difference x[i] - x[j]
    y = M.new_variable(integer=True)    # y coordinate of each square
    ydiff = M.new_variable(binary=True) # sign of each difference y[i] - y[j]
    r_cut = M.new_variable(binary=True) # whether (i, r[i]) is cut
    u_cut = M.new_variable(binary=True) # whether (i, u[i]) is cut

    # Constraints
    M.add_constraint(x[0] == 0)
    M.add_constraint(y[0] == 0)
    for i in range(1, n):
        # absolute bound on our window
        M.add_constraint(x[i] >= -n)
        M.add_constraint(x[i] <= n)
        M.add_constraint(y[i] >= -n)
        M.add_constraint(y[i] <= n)

    for i in range(n):
        # if not r-cut then x-difference is one and y-difference is zero
        M.add_constraint(- 3 * n * r_cut[i] <= x[r[i]] - x[i] - 1)
        M.add_constraint(x[r[i]] - x[i] - 1 <= 3 * n * r_cut[i])
        M.add_constraint(- 3 * n * r_cut[i] <= y[r[i]] - y[i])
        M.add_constraint(y[r[i]] - y[i] <= 3 * n * r_cut[i])

        # if no u-cut then y-difference is one and x-difference is zero
        M.add_constraint(- 3 * n * u_cut[i] <= y[u[i]] - y[i] - 1)
        M.add_constraint(y[u[i]] - y[i] - 1 <= 3 * n * u_cut[i])
        M.add_constraint(- 3 * n * u_cut[i] <= x[u[i]] - x[i])
        M.add_constraint(x[u[i]] - x[i] <= 3 * n * u_cut[i])

    # leave some room for a spanning tree
    M.add_constraint(M.sum(r_cut[i] for i in range(n)) + M.sum(u_cut[i] for i in range(n)) <= n + 1)

    while True:
        try:
            M.solve()
        except MIPSolverException:
            return

        vx = M.get_values(x)
        vx = [int(vx[i]) for i in range(n)]
        vy = M.get_values(y)
        vy = [int(vy[i]) for i in range(n)]

        # remove a cut if any
        vr_cut = M.get_values(r_cut)
        vu_cut = M.get_values(u_cut)
        
        G = Graph(n, multiedges=False, loops=False)
        for i in range(n):
            if not vr_cut[i] and i != r[i]:
                G.add_edge(i, r[i], 'r')
            if not vu_cut[i] and i != u[i]:
                G.add_edge(i, u[i], 'u')
        comps = G.connected_components()
        if len(comps) > 1:
            comp = set(comps[0])
            if verbose:
                print('forbid cut for {}'.format(comp))
            rights = [i for i in comp if r[i] not in comp]
            lefts = [r_inv[i] for i in comp if r_inv[i] not in comp]
            ups = [i for i in comp if u[i] not in comp]
            downs = [u_inv[i] for i in comp if u_inv[i] not in comp]
            cut = M.sum(r_cut[i] for i in rights + lefts) + M.sum(u_cut[i] for i in ups + downs)
            M.add_constraint(cut <= len(rights) + len(lefts) + len(ups) + len(downs) - 1)
            if cut_first:
                continue

        # forbid squares at the same position if any
        # (we avoid doing that too much since this potentially creates 2 n^2 new variables)
        has_overlap = False
        for i in range(n):
            for j in range(i):
                if vx[i] == vx[j] and vy[i] == vy[j]:
                    if verbose:
                        print('forbid overlap of squares {} and {}'.format(i, j))
                    M.add_constraint(x[i] - x[j] + 3 * n * xdiff[i,j] + y[i] - y[j] + 3 * n * ydiff[i,j] >= 1)
                    M.add_constraint(x[j] - x[i] + 3 * n * (1 - xdiff[i,j]) + y[i] - y[j] + 3 * n * ydiff[i,j] >= 1)
                    M.add_constraint(x[i] - x[j] + 3 * n * xdiff[i,j] + y[j] - y[i] + 3 * n * (1 - ydiff[i,j]) >= 1)
                    M.add_constraint(x[j] - x[i] + 3 * n * (1 - xdiff[i,j]) + y[j] - y[i] + 3 * n * (1 - ydiff[i,j]) >= 1)
                    has_overlap = True
                    if not all_overlaps:
                        break
            if has_overlap and not all_overlaps:
                break

        if len(comps) == 1 and not has_overlap:
            G = DiGraph(n, multiedges=False, loops=False)
            for i in range(n):
                if not vr_cut[i] and i != r[i]:
                    G.add_edge(i, r[i], 'r')
                if not vu_cut[i] and i != u[i]:
                    G.add_edge(i, u[i], 'u')
            return G

############
# Plotting #
############

def positions_from_layout(r, u, G):
    assert G.is_connected()
    n = len(r)
    x_pos = [None] * n
    y_pos = [None] * n
    x_pos[0] = y_pos[0] = 0
    todo = [0]
    while todo:
        i = todo.pop()
        for (_, j, lab) in G.outgoing_edges(i, labels=True):
            if lab == 'r':
                assert r[i] == j
            elif lab == 'u':
                assert u[i] == j
            else:
                raise RuntimeError

            if x_pos[j] is None:
                assert y_pos[j] is None
                if lab == 'r':
                    x_pos[j] = x_pos[i] + 1
                    y_pos[j] = y_pos[i]
                else:
                    x_pos[j] = x_pos[i]
                    y_pos[j] = y_pos[i] + 1
                todo.append(j)
        for (j, _, lab) in G.incoming_edges(i, labels=True):
            if lab == 'r':
                assert r[j] == i
            elif lab == 'u':
                assert u[j] == i
            else:
                raise RuntimeError

            if x_pos[j] is None:
                assert y_pos[j] is None
                if lab == 'r':
                    x_pos[j] = x_pos[i] - 1
                    y_pos[j] = y_pos[i]
                else:
                    x_pos[j] = x_pos[i]
                    y_pos[j] = y_pos[i] - 1
                todo.append(j)

    return x_pos, y_pos


def square_tiled_plot(r, u, G):
    r"""
    Return a plot of the origami ``(r, u)`` using the layout determined by the subgraph ``G``.
    """
    n = len(r)
    r_inv = [None] * n
    u_inv = [None] * n
    for i in range(n):
        r_inv[r[i]] = i
        u_inv[u[i]] = i

    x, y = positions_from_layout(r, u, G)

    Gr = Graphics()
    for i in range(n):
        Gr += polygon2d([(x[i],y[i]), (x[i]+1,y[i]), (x[i]+1,y[i]+1), (x[i],y[i]+1)], color=(0,0,1), alpha=0.1)
        Gr += text(str(i+1), (x[i]+.5, y[i]+.5), color=(0,0,0), fontsize=11)

        if (x[r[i]], y[r[i]]) != (x[i] + 1, y[i]):
            # r-cut
            Gr += line2d([(x[i]+1,y[i]), (x[i]+1,y[i]+1)], color=(0,0,0))
        else:
            Gr += line2d([(x[i]+1,y[i]), (x[i]+1,y[i]+1)], color=(.3,.3,.3), linestyle='dashed')

        if (x[u[i]], y[u[i]]) != (x[i], y[i] + 1):
            # u-cut
            Gr += line2d([(x[i],y[i]+1),(x[i]+1,y[i]+1)], color=(0,0,0))
        else:
            Gr += line2d([(x[i],y[i]+1),(x[i]+1,y[i]+1)], color=(.3,.3,.3), linestyle='dashed')

        if (x[r_inv[i]], y[r_inv[i]]) != (x[i] - 1, y[i]):
            # r-cut
            Gr += line2d([(x[i],y[i]), (x[i],y[i]+1)], color=(0,0,0))
        if (x[u_inv[i]], y[u_inv[i]]) != (x[i], y[i] - 1):
            # u-cut
            Gr += line2d([(x[i], y[i]), (x[i]+1, y[i])], color=(0,0,0))

    Gr.axes(False)
    Gr.set_aspect_ratio(1)
    return Gr
