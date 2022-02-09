r"""
The layout problem for square-tiled surfaces
"""

from sage.numerical.mip import MIPSolverException

def soft_layout(r, u):
    r"""
    Return a soft layout as a 4-tuple ``(Graph, x_pos, y_pos, r_cut, u_cut)``
    for the origami ``(r,u)`` or ``None`` if there is no such layout.

    - ``Graph`` - the subgraph used in the layout

    - ``x_pos`` (list of ``n`` integers) - x coordinates of each square

    - ``y_pos`` (list of ``n`` integers) - y coordinates of each square

    - ``r_cut`` (list of ``n`` booleans) - positions of the horizontal cuts

    - ``u_cut`` (list of ``n`` booleans) - positions of the vertical cuts
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
            cut = M.sum(r_cut[i] for i in comp if r[i] not in comp) + \
                  M.sum(r_cut[r_inv[i]] for i in comp if r_inv[i] not in comp) + \
                  M.sum(u_cut[i] for i in comp if u[i] not in comp) + \
                  M.sum(u_cut[u_inv[i]] for i in comp if u_inv[i] not in comp)
            M.add_constraint(cut >= 1)
            continue

        # forbid squares at the same position if any
        has_superposition = False
        for i in range(n):
            for j in range(i):
                if vx[i] == vx[j] and vy[i] == vy[j]:
                    M.add_constraint(x[i] - x[j] + 3 * n * xdiff[i,j] + y[i] - y[j] + 3 * n * ydiff[i,j] >= 1)
                    M.add_constraint(x[j] - x[i] + 3 * n * (1 - xdiff[i,j]) + y[i] - y[j] + 3 * n * ydiff[i,j] >= 1)
                    M.add_constraint(x[i] - x[j] + 3 * n * xdiff[i,j] + y[j] - y[i] + 3 * n * (1 - ydiff[i,j]) >= 1)
                    M.add_constraint(x[j] - x[i] + 3 * n * (1 - xdiff[i,j]) + y[j] - y[i] + 3 * n * (1 - ydiff[i,j]) >= 1)
                    has_superposition = True
                if has_superposition:
                    break
            if has_superposition:
                break

        if not has_superposition:
            vx = M.get_values(x)
            vy = M.get_values(y)
            return G, [int(vx[i]) for i in range(n)], [int(vy[i]) for i in range(n)], [bool(vr_cut[i]) for i in range(n)], [bool(vu_cut[i]) for i in range(n)]

def soft_plot(r, u):
    r"""
    Return a plot of the origami ``(r, u)`` using a soft layout (if any).
    """
    n = len(r)
    r_inv = [None] * n
    u_inv = [None] * n
    for i in range(n):
        r_inv[r[i]] = i
        u_inv[u[i]] = i

    ans = soft_layout(o)
    if ans is None:
        raise ValueError('no soft layout')
    G, x, y, r_cut, u_cut = ans
    Gr = Graphics()
    for i,(xi,yi) in enumerate(zip(x,y)):
        Gr += polygon2d([(xi,yi),(xi+1,yi),(xi+1,yi+1),(xi,yi+1)], color=(0,0,1), alpha=0.1)
        Gr += text(str(i+1), (xi+.5, yi+.5), color=(0,0,0), fontsize=11)
        if r_cut[i]:
            Gr += line2d([(xi+1,yi), (xi+1,yi+1)], color=(0,0,0))
        else:
            Gr += line2d([(xi+1,yi), (xi+1,yi+1)], color=(.3,.3,.3), linestyle='dashed')
        if u_cut[i]:
            Gr += line2d([(xi,yi+1),(xi+1,yi+1)], color=(0,0,0))
        else:
            Gr += line2d([(xi,yi+1),(xi+1,yi+1)], color=(.3,.3,.3), linestyle='dashed')
        if r_cut[r_inv[i]]:
            Gr += line2d([(xi,yi), (xi,yi+1)], color=(0,0,0))
        if u_cut[u_inv[i]]:
            Gr += line2d([(xi,yi), (xi+1,yi)], color=(0,0,0))

    Gr.axes(False)
    Gr.set_aspect_ratio(1)
    return Gr
