from numpy import array, eye, hstack, ones, vstack, zeros

def cvxopt_solve_minmax(n, a, B, x_min=-42, x_max=42, solver=None):
    c = hstack([zeros(n), [1]])

    # cvxopt constraint format: G * x <= h
    # first,  a + B * x[0:n] <= x[n]
    G1 = zeros((n, n + 1))
    G1[0:n, 0:n] = B
    G1[:, n] = -ones(n)
    h1 = -a

    # then, x_min <= x <= x_max
    x_min = x_min * ones(n)
    x_max = x_max * ones(n)
    G2 = vstack([
        hstack([+eye(n), zeros((n, 1))]),
        hstack([-eye(n), zeros((n, 1))])])
    h2 = hstack([x_max, -x_min])

    c = cvxopt.matrix(c)
    G = cvxopt.matrix(vstack([G1, G2]))
    h = cvxopt.matrix(hstack([h1, h2]))
    sol = cvxopt.solvers.lp(c, G, h, solver=solver)
    return array(sol['x']).reshape((n + 1,))
