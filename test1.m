function [] = test1(rho, n)
    randn('state', 1);
    A = randn(n);
    B = randn(n);
    c = randn(n, 1);
    err = @(x) norm(A*x + B*x - c);

    x_dir = (A + B) \ c;
    err(x_dir)

    x = admm2(A, B, c, rho, n);

function [x] = admm(A, B, c, rho, n)

    A = [A; eye(n)];
    B = [B; -eye(n)];
    c = [c; zeros(n, 1)];

    x = zeros(n, 1);
    z = zeros(n, 1);
    u = zeros(2*n, 1);

    L = @(rho, x, z, u) rho/2 * norm(A*x + B*z - c + u)^2;
    err = @(x, z) norm(A*x + B*z - c)

    L(rho, x, z, u)
    for k = 1 : 1e3 
        x = -A \ (B*z - c + u);
        % L(rho, x, z, u)
        z = -B \ (A*x - c + u);
        % L(rho, x, z, u)
        u = A*x + B*z - c + u;
        if mod(k, 13)
            u = 0 * u;
        end
        l(k) = L(rho, x, z, u);
        e(k) = err(x, z);
    end
    semilogy([e; l]', '.-');

function [x] = admm2(A, B, c, rho, n)
    
    x = zeros(n, 1);
    z = zeros(n, 1);
    u = zeros(n, 1);

    L = @(rho, x, z, u) rho/2 * norm(A*x + B*z - c + u)^2;
    err = @(x, z, u) norm(A*x + B*z - c);

    L(rho, x, z, u)

    for k = 1 : 200
        x = x - (A'*A + rho*eye(n)) \ (A'*(A*x + B*z - c) + rho*(x - z + u));
        % L(rho, x, z, u)
        z = z - (B'*B + rho*eye(n)) \ (B'*(A*x + B*z - c) - rho*(x - z + u));
        % L(rho, x, z, u)
        u = x - z + u;
        % L(rho, x, z, u)
        l(k) = L(rho, x, z, u);
        e(k) = err(x, z, u);
        f(k) = norm(x -z);
    end
    semilogy([e; l; f]', '.-');


