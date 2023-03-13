% finite diffrence matrixes for 2D grids in horizontal and vertical
% directions (Dh and Dv)
% ALi Gholami

function [Dh,Dv] = drivop(n,m,d)
if ((d < 0) || (d > 2)), error ('Order d must be 0, 1 or 2'), end
% [n,m] = size(f);
N = n*m;

if min(n,m) == 1
    if d==0
        D = eye(max(n,m));
    end
    if d==1
        D = spdiags([[-ones(n-1,1);0] [0;ones(n-1,1)]], [0 1], N, N);
    end
    if d==2
        D = spdiags([[ones(n-2,1);0;0] [0;-2*ones(n-2,1);0] ...
            [0;0;ones(n-2,1)]], [0 1 2], N, N);
    end

else
    if d==0
        Dv = spdiags(ones(N,1),0,N,N);
        Dh = Dv;
    end
    if d==1
        Dv = spdiags([reshape([-ones(n-1,m); zeros(1,m)],N,1) ...
            reshape([zeros(1,m); ones(n-1,m)],N,1)], [0 1], N, N);
        Dh = spdiags([reshape([-ones(n,m-1) zeros(n,1)],N,1) ...
            reshape([zeros(n,1) ones(n,m-1)],N,1)], [0 n], N, N);
    end
    if d==2
        Dv = spdiags([reshape([ones(n-2,m);zeros(2,m)],N,1) ...
            -2*reshape([zeros(1,m) ;ones(n-2,m); zeros(1,m)],N,1) ...
            reshape([zeros(2,m); ones(n-2,m)],N,1)], [0 1 2], N, N);
        Dh = spdiags([reshape([ones(n,m-2) zeros(n,2)],N,1) ...
            -2*reshape([zeros(n,1) ones(n,m-2) zeros(n,1)],N,1) ...
            reshape([zeros(n,2) ones(n,m-2)],N,1)], [0 n 2*n], N, N);
    end
end

if min(n,m) == 1, Dh = D; Dv = []; end





