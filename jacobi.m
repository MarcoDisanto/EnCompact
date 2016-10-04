function x1 = jacobi(A,b,x0,tol)

d = diag(A);
Aod = A - diag(d);
r = A*x0 - b;
a = .5;
while (norm(r) > tol)
    plot(x0,'o-')
    pause(1)
    x1 = -a*r./d + x0;
    x0 = x1;
    r = A*x0 - b;
    disp(norm(r))
end