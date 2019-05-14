%  Generate model problem and run HSCG, CGCG, and GVCG.  Plot the A-norm
%  of the error for each.  Also plot the bound based on the condition number of A.

%  Set up model problem with parameters supplied by user.
n = input('Enter n: ');
rho = input('Enter rho: ');
itmax = input('Enter number of steps to run: ');

lambda = zeros(n,1);
lambda(1) = 0.001; lambda(n) = 1;
for i=2:n-1, lambda(i) = lambda(1) + ((i-1)/(n-1))*(lambda(n)-lambda(1))*rho^(n-i); end;
Lambda = diag(lambda);
[U,R] = qr(randn(n)); A = U*Lambda*U';
for i=1:n-1, for j=i+1:n, A(i,j) = A(j,i); end; end;  % Make sure A is perfectly symmetric.
x0 = zeros(n,1);
x_true = randn(n,1);
flag = 1;
b = A*x_true;

[resid, resest, Zk, Tk, fknorms, inprods, xkdiff, errA, errAest] = gvcg(A, b, x0, itmax, flag, x_true);
semilogy([0:itmax], errA/errA(1), '-b', 'LineWidth', 2); hold on
%semilogy([0:itmax], errAest/errAest(1), '.c'), shg, pause(1)

gvcgrr
semilogy([0:itmax], errA/errA(1),'--r', 'LineWidth', 2)

[resid, resest, Zk, Tk, fknorms, inprods, xkdiff, errA, errAest] = gvcgwr(A, b, x0, itmax, flag, x_true, replacesteps);
semilogy([0:itmax], errA/errA(1),'-.m', 'LineWidth', 2)
[resid, resest, Zk, Tk, fknorms, inprods, xkdiff, errA, errAest] = gvcgwr(A, b, x0, itmax, flag, x_true, [1:itmax]);
semilogy([0:itmax], errA/errA(1),'-k')
xlabel('Iteration'), ylabel('A-norm of Error')
title('GVCG (solid), with r replace (dashed), with w replace (dash-dot, thin solid)')
