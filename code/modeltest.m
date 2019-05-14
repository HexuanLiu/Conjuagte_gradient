%  Generate model problem and run HSCG, CGCG, and GVCG.  Plot the A-norm
%  of the error for each.  Also plot the bound based on the condition number of A
%  and that based on the size of the minimax polynomial on intervals of width 1.e-7
%  or 1.e-14 about the eigenvalues of A.

%  Set up model problem with parameters supplied by user.
n = input('Enter n: ');
rho = input('Enter rho: ');
itmax = input('Enter number of steps to run: ');

%    Set eigenvalues.
lambda = zeros(n,1);
lambda(1) = 0.001; lambda(n) = 1;
for i=2:n-1, lambda(i) = lambda(1) + ((i-1)/(n-1))*(lambda(n)-lambda(1))*rho^(n-i); end;
Lambda = diag(lambda);
%    Choose random orthonormal eigenvectors for A.
[U,R] = qr(randn(n)); A = U*Lambda*U';
for i=1:n-1, for j=i+1:n, A(i,j) = A(j,i); end; end;  % Make sure A is perfectly symmetric.

x0 = zeros(n,1);
flag = 1;
x_true = randn(n,1);
b = A*x_true;

[resid, resest, Tk, Zk, fknorms, inprods, xkdiff, errA, errAest] = hscg(A, b, x0, itmax, flag, x_true);
figure(1)
semilogy([0:itmax], errA/errA(1), '-b', 'LineWidth', 2); hold on
xlabel('Iteration'), ylabel('A-norm of Error'), shg, pause(1)

[resid, resest, Tk, Zk, fknorms, inprods, xkdiff, errA, errAest] = cgcg(A, b, x0, itmax, flag, x_true);
semilogy([0:itmax], errA/errA(1), '--r', 'LineWidth', 2); shg, pause(1)

[resid, resest, Tk, Zk, fknorms, inprods, xkdiff, errA, errAest] = gvcg(A, b, x0, itmax, flag, x_true);
semilogy([0:itmax], errA/errA(1), '-.m', 'LineWidth', 2); shg, pause(1)

kappa = cond(A);
rtkappa = sqrt(kappa);
kappabnd = zeros(itmax,1);
kappabnd(1) = 2;
for k=1:itmax,
  kappabnd(k+1) = 2*( (rtkappa-1)/(rtkappa+1) )^k;
end;
semilogy([0:itmax], kappabnd, '-k')
title('HSCG (solid), CGCG (dashed), GVCG (dash-dot)'), shg, pause(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Now run remez to compute upper bounds based on eigenvalues in small
%  intervals about the eigenvalues of A.

npts = 21;
delta1 = 1.e-7;
[intbnd1,eigAhat1,indices1] = remez(lambda,delta1,npts,itmax);
figure(1)
semilogy([0:itmax],[1;intbnd1],'o'), shg, pause(1)
delta2 = 1.e-14;
[intbnd2,eigAhat2,indices2] = remez(lambda,delta2,npts,itmax);
figure(1)
semilogy([0:itmax],[1;intbnd2],'+'), shg, hold off

%  Now check by running remez_vpa with indices1 and indices2 as the initial
%  guesses for where the minimax polynomials attain their max absolute
%  values.  Could do the whole thing in multiprecision, but that is slow.
%  Hopefully, things won't change much.
%[intbnd1_chk] = remez_vpa(lambda, delta1, npts, itmax, indices1);
%err1 = max(abs((intbnd1 - intbnd1_chk)./intbnd1_chk))
%[intbnd2_chk] = remez_vpa(lambda, delta2, npts, itmax, indices2);
%err2 = max(abs((intbnd2 - intbnd2_chk)./intbnd2_chk))
