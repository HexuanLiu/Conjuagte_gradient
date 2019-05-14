%  Load bcsstk03 and run HSCG, CGCG, and GVCG.  Plot the A-norm
%  of the error for each.  Also plot the bound based on
%  the condition number of A.

%  Load input matrix.

load bcsstk03
A = Problem.A;
A = full(A);
A = A/norm(A);  % To avoid scale factors, normalize A.

[n,n] = size(A);
x0 = zeros(n,1);
x_true = randn(n,1);   % Set random solution.
b = A*x_true;          % Compute right-hand side vector.
flag = 1;
itmax = 1200;

%  HSCG

[resid, resest, Tk, Zk, fknorms, inprods, xkdiff, errA, errAest] = hscg(A, b, x0, itmax, flag, x_true);

figure(1)
semilogy([0:itmax], errA/errA(1), '-b', 'LineWidth', 2); hold on
xlabel('Iteration'), ylabel('A-norm of Error'), shg, pause(1)
%figure(2)
%semilogy([1:itmax], fknorms, '.b'); hold on
%xlabel('Iteration'), ylabel('fknorms'), shg, pause(1)
%figure(3)
%semilogy([1:itmax-1], inprods, '.b'); hold on
%xlabel('Iteration'), ylabel('inprods'), shg, pause(1)
%figure(4)
%semilogy([1:itmax], xkdiff, '.b'); hold on
%xlabel('Iteration'), ylabel('xkdiff'), shg, pause(1)

%  CGCG

[resid, resest, Tk, Zk, fknorms, inprods, xkdiff, errA, errAest] = cgcg(A, b, x0, itmax, flag, x_true);
figure(1)
semilogy([0:itmax], errA/errA(1), '--r', 'LineWidth', 2); hold on
shg, pause(1)
%figure(2)
%semilogy([1:itmax], fknorms, '.r'); shg, pause(1)
%figure(3)
%semilogy([1:itmax-1], inprods, '.r'); shg, pause(1)
%figure(4)
%semilogy([1:itmax], xkdiff, '.r'); shg, pause(1)

%  GVCG

[resid, resest, Tk, Zk, fknorms, inprods, xkdiff, errA, errAest] = gvcg(A, b, x0, itmax, flag, x_true);
figure(1)
semilogy([0:itmax], errA/errA(1), '-.m', 'LineWidth', 2); hold on
shg, pause(1)
%figure(2)
%semilogy([1:itmax], fknorms, '.m'); shg, pause(1), hold off
%figure(3)
%semilogy([1:itmax-1], inprods, '.m'); shg, pause(1), hold off
%figure(4)
%semilogy([1:itmax], xkdiff, '.m'); shg, pause(1), hold off

kappa = cond(A);
rtkappa = sqrt(kappa);
kappabnd = zeros(itmax,1);
kappabnd(1) = 2;
for k=1:itmax,
  kappabnd(k+1) = 2*( (rtkappa-1)/(rtkappa+1) )^k;
end;
figure(1)
semilogy([0:itmax], kappabnd, '-k')
title('HSCG (solid), CGCG (dashed), GVCG (dash-dot)'), shg, hold off
