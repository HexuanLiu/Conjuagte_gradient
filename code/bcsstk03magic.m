%  Run the CG algorithm on the bcsstk03 problem using finite precision arithmetic.
%  Then extend the tridiagonal matrix to find one with eigenvalues close to 
%  eigenvalues of A, for which exact CG behaves like the finite precision
%  computation.

%  Load input matrix.

load bcsstk03
A = Problem.A;
A = full(A);
A = A/norm(A);  % To avoid scale factors, normalize A.

[n,n] = size(A);
x0 = zeros(n,1);
x_true = randn(n,1);   % Set random solution.
b = A*x_true;          % Compute right-hand side vector.

%  Make sure we have the correct x_true by solving the linear system using multiprecision.
ndigits = 64;
digits(ndigits)
A_vpa = vpa(A); b_vpa = vpa(b);
x_vpa = A_vpa\b_vpa;
x_true = double(x_vpa);

itmax = input('Enter number of steps to run: ');

%  Run itmax steps of CG algorithm for A
%  with solution vector x_true.
flag = 1;

[resid, resest, Tk, Zk, fknorms, inprods, xkdiff, errA, errAest] = hscg(A, b, x0, itmax, flag, x_true);
%[resid, resest, Tk, Zk, fknorms, inprods, xkdiff, errA, errAest] = cgcg(A, b, x0, itmax, flag, x_true);
%[resid, resest, Tk, Zk, fknorms, inprods, xkdiff, errA, errAest] = gvcg(A, b, x0, itmax, flag, x_true);
%[resid, resest, Tk, Zk, fknorms, inprods, xkdiff, errA, errAest] = gvcgwr(A, b, x0, itmax, flag, x_true, 10);

%  Plot the 2-norm of the residuals and the A-norm of the errors.
figure(2)
semilogy([0:itmax],resid/resid(1),'-k', [0:itmax],resest/resest(1),'--k', 'LineWidth',2)
xlabel('Iteration'), ylabel('2-norm of Residual'), hold on, shg, pause(1)
figure(3)
semilogy([0:itmax],errA/errA(1),'-k', [0:itmax],errAest/errAest(1),'--k', 'LineWidth',2)
xlabel('Iteration'), ylabel('A (or T)-norm of Error'), hold on; shg, pause(1)

%  Now call extendT to extend the tridiagonal matrix Tk returned by CG.
ndigits = 64;   % Use 64 decimal digits for extension.
[T_vpa, Q_vpa] = extendT(A, Tk, Zk, ndigits);

%  Now run exact CG for T_vpa, with r0_vpa equal to the first unit vector.
digits(ndigits);
[Np1,N] = size(T_vpa);
r0_vpa = vpa(zeros(N,1)); r0_vpa(1) = vpa(1);
[residx_vpa, errT_vpa] = cg_vpa(T_vpa(1:N,1:N), r0_vpa, itmax, ndigits);
figure(2)
semilogy([0:itmax],double(residx_vpa/residx_vpa(1)),'or'), shg, pause(1), hold off
figure(3)
semilogy([0:itmax],double(errT_vpa/errT_vpa(1)),'or'), shg, pause(1), hold off
