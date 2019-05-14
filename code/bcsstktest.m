%  Read in the bcsstk matrices and run HSCG, CGCG, and GVCG on the 
%  diagonally scaled matrices.  Normally the diagonal would be used
%  as a preconditioner, but to avoid additional rounding errors
%  associated with preconditioned variants, just prescale by the 
%  diagonal.  Plot the results.

for kount=1:6,

%    Load matrix

  if kount==1, load bcsstk14.mat; end;
  if kount==2, load bcsstk15.mat; end;
  if kount==3, load bcsstk16.mat; end;
  if kount==4, load bcsstk17.mat; end;
  if kount==5, load bcsstk18.mat; end;
  if kount==6, load bcsstk27.mat; end;

  A = Problem.A;

%    Prescale by the diagonal.

  M = sparse(diag(diag(A)));
  d = diag(A); rtdinv = 1 ./ sqrt(d);
  MinvA = sparse(diag(rtdinv))*A*sparse(diag(rtdinv));
  [n,n] = size(A);
  MinvAfull = full(MinvA); normMinvA = norm(MinvAfull); condMinvA = cond(MinvAfull);
  A = MinvA/normMinvA; M = sparse(eye(n));  % Avoid questions about preconditioned implementations.

  x0 = zeros(n,1);
  x_true = randn(n,1);   % Set random solution.
  b = A*x_true;
  if kount==1, itmax=800; end;
  if kount==2, itmax=800; end;
  if kount==3, itmax=400; end;
  if kount==4, itmax=3000; end;
  if kount==5, itmax=2000; end;
  if kount==6, itmax=400; end;
  flag = 1;

  subplot(3,2,kount)
  kappabound = ones(itmax+1,1);
  rtkappa = sqrt(condMinvA);
  factor = (rtkappa-1)/(rtkappa+1);
  for k=1:itmax, kappabound(k+1) = 2*factor^k; end;
  semilogy([0:itmax],kappabound,'-k'); 
  axis([0 itmax 1.e-16 1.e1]); 
  ax = gca; ax.YTick = [1.e-16 1.e-12 1.e-8 1.e-4 1.e0]; hold on

%    HSCG

  [resid, resest, Tk, Zk, fknorms, inprods, xkdiff, errA, errAest] = hscg(A, b, x0, itmax, flag, x_true);
  semilogy([0:itmax],errA/errA(1),'-b','LineWidth',2), shg, pause(1), hold on

%    CGCG

  [resid, resest, Tk, Zk, fknorms, inprods, xkdiff, errA, errAest] = cgcg(A, b, x0, itmax, flag, x_true);
  semilogy([0:itmax],errA/errA(1),'--r','LineWidth',2), shg, pause(1)

%    GVCG

  [resid, resest, Tk, Zk, fknorms, inprods, xkdiff, errA, errAest] = gvcg(A, b, x0, itmax, flag, x_true);
  semilogy([0:itmax],errA/errA(1),'-.m','LineWidth',2), shg, pause(1)
  eigTk = eig(Tk(1:itmax,1:itmax));
  mineigTk = min(eigTk), maxeigTk = max(eigTk), condTk = maxeigTk/mineigTk % See if evals of Tk are negative
                                                                           % or if cond(Tk) is greater than that of A.

  xlabel('Iteration'), ylabel('A-norm of Error')
  if kount==1, title(['BCSSTK14: ', 'n = ', int2str(n), ', \kappa = ', num2str(condMinvA,2)]), end;
  if kount==2, title(['BCSSTK15: ', 'n = ', int2str(n), ', \kappa = ', num2str(condMinvA,2)]), end;
  if kount==3, title(['BCSSTK16: ', 'n = ', int2str(n), ', \kappa = ', num2str(condMinvA,2)]), end;
  if kount==4, title(['BCSSTK17: ', 'n = ', int2str(n), ', \kappa = ', num2str(condMinvA,2)]), end;
  if kount==5, title(['BCSSTK18: ', 'n = ', int2str(n), ', \kappa = ', num2str(condMinvA,2)]), end;
  if kount==6, title(['BCSSTK27: ', 'n = ', int2str(n), ', \kappa = ', num2str(condMinvA,2)]), end;
  pause(1)
end;
