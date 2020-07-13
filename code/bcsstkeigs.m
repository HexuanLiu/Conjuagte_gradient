%  Read in the bcsstk matrices and compute and plot their eigenvalues. 

for kount=0:7,

%    Load matrix
  if kount==0, load bcsstk03.mat; end;
  if kount==1, load bcsstk14.mat; end;
  if kount==2, load bcsstk15.mat; end;
  if kount==3, load bcsstk16.mat; end;
  if kount==4, load bcsstk17.mat; end;
  if kount==5, load bcsstk18.mat; end;
  if kount==6, load bcsstk27.mat; end;
  if kount==7,
    n = 48; rho = 0.8;
    lambda1 = .001; lambdan = 1;
    lambda = lambda1*ones(n,1);
    for i=2:n, lambda(i) = lambda(1) + ((i-1)/(n-1))*(lambdan-lambda1)*rho^(n-i); end;
    Lambda = diag(lambda);
    [U,R] = qr(randn(n,n));
    A = U*Lambda*U';
  end; 

  if kount < 7,
  A = Problem.A;
  end;

%    Prescale by the diagonal.
  if kount > 0 & kount < 7,
  M = sparse(diag(diag(A)));
  d = diag(A); rtdinv = 1 ./ sqrt(d);
  MinvA = sparse(diag(rtdinv))*A*sparse(diag(rtdinv));
  else
  MinvA = A;
  end;
  [n,n] = size(A);
  MinvAfull = full(MinvA); normMinvA = norm(MinvAfull); condMinvA = cond(MinvAfull);
  A = MinvA/normMinvA; M = sparse(eye(n));  % Avoid questions about preconditioned implementations.
  eigMinvA = sort(eig(MinvAfull),'ascend')/normMinvA;
  if kount==0, subplot(4,2,1); end;
  if kount > 0, subplot(4,2,kount+1); end;
  if kount==0, semilogy([1:n],eigMinvA,'.', [n-3:n], eigMinvA(n-3:n), 'o'); end;
  if kount > 0, semilogy(eigMinvA,'-', 'LineWidth',2); end;
  ylabel('Evals of A')
  if kount==0, title(['BCSSTK3: ', 'n = ', int2str(n), ', \kappa = ', num2str(condMinvA,2)]), end;
  if kount==1, title(['BCSSTK14: ', 'n = ', int2str(n), ', \kappa = ', num2str(condMinvA,2)]), end;
  if kount==2, title(['BCSSTK15: ', 'n = ', int2str(n), ', \kappa = ', num2str(condMinvA,2)]), end;
  if kount==3, title(['BCSSTK16: ', 'n = ', int2str(n), ', \kappa = ', num2str(condMinvA,2)]), end;
  if kount==4, title(['BCSSTK17: ', 'n = ', int2str(n), ', \kappa = ', num2str(condMinvA,2)]), end;
  if kount==5, title(['BCSSTK18: ', 'n = ', int2str(n), ', \kappa = ', num2str(condMinvA,2)]), end;
  if kount==6, title(['BCSSTK27: ', 'n = ', int2str(n), ', \kappa = ', num2str(condMinvA,2)]), end;
  if kount==7, title(['Model Problem:  n = 48, \kappa = 1000']), end;
  shg, pause(1)
end;
