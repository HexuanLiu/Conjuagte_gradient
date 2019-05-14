%  Put eigenvalues in intervals of different widths about the
%  eigenvalues of bcsstk matrices and compare convergence rates
%  using full reorthogonalization.

for kase=1:5,

%    Load bcsstk matrix.

  if kase==1,
    load bcsstk03
    itmax = 820;
    A = Problem.A;
    A = full(A);
    A = A/norm(A);
  end;
  if kase==2, load bcsstk14; itmax = 402; end;
  if kase==3, load bcsstk15; itmax = 554; end;
  if kase==4, load bcsstk16; itmax = 181; end;
  if kase==5, load bcsstk27; itmax = 234; end;
  if kase > 1,
    A = Problem.A;
    M = sparse(diag(diag(A)));
    d = diag(A); rtdinv = 1 ./ sqrt(d);
    MinvA = sparse(diag(rtdinv))*A*sparse(diag(rtdinv));
    A = MinvA;
    A = full(A);
  end;
  [n,n] = size(A);
  npts = 11;
  nhat = n*npts;
  x0 = zeros(nhat,1);
  x_true_n = randn(n,1);
  x_true = zeros(nhat,1);
  for i=1:n, x_true((i-1)*npts+1:i*npts) = x_true_n(i); end;

%    Compute its eigenvalues.

  eigA = eig(A);
    
%    Put eigenvalues in intervals of width delta about the eigenvalues of A.

  for intwidth=1:2

    if intwidth==1,
      delta = 1.e-14;
    else
      delta = 1.e-7;
    end;

    eigAhat = zeros(nhat,1);
    nptsm1o2 = (npts-1)/2;
    for j=1:n,
      eigAhat((j-1)*npts+1:j*npts) = eigA(j) + [-nptsm1o2:nptsm1o2]'*(delta/(npts-1));
    end;
    eigAhat = sort(eigAhat);  % In case intervals overlap, sort the points.
    b = eigAhat.*x_true;

%      Run CG with full reorthogonalization.
    resid = zeros(itmax+1,1); errA = zeros(itmax+1,1);
    Qk = zeros(nhat,itmax+1);
    xk = x0;
    rk = b - eigAhat.*xk;
    rkrk = rk'*rk;
    resid(1) = sqrt(rkrk);
    Qk(:,1) = rk/resid(1);
    err = x_true - xk;
    errA(1) = sqrt(err'*(eigAhat.*err));
    pk = rk;

    for k=1:itmax,
      Apk = eigAhat.*pk;
      aknum = rkrk; akden = Apk'*pk;
      ak = aknum/akden;
      xk = xk + ak*pk;
      rk = rk - ak*Apk;
      for kount=1:2,
      for kk=1:k,
        rk = rk - (Qk(:,kk)'*rk)*Qk(:,kk);
      end;
      end;
      Qk(:,k+1) = rk/sqrt(rk'*rk);
      rk_true = b - eigAhat.*xk;
      resid(k+1) = sqrt(rk_true'*rk_true);
      err = x_true - xk;
      errA(k+1) = sqrt(err'*(eigAhat.*err));
      rkrk = rk'*rk;
      bk = rkrk/aknum;
      pk = rk + bk*pk;
    end;

    subplot(3,2,kase)
    if intwidth==1,
      semilogy([0:itmax],double(errA/errA(1)),'-b','LineWidth',2), 
      xlabel('Iteration'), ylabel('Ahat-norm of Error'), shg, pause(1), hold on
    else
      semilogy([0:itmax],double(errA/errA(1)),'--r','LineWidth',2), 
      if kase==1, title('BCSSTK03'); end;
      if kase==2, title('BCSSTK14'); end;
      if kase==3, title('BCSSTK15'); end;
      if kase==4, title('BCSSTK16'); end;
      if kase==5, title('BCSSTK27'); end;
      axis([0 itmax 1.e-16 1.e2])
      ax = gca; ax.YTick = [1.e-16 1.e-12 1.e-8 1.e-4 1.e0];
      shg, pause(1), hold off
    end;
  end;
end;
