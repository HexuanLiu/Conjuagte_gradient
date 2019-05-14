function [resid_vpa, errT_vpa] = cg_vpa(T_vpa, r0_vpa, itmax, ndigits)
global Zk_vpa

%  Runs itmax steps of cg using ndigits arithmetic, with matrix T_vpa and initial residual
%  r0_vpa.  Uses full reorthogonalization of residual vectors.  
%  Returns residual norms at each step in resid_vpa and T-norm of error at each step
%  in errT_vpa.
%  T_vpa can be the matrix returned by extendT, with r0_vpa equal to the first unit vector.

digits(ndigits)   % Use ndigits decimal digits for computation.
[N,N] = size(T_vpa);

resid_vpa = vpa(zeros(itmax+1,1));  errT_vpa = vpa(zeros(itmax+1,1));
xk_vpa = vpa(zeros(N,1)); rk_vpa = r0_vpa; b_vpa = r0_vpa;
x_true_vpa = T_vpa\b_vpa;
resid_vpa(1) = norm(rk_vpa);
diff_vpa = x_true_vpa - xk_vpa;
errT_vpa(1) = sqrt(diff_vpa'*T_vpa*diff_vpa);
pk_vpa = rk_vpa;
bkden_vpa = rk_vpa'*rk_vpa;
Zk_vpa = vpa(zeros(N,itmax));
Zk_vpa(:,1) = rk_vpa/sqrt(rk_vpa'*rk_vpa);

%  Iteration
for k=1:itmax,
  Tpk_vpa = T_vpa*pk_vpa;
  akden_vpa = Tpk_vpa'*pk_vpa;  ak_vpa = bkden_vpa/akden_vpa;
  xk_vpa = xk_vpa + ak_vpa*pk_vpa;
  rk_vpa = rk_vpa - ak_vpa*Tpk_vpa;

%    Full reorthogonalization
  for kount=1:2,
    for kk=1:k
      rk_vpa = rk_vpa - (Zk_vpa(:,kk)'*rk_vpa)*Zk_vpa(:,kk);
    end;
  end;

  Zk_vpa(:,k+1) = rk_vpa/sqrt(rk_vpa'*rk_vpa);
  resid_vpa(k+1) = norm(b_vpa - T_vpa*xk_vpa);
  diff_vpa = x_true_vpa - xk_vpa; errT_vpa(k+1) = sqrt(diff_vpa'*T_vpa*diff_vpa);
  bknum_vpa = rk_vpa'*rk_vpa; bk_vpa = bknum_vpa/bkden_vpa; bkden_vpa = bknum_vpa;
  pk_vpa = rk_vpa + bk_vpa*pk_vpa;
end;

