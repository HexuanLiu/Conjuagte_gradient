function [resid, resest, Tk, Zk, fknorms, inprods, xkdiff, errA, errAest] = gvcg(A, b, x0, itmax, flag, x_true, replacesteps)

%  Pipelined CG

%  User enters HPD matrix A, right-hand side vector b, initial guess x0, and 
%  max no. of iterations itmax.
%  Returns residual norms in array resid, norms of updated residuals in array resest.
%  If user sets flag=1 and supplies the true solution x_true, will also return 
%  the A-norm of the error at each step in array errA and the quantity sqrt(rk'*(A\rk))
%  in array errAest.
% 
%  Looks at recurrence for zk = (-1)^k rk/||rk|| to see how closely it is satisfied
%  and how nearly orthogonal successive vectors are.  This determines size of
%  intervals in G. paper.  Returns arrays Tk (of size itmax by itmax+1) and Zk (of size n by itmax+1).
%  Also returns array fknorms (of size itmax by 1) containing the norms of each column of
%  A*Zk(:,1:itmax) - Zk*Tk, array inprods (of size itmax-1 by 1) containing the absolute values
%  of successive inner products Tk(k+1,k)*Zk(:,k)'*Zk(:,k+1), and array xkdiff (of size itmax by 1)
%  containing the values norm(xk - (x0 + Zk(:,1:k)*(Tk(1:k,1:k)\(resid(1)*e1(1:k))))).
%  Prints out:  eps1 = norm of largest column of A*Zk(:,1:itmax) - Zk*Tk
%               eps2 = max_k | Tk(k+1,k)*Zk(:,k)'*Zk(:,k+1) |
%               eps3 = max_k norm(xk - (x0 + Zk(:,1:k)*( Tk(1:k,1:k)\(resid(1)*e1(1:k)) )))
%  In exact arithmetic, all of these would be 0.


n = length(b);
resid = zeros(itmax+1,1); resest = zeros(itmax+1,1);
xkdiff = zeros(itmax,1); e1 = zeros(itmax,1); e1(1) = 1;
if flag==1, errA = zeros(itmax+1,1); errAest = zeros(itmax+1,1); end;
if flag~=1, errA = []; errAest = []; end;
Tk = zeros(itmax+1,itmax); Zk = zeros(n,itmax+1);

%  Initialization
xk = x0;
rk = b - A*xk;
resid(1) = norm(rk); resest(1) = resid(1);
if flag==1, diff = x_true - x0; errA(1) = sqrt(diff'*A*diff); errAest(1) = sqrt(rk'*(A\rk)); end;
Zk(:,1) = rk/resest(1);
pk = rk;
sk = A*pk; 
wk = sk; zk = A*wk; 
aknum = rk'*rk;
akden = sk'*pk; ak = aknum/akden; aksave(1) = ak;
Tk(1,1) = 1/ak;

%  Iteration
step = 1;
for k=1:itmax,
  xk = xk + ak*pk;
%      Check to see if xk = x0 + Zk*(Tk\(resid(1)*e1)).
    xkchk = x0 + Zk(:,1:k) * (Tk(1:k,1:k)\(resid(1)*e1(1:k)));
    xkdiff(k) = norm(xk-xkchk);

  rk = rk - ak*sk;
  if k==replacesteps(step),
    if step < length(replacesteps), step = step+1; end;
    wk = A*rk;
  else
    wk = wk - ak*zk;
  end;
  resest(k+1) = norm(rk); resid(k+1) = norm(b-A*xk);
  if flag==1, diff = x_true - xk; errA(k+1) = sqrt(diff'*A*diff); errAest(k+1) = sqrt(rk'*(A\rk)); end;
  Zk(:,k+1) = (-1)^k*rk/resest(k+1);
  Tk(k+1,k) = resest(k+1)/(ak*resest(k));
  if k < itmax, Tk(k,k+1) = Tk(k+1,k); end;
  qk = A*wk;
  bknum = rk'*rk; bk = bknum/aknum; aknum = bknum;
  akm1 = ak;
  ak = aknum/((rk'*wk) - (bk/ak)*aknum);
  if k < itmax, Tk(k+1,k+1) = (1/ak + bk/akm1); end;
  pk = rk + bk*pk;
  sk = wk + bk*sk;
  zk = qk + bk*zk;
end;

fknorms = zeros(itmax,1); inprods = zeros(itmax-1,1);
Gk = A*Zk(:,1:itmax) - Zk*Tk;
gknorm = 0;
for j=1:itmax,
  colnorm = norm(Gk(:,j));
  fknorms(j) = colnorm;
  if colnorm > gknorm, gknorm = colnorm; jmax = j; end;
end;
maxinprod = 0;
for k=1:itmax-1
  inprod = abs(Tk(k+1,k)*Zk(:,k)'*Zk(:,k+1));
  inprods(k) = inprod;
  if inprod > maxinprod, maxinprod = inprod; end;
end;
normA = norm(full(A));
eps1 = gknorm/normA, eps2 = maxinprod/normA
eps3 = max(xkdiff);
if flag==1,
  eps3 = eps3/norm(x_true)
else
  eps3 = eps3/norm(xk)
end;
