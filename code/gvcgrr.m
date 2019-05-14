%  Pipelined CG

%  User enters HPD matrix A, b, x0, and itmax.
%  Returns residual norms in array resid, norms of updated residuals in array resest.
%  If user sets flag=1 and supplies the true solution x_true, will also return 
%  the A-norm of the error at each step in array errA.
%  Looks at recurrence for zk = (-1)^k rk/||rk|| to see how closely it is satisfied
%  and how nearly orthogonal successive vectors are.  This determines size of
%  intervals in G. paper.  Returns arrays Tk and Zk.
%  
%  This version does residual replacement, as described by Cools, et al.

n = length(b);
resid = zeros(itmax+1,1); resest = zeros(itmax+1,1);
xkdiff = zeros(itmax,1); e1 = zeros(itmax,1); e1(1) = 1;
if flag==1, errA = zeros(itmax+1,1); errAest = zeros(itmax+1,1); end;

%  Initialization
xk = x0;
rk = b - A*xk;
resid(1) = norm(rk); resest(1) = resid(1);
if flag==1, diff = x_true - x0; errA(1) = sqrt(diff'*A*diff); errAest(1) = sqrt(rk'*(A\rk)); end;
uk = rk; wk = A*uk; zeta = norm(b); tau = sqrt(eps); 
theta = sqrt(n)*norm(A,Inf); 
mu = 0;
for i=1:n, 
  mui = 0;
  for j=1:n,
    if A(i,j) ~= 0, mui = mui+1; end;
  end;
  if mui > mu, mu = mui; end;
end;
replace = 0;
zkm1 = zeros(n,1); qkm1 = zeros(n,1); skm1 = zeros(n,1); pkm1 = zeros(n,1);
gkm1 = zeros(n,1); jkm1 = zeros(n,1); fk = 0; pik = 0; sigmak = 0; phik = 0; psik = 0; hk = 0;
replacesteps = [];

%  Iteration
for k=0:itmax-1,
  gammak = uk'*rk; delta = uk'*wk; rhokp1 = norm(rk);
  if k > 0,
    chik = norm(xkm1); pik = norm(pkm1); sigmak = norm(skm1); xik = norm(ukm1);
    omegak = norm(wkm1); phik = norm(qkm1); psik = norm(zkm1); nuk = norm(mkm1);
  end;
  mk = wk;
  vk = A*mk;
  if k > 0,
    betak = gammak/gammakm1; alphak = 1/(delta/gammak - betak/alphakm1);
  else
    betak = 0; alphak = gammak/delta;
  end;
  zk = vk + betak*zkm1;
  qk = mk + betak*qkm1;
  sk = wk + betak*skm1;
  pk = uk + betak*pkm1;
  xkp1 = xk + alphak*pk;
  rkp1 = rk - alphak*sk;
  resest(k+2) = norm(rkp1); resid(k+2) = norm(b-A*xkp1);
  if flag==1, diff = x_true - xkp1; errA(k+2) = sqrt(diff'*A*diff); errAest(k+2) = sqrt(rkp1'*(A\rkp1)); end;
  ukp1 = uk - alphak*qk;
  wkp1 = wk - alphak*zk;
  if k > 0,
    ekm1f = theta*chik + 2*abs(alphakm1)*theta*pik + rhok + 2*abs(alphakm1)*sigmak;
    ekm1h = theta*xik + 2*abs(alphakm1)*theta*phik + omegak + 2*abs(alphakm1)*psik;
    if k > 1,
      ekm1g = theta*xik + 2*abs(betakm1)*theta*pikm1 + omegak + 2*abs(betakm1)*sigmakm1;
      ekm1j = (mu*sqrt(n)+2)*theta*nuk + 2*abs(betakm1)*theta*phikm1 + 2*abs(betakm1)*psikm1;
    end;
    if k==1 | replace==1,
      fk = eps*sqrt((mu*sqrt(n)+1)*theta*chik + zeta) + eps*sqrt(abs(alphakm1)*mu*sqrt(n)*theta*pik) + ...
           eps*sqrt(ekm1f);
      gkm1 = eps*sqrt(mu*sqrt(n)*theta*pik);
      hk = eps*sqrt(mu*sqrt(n)*theta*xik) + eps*sqrt(abs(alphakm1)*mu*sqrt(n)*theta*phik) + eps*sqrt(ekm1h);
      jkm1 = eps*sqrt(mu*sqrt(n)*theta*phik);
      replace = 0;
    else
      fk = fkm1 + abs(alphakm1)*abs(betakm1)*gkm2 + abs(alphakm1)*hkm1 + eps*sqrt(ekm1f) + ...
           eps*abs(alphakm1)*sqrt(ekm1g);
      gkm1 = abs(betakm1)*gkm2 + hkm1 + eps*sqrt(ekm1g);
      hk = hkm1 + abs(alphakm1)*abs(betakm1)*jkm2 + eps*sqrt(ekm1h) + eps*abs(alphakm1)*sqrt(ekm1j);
      jkm1 = abs(betakm1)*jkm2 + eps*sqrt(ekm1j);
    end;
    if fkm1 <= tau*rhok & fk > tau*rhokp1,
      sk = A*pk; qk = sk; zk = A*qk;
      rkp1 = b - A*xkp1; ukp1 = rkp1; wkp1 = A*ukp1;
      replace = 1; replacesteps = [replacesteps; k+1];
    end;
  end;
  zkm1 = zk; qkm1 = qk; skm1 = sk; pkm1 = pk; xkm1 = xk; xk = xkp1; rkm1 = rk; rk = rkp1;
  ukm1 = uk; uk = ukp1; wkm1 = wk; wk = wkp1; gkm2 = gkm1; jkm2 = jkm1; mkm1 = mk; 
  gammakm1 = gammak; alphakm1 = alphak; betakm1 = betak; rhok = rhokp1; fkm1 = fk;
  pikm1 = pik; sigmakm1 = sigmak; phikm1 = phik; psikm1 = psik; hkm1 = hk;
end;
