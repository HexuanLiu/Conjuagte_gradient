function [intbnd,eigAhat,allindices] = remez(eigA, delta, npts, kmax)

%  Find the minimax polynomials of degree 1 to kmax on the union
%  of intervals of width delta about each value in eigA; more precisely,
%  distribute npts points in each of these intervals and find the minimax
%  polynomial on this set of points.  (Note:  npts must be odd.)
%  Return the max norm of each polynomial in array intbnd (of length kmax),
%  the points over which the max is taken in array eigAhat (of length n*npts),
%  and the indices of the points at which each polynomial achieves its max 
%  absolute value in array allindices (of length 2+3+...+kmax+1 = kmax*(kmax+1)/2.

n = length(eigA);
nhat = n*npts;
eigAhat = zeros(nhat,1);
nptsm1o2 = (npts-1)/2;
for j=1:n,
  eigAhat((j-1)*npts+1:j*npts) = eigA(j) + [-nptsm1o2:nptsm1o2]'*(delta/(npts-1));
end;
eigAhat = sort(eigAhat);  % In case intervals overlap, sort the points.
%figure
%semilogx(real(eigA),imag(eigA),'xk'), shg, hold on;
%semilogx(real(eigAhat),imag(eigAhat),'.r'), shg, pause(.1)

%  Loop over polynomial degree k.

errAk = zeros(kmax,1);
intbnd = zeros(kmax,1);
allindices = [];

for k=1:kmax,

%    Initial guess for indices of points at which Bk attains its max abs value.
  if k==1, indices(1) = 1; indices(2) = nhat; end;
  if k > 1,    % Use result from previous degree, adding one more point.
    for kk=k:-1:2,
      indices(kk+1) = indices(kk);
    end;
    indices(2) = fix((indices(3)+1)/2); if indices(2) <= 1, indices(2) = 2; end;
  end;

%    Evaluate Bk at each eigenvalue of Ahat and find where its absolute value
%    is maximal.  If it is at one of the eigenvalues in indices, we are done.
%    Otherwise make an exchange and continue.

  flag = 0;      % Indicates we have not yet found the indices of points where Bk is max.
  while flag==0,

%      Evaluate weights associated with the first form of the barycentric interpolation
%      formula.
    w = ones(k+1,1);
    for i=1:k+1,
      for j=1:k+1,
        if j~=i, w(i) = w(i)*(eigAhat(indices(i)) - eigAhat(indices(j))); end;
      end;
      w(i) = ((-1)^(i-1))/w(i);
    end;

    Bkmax = 0; indxmax = 0; signBkx = zeros(nhat,1); Bkxvals = zeros(nhat,1);
    for m=1:nhat,
      x = eigAhat(m);
      l = 1;
      while l <= k+1 & m~=indices(l),
        l = l+1;
      end;
      if l <= k+1, 
        Bkx = (-1)^(l-1);
      else
        phiofx = 1; 
        for j=1:k+1, phiofx = phiofx*(x-eigAhat(indices(j))); end;
        Bkx = 0;
        for i=1:k+1, Bkx = Bkx + w(i)/(x - eigAhat(indices(i))); end;
        Bkx = phiofx*Bkx;
      end;
%      Bkxchk = 0;
%      for j=1:k+1,
%        termj = (-1)^(j-1);
%        for i=1:k+1,
%          if i~=j,
%            termj = termj*(eigAhat(indices(i))-x)/(eigAhat(indices(i))-eigAhat(indices(j)));
%          end;
%        end;
%        Bkxchk = Bkxchk + termj;
%      end;
%      if abs(Bkx-Bkxchk) > 1.e-8,
%        k, m, indices, x, Bkx, Bkxchk, errBkx = Bkx - Bkxchk, pause
%      end;
      signBkx(m) = sign(Bkx);
      Bkxvals(m) = Bkx;
      Bkx = abs(Bkx);
      if Bkx > Bkmax*(1+1.e-15), Bkmax = Bkx; indxmax = m; end;
    end;

    l = 1;
    while indices(l) < indxmax & l < k+1, l = l+1; end;
    if indices(l)==indxmax,    % We have found the right indices.  Evaluate Bkden.
      phiof0 = 1;
      for j=1:k+1, phiof0 = phiof0*(-eigAhat(indices(j))); end;
      Bkden = 0;
      for i=1:k+1, Bkden = Bkden + w(i)/(-eigAhat(indices(i))); end;
      Bkden = phiof0*Bkden;
%      Bkden = 0;
%      for j=1:k+1,
%        termj = (-1)^(j-1);
%        for i=1:k+1,
%          if i~=j,
%            termj = termj*eigAhat(indices(i))/(eigAhat(indices(i))-eigAhat(indices(j)));
%          end;
%        end;
%        Bkden = Bkden + termj;
%      end;
%      semilogx(eigAhat(indices),k*ones(k+1,1),'o'), shg, pause(.1)
      flag = 1; errk = 1/Bkden; errAk(k) = errk; intbnd(k) = errk;
      allindices = [allindices; indices'];
    else
      if l==1,
        indices(l) = indxmax;
      else
        if signBkx(indices(l-1))==signBkx(indxmax),
          indices(l-1) = indxmax;
        else
          indices(l) = indxmax;
        end;
      end;
    end;
  end;
end;
%hold off
