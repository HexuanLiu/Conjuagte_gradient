function [T_vpa,Q_vpa] = extendT(A, TJ, QJ, ndigits)

%  Given a Hermitian positive definite matrix A, a Hermitian tridiagonal matrix TJ, 
%  and a set of normalized vectors QJ, this routine uses multiple precision arithmetic 
%  to extend TJ to a larger Hermitian tridiagonal matrix, by constructing new vectors
%  that are orthogonal to each other and to the unconverged Ritz vectors and unconverged
%  cluster vectors corresponding to a perturbed Lanczos recurrence of the form:  
%  A QJ = QJ TJ + betaJ q_{J+1} e_J^T + F_J. 
%
%  Input:
%  A - an n by n Hermitian positive definite matrix.
%  QJ - an n by J+1 array of vectors, each with norm 1 (assumed to come from a finite precision
%       Lanczos computation).
%  TJ - a J+1 by J Hermitian tridiagonal matrix (also assumed to come from a finite precision
%       Lanczos computation).  That is, TJ(1:J,1:J) is a Hermitian tridiagonal matrix and
%       TJ(J+1,J) is the next coefficient in the Lanczos recurrence.
%  ndigits - number of decimal digits to use in computation.  Recommended:  32 or 64.
%
%  Output:  (These are multiprecision arrays.  Use double() to convert to double precision.)
%  T_vpa - a larger Hermitian tridiagonal matrix, whose upper left J by J block is TJ(1:J,1:J).
%          The eigenvalues of T should be close to eigenvalues of A, if parameters
%          eps1 and eps2 (which are printed out) are small.  T will have dimension J+n-m,
%          where m is the number of vectors identified as unconverged Ritz vectors and
%          unconverged cluster vectors.
%      
%  Q_vpa - an n by J+n-m array with orthonormal columns satisfying 
%          A Q_vpa = Q_vpa T_vpa + G, where the norm of G should be small if eps1 and eps2
%          are small.

[n,n] = size(A);
[n,Jp1] = size(QJ);
[Jp1,J] = size(TJ);

digits(ndigits)     % Use ndigits decimal digits.
A_vpa = vpa(A); T_vpa = vpa(TJ);  Q_vpa = vpa(QJ);

%  Make sure that A_vpa and T_vpa are perfectly Hermitian.
for i=1:n-1, for j=i+1:n, A_vpa(i,j) = conj(A_vpa(j,i)); end; end; 
for j=1:J-1, T_vpa(j,j+1) = conj(T_vpa(j+1,j)); end;

%  Make sure that columns of Q_vpa have norm 1.
for j=1:J+1, Q_vpa(:,j) = Q_vpa(:,j)/sqrt(Q_vpa(:,j)'*Q_vpa(:,j)); end;

%  First show how nonorthogonal the columns of QJ are and print out the parameters eps1 and eps2.

ImQpQ_vpa = vpa(eye(J)) - Q_vpa(:,1:J)'*Q_vpa(:,1:J);
nonorth_vpa = sqrt(trace(ImQpQ_vpa'*ImQpQ_vpa));  
nonorth = double(nonorth_vpa);
fprintf('Norm of I - QJtrans*QJ = %12.4d\n',nonorth)

FJ_vpa = A_vpa*Q_vpa(:,1:J) - Q_vpa(:,1:J+1)*T_vpa(1:J+1,1:J);
eps1_vpa = norm(FJ_vpa(:,1));
for j=2:J, 
  normfj_vpa = norm(FJ_vpa(:,j)); 
  if normfj_vpa > eps1_vpa, eps1_vpa = normfj_vpa; end; 
end;
normA = norm(A);
eps1 = double(eps1_vpa)/normA;
fprintf('eps1 = %12.4d\n',eps1)

eps2_vpa = abs(T_vpa(2,1)*Q_vpa(:,1)'*Q_vpa(:,2));
for j=2:J,
  inprod_vpa = abs(T_vpa(j+1,j)*Q_vpa(:,j)'*Q_vpa(:,j+1));
  if inprod_vpa > eps2_vpa, eps2_vpa = inprod_vpa; end;
end;
eps2 = double(eps2_vpa)/normA;
fprintf('eps2 = %12.4d\n',eps2)

%  Plot evals of TJ in bins about evals of A.  Ask user for bin width.

eigA_vpa = sort(eig(A_vpa)); eigT_vpa = sort(eig(T_vpa(1:J,1:J)));
minspace_vpa = min(eigA_vpa(2:n) - eigA_vpa(1:n-1));
minspace = double(minspace_vpa);
fprintf('Min space between evals of A = %12.4d\n',minspace)
width = input('Enter bin width: ');
width_vpa = vpa(width);

%    If bin width is greater than half of min space between eigenvalues,
%    will have to group some evals of A together in a bin.

egroup_vpa(1,1) = eigA_vpa(1); egroup_vpa(1,2) = eigA_vpa(1); ngroups = 1;
for i=2:n,
  if eigA_vpa(i) - width_vpa <= egroup_vpa(ngroups,2)+width_vpa,
    egroup_vpa(ngroups,2) = eigA_vpa(i);
  else
    ngroups = ngroups+1;
    egroup_vpa(ngroups,1) = eigA_vpa(i); egroup_vpa(ngroups,2) = eigA_vpa(i);
  end;
end;

edges_vpa = vpa(zeros(2*ngroups+2,1));
edges_vpa(1) = vpa(-Inf); edges_vpa(2*ngroups+2) = vpa(Inf);
for i=1:ngroups,
  edges_vpa(2*i) = egroup_vpa(i,1) - width_vpa;
  edges_vpa(2*i+1) = egroup_vpa(i,2) + width_vpa;
end;

NT = histc(double(eigT_vpa), double(edges_vpa));
figure(1)
subplot(2,1,1)
bar(NT)
title(['Step ',int2str(J)]), shg, pause(1)

%  Print out lower bound on distance of eigenvalues of extended T to eigenvalues of A
%  based on interlacing of roots of orthogonal polynomials. 

lower_bound_vpa = vpa(0);
i = 1;
for l=1:J,      %  Loop over eigenvalues of TJ  (Could do the same for all submatrices of TJ)

  if i==1 & eigT_vpa(l) < eigA_vpa(1),      % If TJ has an eval < smallest eval of A, so will any extension T.
    dist_vpa = eigA_vpa(1)-eigT_vpa(l);
    if dist_vpa > lower_bound_vpa, lower_bound_vpa = dist_vpa; end;
  elseif i==n & eigT_vpa(l) > eigA_vpa(n),  % If TJ has an eval > largest eval of A, so will any extension T.
    dist_vpa = eigT_vpa(l)-eigA_vpa(n);
    if dist_vpa > lower_bound_vpa, lower_bound_vpa = dist_vpa; end;
  else
    while i < n & eigA_vpa(i) <= eigT_vpa(l), i=i+1; end;  % eigA(i-1) <= eigTJ(l) < eigA(i).
    if l < J & eigT_vpa(l+1) < eigA_vpa(i),                % eigA(i-1) <= eigTJ(l) <= eigTJ(l+1) < eigA(i).
      dist1_vpa = eigT_vpa(l) - eigA_vpa(i-1); dist2_vpa = eigA_vpa(i) - eigT_vpa(l+1);
      dist_vpa = min([dist1_vpa; dist2_vpa]);
      if dist_vpa > lower_bound_vpa, lower_bound_vpa = dist_vpa; end;
      i = i-1;              % Could be more evals of TJ between this pair of evals of A.
    end;
  end;
end;
lower_bound = double(lower_bound_vpa);
fprintf('Lower bound on maxdist = %12.4d\n\n',lower_bound)

%%%%%%  Now start determining how to extend T.  %%%%%%

%  Compute Ritz vectors.

[S_vpa, Theta_vpa] = eig(T_vpa(1:J,1:J));
[eigT_vpa,indx] = sort(diag(Theta_vpa));
S_vpa = S_vpa(:,indx);
for j=1:J, 
  if S_vpa(J,j) < 0, S_vpa(:,j) = -S_vpa(:,j); end;  % Make bottom entries of S positive.
end;
Y_vpa = Q_vpa(:,1:J)*S_vpa;

%-----------------------------------------------------------------------------------------------
%  Try to identify the unconverged Ritz vectors and unconverged cluster vectors.
%  T(J+1,J)*Q(:,J+1) should be almost orthogonal to these and T(J+1,J)*Q(:,J)
%  should be almost equal to a linear combination of these.
%
%  This is the tricky part of the computation.  I don't know if we find the best 
%  set of vectors from the finite precision computation to orthogonalize the new 
%  exact arithmetic Lanczos vectors against.

%    Replace clustered Ritz vectors by cluster vectors.  Save only the unconverged
%    Ritz vectors and cluster vectors.

cluster_width = sqrt(eps)*normA   % Define a cluster based on square root of machine precision.
conv_tol = cluster_width
clusterdef = 1; convdef = 1;
while clusterdef==1 | convdef==1,          % Allow user to change cluster definition or convergence definition
                                           % if not happy with results.
JJ = 0; clusterflag = 0; 
for l=1:J,
  if clusterflag==0,
    if l==J | eigT_vpa(l+1)-eigT_vpa(l) > cluster_width,        % Not part of a cluster
      betaJl_vpa = T_vpa(J+1,J)*S_vpa(J,l);                     % Test to see if converged.
      if betaJl_vpa > conv_tol,
        JJ = JJ+1; YY_vpa(:,JJ) = Y_vpa(:,l);                   % If unconverged, keep it.
      end;
    else
      l1 = l; clusterflag = 1;                          % Start of a cluster
    end;
  else
    if l==J | eigT_vpa(l+1)-eigT_vpa(l) > cluster_width,        % End of a cluster
      wc_vpa = sqrt(S_vpa(J,l1:l)*S_vpa(J,l1:l)');
      if T_vpa(J+1,J)*wc_vpa > conv_tol,
        JJ = JJ+1;
        YY_vpa(:,JJ) = Y_vpa(:,l1:l)*S_vpa(J,l1:l)'/sqrt(S_vpa(J,l1:l)*S_vpa(J,l1:l)');
      end;
      clusterflag = 0;
    end;
  end;
end;
fprintf('No. of unconverged Ritz/cluster vectors = %3i\n',JJ)

%    Check how nearly orthonormal the columns of YY_vpa are.
nonorthYY_vpa = norm(vpa(eye(JJ)) - YY_vpa'*YY_vpa,'fro');
nonorthYY = double(nonorthYY_vpa)
%    Check how nearly orthogonal Q_vpa(:,J+1) is to the columns of YY_vpa.
nonorthQJp1_vpa = norm(Q_vpa(:,J+1)'*YY_vpa);
nonorthQJp1 = double(nonorthQJp1_vpa)
%    Check how close Q_vpa(:,J) is to being a linear combination of columns of YY_vpa.
noncombQJ_vpa = norm( Q_vpa(:,J) - YY_vpa*(YY_vpa'*Q_vpa(:,J)));
noncombQJ = double(noncombQJ_vpa)
%    Ask user if he wants to try a different value for cluster_width or convergence tolerance.
clusterdef = input('Enter new cluster_width?  (0=no,1=yes): ');
if clusterdef==1, cluster_width = input('New cluster_width: '); JJ = 0; YY_vpa = []; end;
convdef = input('Enter new convergence tolerance?  (0=no,1=yes): ');
if convdef==1, conv_tol = input('New convergence tolerance: ');  JJ = 0; YY_vpa = []; end;

end;

%  Form an orthonormal basis for the unconverged Ritz/cluster vectors.

m = JJ;
[W_vpa,Sigma_vpa,X_vpa] = svd(YY_vpa(:,1:m),0);

%-----------------------------------------------------------------------------------------------
  
%  Continue 3-term recurrence with perturbation terms designed to make new vectors
%  orthogonal to each other and to the columns of W.
  
N = n+J-m;
F_vpa = vpa(zeros(n,N));
F_vpa(:,1:J) = A_vpa*Q_vpa(:,1:J) - Q_vpa(:,1:J+1)*T_vpa(1:J+1,1:J);   % Rounding errors from fp computation.

%    Orthogonalize Q(:,J+1) against columns of W and renormalize.

vJp1_vpa = A_vpa*Q_vpa(:,J) - T_vpa(J,J)*Q_vpa(:,J) - T_vpa(J-1,J)*Q_vpa(:,J-1);
if m > 0, vJp1_vpa = vJp1_vpa - W_vpa*(W_vpa'*vJp1_vpa); end;
betaJp1_vpa = sqrt(vJp1_vpa'*vJp1_vpa); 
if m==0, 
  F_vpa(:,J) = vpa(zeros(n,1));
else
  F_vpa(:,J) = W_vpa*W_vpa'*(T_vpa(J+1,J)*Q_vpa(:,J+1)+F_vpa(:,J));
end;
Q_vpa(:,J+1) = vJp1_vpa/betaJp1_vpa;
T_vpa(J,J+1) = betaJp1_vpa; T_vpa(J+1,J) = betaJp1_vpa;

betaj_vpa = T_vpa(J+1,J);
for j=J+2:N+1,
  Q_vpa(:,j) = A*Q_vpa(:,j-1);
  Q_vpa(:,j) = Q_vpa(:,j) - betaj_vpa*Q_vpa(:,j-2);
  alphaj_vpa = Q_vpa(:,j)'*Q_vpa(:,j-1);
  T_vpa(j-1,j-1) = alphaj_vpa;
  Q_vpa(:,j) = Q_vpa(:,j) - alphaj_vpa*Q_vpa(:,j-1);
  if j==J+2,
    for l=1:m,
      F_vpa(:,j-1) = F_vpa(:,j-1) + (Q_vpa(:,j-1)'*A_vpa*W_vpa(:,l) - T_vpa(J+1,J)*Q_vpa(:,j-2)'*W_vpa(:,l))*W_vpa(:,l);
    end;
  else
    for l=1:m,
     F_vpa(:,j-1) = F_vpa(:,j-1) + (Q_vpa(:,j-1)'*A_vpa*W_vpa(:,l))*W_vpa(:,l);
    end;
    F_vpa(:,j-1) = F_vpa(:,j-1) + Q_vpa(:,J+1)*(T_vpa(J+1,J)*Q_vpa(:,j-1)'*Q_vpa(:,J));
  end;
  Q_vpa(:,j) = Q_vpa(:,j) - F_vpa(:,j-1);

%    At this point, the new Lanczos vector should be perfectly orthogonal to W
%    and to the other new Lanczos vectors.  If it is not because we used only
%    64 decimal place arithmetic, orthogonalize one more time to be sure.

  if m > 0,
    newpert_vpa = W_vpa*W_vpa'*Q_vpa(:,j);
  else
    newpert_vpa = vpa(zeros(n,1));
  end;
  Q_vpa(:,j) = Q_vpa(:,j) - newpert_vpa; F_vpa(:,j-1) = F_vpa(:,j-1) + newpert_vpa;
  newpert_vpa = Q_vpa(:,J+1:j-1)*Q_vpa(:,J+1:j-1)'*Q_vpa(:,j);
  Q_vpa(:,j) = Q_vpa(:,j) - newpert_vpa; F_vpa(:,j-1) = F_vpa(:,j-1) + newpert_vpa;

  betaj_vpa = sqrt(Q_vpa(:,j)'*Q_vpa(:,j));
  T_vpa(j,j-1) = betaj_vpa;
  if j < N+1, T_vpa(j-1,j) = betaj_vpa; end;
  if j < N+1, Q_vpa(:,j) = Q_vpa(:,j)/betaj_vpa; end;

%    Plot eigenvalues of T in bins about eigenvalues of A.

  eigT_vpa = sort(eig(T_vpa(1:j-1,1:j-1)));
  subplot(2,1,2)
  NTnew = zeros(length(edges_vpa),1);
  k = 1;
  for i=1:length(eigT_vpa),
    while k < length(edges_vpa) & edges_vpa(k) <= eigT_vpa(i),
      k = k+1;
    end;
    if edges_vpa(k) > eigT_vpa(i) & k > 1, k = k-1; end;
    if edges_vpa(k) <= eigT_vpa(i) & k < length(edges_vpa) & eigT_vpa(i) < edges_vpa(k+1),
      NTnew(k) = NTnew(k)+1;
    end;
  end;
  bar(NTnew),
  title(['Step ',int2str(j-1)])
  pause(1),
end;
title('Eigenvalues of Matrix in Equivalent Exact Arithmetic Computation')

%  Check.

Qorth_vpa = [W_vpa Q_vpa(:,J+1:N)];
nonorth = norm(double(eye(n,n) - Qorth_vpa'*Qorth_vpa),'fro');
fprintf('Norm of I - [W,Q]trans*[W,Q] = %12.4d\n',nonorth)
normF = norm(double(F_vpa),'fro');
fprintf('Norm of perturbations = %12.4d\n',normF)

%  Check distance from eigenvalues of T to those of A.
maxdist_vpa = vpa(0); lmax = 0;
for l=1:length(eigT_vpa),
  distl_vpa = min(abs(eigA_vpa - eigT_vpa(l)));
  if distl_vpa > maxdist_vpa, maxdist_vpa = distl_vpa; lmax = l; end;
end;
fprintf('Maxdist from eval of T to nearest eval of A = %12.4d\n\n',double(maxdist_vpa))
return

%  Run exact CG with matrix Theta and initial vector equal to the first column of S'.

%resid_vpa = vpa(zeros(itmax+1,1));  errT_vpa = vpa(zeros(itmax+1,1));
%[S,Theta] = eig(T(1:N,1:N));
%xk_vpa = vpa(zeros(N,1));
%b_vpa = S(1,:)';
%x_true_vpa = Theta\b_vpa;
%rk_vpa = b_vpa - Theta*xk_vpa;
%resid_vpa(1) = norm(rk_vpa);
%diff_vpa = x_true_vpa - xk_vpa; errA_vpa(1) = sqrt(diff_vpa'*Theta*diff_vpa);
%pk_vpa = rk_vpa;
%bkden_vpa = rk_vpa'*rk_vpa;
%Zk_vpa = vpa(zeros(N,itmax));
%Zk_vpa(:,1) = rk_vpa/sqrt(rk_vpa'*rk_vpa);

%  Iteration
%for k=1:itmax,
%  Apk_vpa = Theta*pk_vpa;
%  akden_vpa = Apk_vpa'*pk_vpa;  ak_vpa = bkden_vpa/akden_vpa;
%  xk_vpa = xk_vpa + ak_vpa*pk_vpa;
%  rk_vpa = rk_vpa - ak_vpa*Apk_vpa;
%  for kount=1:2,
%    for kk=1:k
%      rk_vpa = rk_vpa - (Zk_vpa(:,kk)'*rk_vpa)*Zk_vpa(:,kk);
%    end;
%  end;
%  Zk_vpa(:,k+1) = rk_vpa/sqrt(rk_vpa'*rk_vpa);
%  resid_vpa(k+1) = norm(b_vpa - Theta*xk_vpa);
%  diff_vpa = x_true_vpa - xk_vpa; errA_vpa(k+1) = sqrt(diff_vpa'*Theta*diff_vpa);
%  bknum_vpa = rk_vpa'*rk_vpa; bk_vpa = bknum_vpa/bkden_vpa; bkden_vpa = bknum_vpa;
%  pk_vpa = rk_vpa + bk_vpa*pk_vpa;
%end;

%figure(2)
%subplot(2,1,1)
%semilogy([0:itmax]',resid,'-k', [0:itmax]',resest,'--b', [0:itmax]',double(resid_vpa),'or')
%xlabel('Iteration'), ylabel('2-norm of Residual'), 
%title('FPCG (solid, dashed) and Equivalent Exact CG (o)')
%shg, pause

%figure(3)
%subplot(2,1,2)
%semilogy([0:itmax]',errA/errA(1),'-k', [0:itmax]',errAest/errA(1),'--k', ...
%[0:itmax]',double(errA_vpa/errA_vpa(1)),'or')
%xlabel('Iteration'), ylabel('A- (or T-) Norm of Error'), 
%title('FPCG (solid, dashed) and Equivalent Exact CG (o)')
%shg, pause
