function [Q,D] = jdiag(A, B, evaOption)
% Joint diagonalization function JDIAG
% JDIAG returns the eigenvectors and the eigenvalues from Aq = dBq
% where q is an eigenvector and d is an eigenvalue, respectively.
% Both are in a range of 1 <= d, q <= dim(A or B).
% Q gives you the joint diagonalization property such that inv(B)*A*Q = Q*D
%                                Q'*A*Q = D
%                                Q'*B*Q = I
% where, I is the identity matrix.
% Although this gives you a similar solution from [Qeve,Qeva] = eig(B\A),
% this can have a relationship described as follows:
%           diag(Q'*A*Q) = diag(Qeve'*A*Qeve)./diag(Qeve'*B*Qeve)
% However, the order of the eigenvalues can be different from each other.
% 
% JDIAG input arguments:
% A                              - a (semi) positive definite matrix
% B                              - a positive definite matrix
% evaOption                      - 'vector' returns D as a vector, diag(D)
%                                - 'matrix' returns D as a diag. matrix
% 
% Latest update   :     6th/December-2018
%
%
%
% For example,
%  A = gsmat(3,1,[3 4 5]);
%  B = gsmat(3,1,[10 20 30]);
%  [Q,D] = JDIAG(A,B);
%
% Q'*A*Q                                Q'*B*Q
% ans =                                 ans = 
%     0.3000   -0.0000    0.0000            1.0000   -0.0000   -0.0000
%    -0.0000    0.2000   -0.0000           -0.0000    1.0000   -0.0000
%     0.0000   -0.0000    0.1667           -0.0000   -0.0000    1.0000
% 
% [Qeve,Qeva] = eig(B\A);
% Qeve'*A*Qeve                          Qeve'*B*Qeve
% ans =                                 ans =
%     3.0000   -0.0000   -0.0000           10.0000    0.0000   -0.0000
%    -0.0000    5.0000   -0.0000            0.0000   30.0000   -0.0000
%    -0.0000    0.0000    4.0000           -0.0000   -0.0000   20.0000
%
% diag(Qeve'*A*Qeve)./diag(Qeve'*B*Qeve)
% ans =
%     0.3000
%     0.1667
%     0.2000
%
if nargin < 3
    evaOption = 'matrix';
end

[Bc,pd] = chol(B,'lower');  % B = Bc*tranpose(Bc)
argname = char(inputname(2));

if pd ~= 0
    error(['Matrix ',argname,' is NOT a positive definite.']);
elseif pd == 0
    % Matrix B is a Positive definite.
    C = Bc\A/transpose(Bc);   % C = inv(Bc)*A*inv(transpose(Bc))
    [U,T] = schur(C);
    X = transpose(Bc)\U;
    
    [dd,dind] = sort(diag(T),'descend');
    D = diag(dd);
    Q = X(:,dind);
end

switch lower(evaOption)
    case 'vector'
        D = diag(D);
    otherwise
end

end