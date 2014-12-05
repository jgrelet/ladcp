function [m,dm,c]=lesqchol(d,g)
% function [m,dm,c]=lesqcholw(d,g)

% fit least squares method to linear problem 
% Use Cholesky transform
% 
%input parameters:
%  d:= data vector ;  g:= model matrix 
% output parameters:
% m=model factors; dm= model data, c=correlation 

n = length(d);
[i,j] = size(g);
if i~=n 
  error(' wrong arguments')
end
[r,b] = chol( g.' * g);
if b~=0 
  m = g(1,:)'+NaN; 
  dm = d+NaN; 
  c = NaN; 
  return, 
end
y = forwardsub(r.' , g.' * d);
m  = backsub(r,y);
if nargout<2 
  return
end
dm = g * m;
if nargout<3 
  return 
end
co = cov([d,dm]);
c  = co(1,2) / sqrt( co(1,1)*co(2,2) );
 
%-------------------------------------------------------------------
function  X = backsub(A,B)

% X = BACKSUB(A,B)  Solves AX=B where A is upper triangular.
% A is an nxn upper-triangular matrix (input)
% B is an nxp matrix (input)
% X is an nxp matrix (output)

[n,p] = size(B);
X = zeros(n,p);

X(n,:) = B(n,:)/A(n,n);

for i = n-1:-1:1,
  X(i,:) = (B(i,:) - A(i,i+1:n)*X(i+1:n,:))/A(i,i);
end
%-------------------------------------------------------------------
function  X = forwardsub(A,B)

% X = FORWARDSUB(A,B))  Solves AX=B where A is lower triangular.
% A is an nxn lower-triangular matrix, input.
% B is an nxp matrix (input)
% X is an nxp matrix (output)

[n,p] = size(B);
X = zeros(n,p);

X(1,:) = B(1,:)/A(1,1);

for i = 2:n,
  X(i,:) = (B(i,:) - A(i,1:i-1)*X(1:i-1,:))/A(i,i);
end
