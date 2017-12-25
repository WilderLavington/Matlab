function y = ShiftedInversePowerMethod(A,lambda)
B = A-lambda*eye(size(A));
[L,U] = lu(B);
[r,c] = size(U);
v = lambda*rand(1,r);
z = L\v';
y = U\(z);
end