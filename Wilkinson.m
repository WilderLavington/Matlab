function y = Wilkinson(a,n)
htr = .5*(a(n-1,n-1)+a(n,n));
dscr = sqrt((.5*(a(n-1,n-1)-a(n,n)))^2+a(n,n-1)^2);
if htr<0
    dscr = -dscr;
end
root1 = htr+dscr;
if root1 == 0;
    root2 = 0;
else
    det = a(n-1,n-1)*a(n,n)-a(n,n-1)^2;
    root2 = det/root1;
end
if abs(a(n,n)-root1) < abs(a(n,n)-root2)
    shift = root1;
else
    shift = root2;
end
y = shift;
end