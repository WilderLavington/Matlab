function P = system_solve_mc(l,r,s,d)
    ptm = zeros(length(s));
    ptm(1,1) = s(1);
    ptm(1,2) = l(1);
    %must pad drop-out
    d = [0; d; 0]; 
    for ii = 2:length(s)-1 %rows
        ptm(ii,1) = d(ii-1); 
        ptm(ii,ii) = s(ii);
        ptm(ii,ii-1) = r(ii);
        ptm(ii,ii+1) = l(ii);
    end
    ptm(length(s),length(s)-1) = r(length(s));
    ptm(length(s),length(s)) = s(length(s));
    Q = ptm(2:end-1,2:end-1);
    T = [ptm(2:end-1,1) ptm(2:end-1,end)];
    N = eye(length(Q(:,1)))-Q;
    P = linsolve(N,T);
    %%%%=========================================
    P = P(1,2);
end