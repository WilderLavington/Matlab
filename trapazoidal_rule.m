function trap = trapazoidal_rule(f,a,b,n)
    trap = 0;
    sub_intervals = linspace(a,b,n+1);
    for i = 1:n
        h = sub_intervals(i+1)-sub_intervals(i);
        trap = trap + h*(f(sub_intervals(i+1))+f(sub_intervals(i)))/2;
    end
end