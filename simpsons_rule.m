function trap = simpsons_rule(f,a,b,n)
    trap = 0;
    sub_intervals = linspace(a,b,n+1);
    for i = 1:2:n
        h = sub_intervals(i+1)-sub_intervals(i);
        x_0 = sub_intervals(i);
        x_1 = (sub_intervals(i) + sub_intervals(i+1))/2;
        x_2 = sub_intervals(i+1);
        trap = trap + h*(f(x_0) + 4*f(x_1) + f(x_2))/3;
    end
end