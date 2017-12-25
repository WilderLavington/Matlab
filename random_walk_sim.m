%%%% random whalk simulation to display stationary and non stationary time
%%%% series
clear all
%%%% vanilla random walk 
timestep = 100;
s = rng;
vanilla = rand(timestep,1);
for ii = 2:timestep
    vanilla(ii) = vanilla(ii)*(-1)^randi(100) ;
end
new = rand(timestep,1);
for ii = 2:timestep
    new(ii) = sum(vanilla(1:ii-1));
end
%%difference
check = min(new);
if check < 0
    new = new - check;
end
new = new+1;
newfts = diff(log(new));
figure(1) 
plot(new)
title('Random Walk')
hold on
figure(2) 
plot(log(new))
hold on
title('log scale Random Walk')
hold on 
figure(3)
plot(newfts)
title('Differenced log scale Random Walk')
hold on 
%%%% vanilla random walk with drift
