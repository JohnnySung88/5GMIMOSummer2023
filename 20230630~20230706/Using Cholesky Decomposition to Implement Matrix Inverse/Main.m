clc
clear

% a=rand(1000)+i*rand(1000);
load('a.mat');
b=a'*a;
c=inv(b);
g = @() inv(b)
timeit(g)
f = @() chol(b)
timeit(f)
e = @() (inv((chol(b))'*chol(b)))
timeit(e)
d = chol(b);
d = inv(d'*d);
sum(sum(abs(c-d)^2, 2))


b=a'*a;
c=inv(b);
e = @() CDI(b)
timeit(e)
d = CDI(b);
sum(sum(abs(c-d)^2, 2))


b=a'*a;
c=inv(b);
e = @() CD(b)
timeit(e)
d = CD(b);
sum(sum(abs(c-d)^2, 2))