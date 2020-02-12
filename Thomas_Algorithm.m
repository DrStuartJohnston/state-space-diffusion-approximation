function x = Thomas_Algorithm(a,b,c,F)

%% Solve tridiagonal system using the Thomas Algorithm

n = length(b);                                              %Length of solution vector

c(1) = c(1) / b(1);                                         %Rescale by diagonal               
F(1) = F(1) / b(1);                                         %Rescale by diagonal

%% Scale non-boundary values

for i = 2:n-1
    gamma = b(i) - a(i) * c(i-1);
    c(i) = c(i) /gamma;
    F(i) = (F(i) - a(i) * F(i-1))/gamma;
end

%% Scale boundary values

F(n) = (F(n) - a(n) * F(n-1))/(b(n) - a(n) * c(n-1));

%% Solve system
x(n) = F(n);
for i = n-1:-1:1
    x(i) = F(i)-c(i)*x(i+1);
end