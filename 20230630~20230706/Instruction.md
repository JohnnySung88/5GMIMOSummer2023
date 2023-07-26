1. Using C environment
>> mex -setup C

2. Compile the file
>> mex CDI.c

Using Cholesky Decomposition to Implement Matrix Inverse:
>> a=rand(1000)+i*rand(1000);
>> b=a'*a;
# To create an Hermitian matrix [b = b'] -> b = a'a = (a'a)' = a'a = b'; but (a'a) ≠ (aa')
>> c=inv(b);
>> d=CDI(b);
>> sum(sum(abs(c-d)^2, 2))
# The result should be 0 (or very close to 0)

Time Functions: To measure the time required to run a function, use the timeit function.
>> B = @() Function(A); % handle to function
>> timeit(B)
# @()語法表示創建一個匿名函數（匿名函數），該匿名函數會調用Function(A)。
# https://www.mathworks.com/help/matlab/matlab_prog/measure-performance-of-your-program.html