n = 6;
A = rand(n,n);
B = A;
for k = 1:n-1
    for i = k+1:n
        B(i,k) = B(i,k)/B(k,k);
        for j = k+1:n
            B(i,j) = B(i,j) - B(i,k)*B(k,j);
        end
    end
end