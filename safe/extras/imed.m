function dist = imed2(A,D)

dist = nan(size(A,2));

D = exp(-D.^2/2);

for i = 1 : size(A,2)
    
    for j = i+1 : size(A,2)
        R = A(:,i) - A(:,j);
        t1 = (R .* A(:,i))'*D;
        t2 = (R .* A(:,j))'*D;

        dist(i,j) = (sum(t1) - sum(t2))/(2*pi);        
    end
    i
end
