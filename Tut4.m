m  = input("Enter the no. of constraints : ");
n  = input("Enter the no. of Variables  : ");

X = input("Enter the Constraints parameters : ");
B = input("Enter the Constraints : ");
z = input("Enter the objective function parameters : ");


A = zeros(m , n + m);
Z = zeros(1 , n+m);

for i = 1:m
    for j = 1:n
        A(i , j) = X(i, j);
    end
end

for i = 1:m
    A(i , n+i) = 1;
end

for i = 1:n
    Z(1 , i) = z(1 , i);
end

for i = n+1:n+m
    Z(1 , i) = 0;
end


%disp(A)
%disp(Z)
%disp(B)

table (m , n+m+1);

for i = 1:m
    for j = 1:n+m
        table(i , j) = A(i , j);
    end
    table(i , n+m+1) = B(i , 1);
end

temp = 0;


y  = zeros(1 , n+m);
cb = zeros(1 , m);

for j = 1:n+m
    temp = 0;
    for i = 1:m
        temp = temp + table(i , j)*cb(1 , i);
    end
    temp = temp - Z(j);
    y(1 , j) = temp;
end

disp(y);



