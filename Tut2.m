A = input("Enter matrix A : ");
C = input("Enter matrix C : ");
b = input("Enter matrix b : ");

m = input("Enter M : ");
n = input("Enter N : ");

array = 1:n;
combinations = nchoosek(array , m);
[ways , chosen] = size(combinations);

ans = 0;

for i = 1 : ways
    basis_matrix = zeros(m , m);
    for j = 1:m  % chosen colm from permutation 
        for k = 1:m % traversing in colm
            basis_matrix(k , j) = A(k ,   combinations(i , j) );
        end
    end

   if(det(basis_matrix) == 0)
       continue;
   end

    basis_solution = basis_matrix\b';
    flag = 0;
    for j = 1:m
        if basis_solution(j) < 0
            flag = 1;
            break
        end
    end
    

    if flag == 1
        continue
    end
       
    temp = 0;

    for j = 1:m
       temp = temp + C( combinations(i , j) )*basis_solution(j);
    end

    if temp > ans
        ans = temp;
    end

end

disp(ans)
