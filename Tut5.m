function x = simplex (A_leq, A_eq, A_geq, b_leq, b_eq, b_geq, cT)
    % count = size(A, 1) for {leq, eq, geq}
    count_leq = size(A_leq, 1);
    count_eq = size(A_eq, 1);
    count_geq = size(A_geq, 1);

    A_leq = [A_leq zeros(count_leq, count_geq) eye(count_leq) zeros(count_leq, count_eq+count_geq)];
    A_eq = [A_eq zeros(count_eq, count_leq+count_geq) eye(count_eq) zeros(count_eq, count_geq)];
    A_geq = [A_geq -eye(count_geq) zeros(count_geq, count_leq+count_eq) eye(count_geq)];

    sz = size(cT, 2);

    A = [A_leq; A_eq; A_geq];
    cT = [cT zeros(1, count_leq+count_eq+2*count_geq)];
    b = [b_leq' b_eq' b_geq']';
    [m, n] = size(A);
    
    % Big M method for processing artificial variables
    M = 0;
    for i = 1: n
        M = max(M, 1000*abs(cT(i)));
    end
    for i = n-count_geq-count_eq+1:n
        cT(i) = -M;
    end

    v = n-m+1:n;
    xb = b;
    y = A;

    diff = zeros(n, 1);
    x = zeros(sz, 1);
    
    while true
        smalldiff = 0;
        entering_idx = 0;
        for i = 1:n
            diff(i) = cT(v) * y(:, i) - cT(i); 
            
            if diff(i) < smalldiff
                smalldiff = diff(i);
                entering_idx = i;
            end
        end
        if smalldiff == 0
            break;
        end
        
        ratios = zeros(m, 1);
        smallval = Inf;
        leaving_idx = 0;
        for j = 1:m
            ratios(j) = xb(j) / y(j, entering_idx);
            if ratios(j) < smallval && 0 < ratios(j)
                smallval = ratios(j);
                leaving_idx = j;
            end
        end
        if leaving_idx == 0
            disp("The system has an unbounded solution, returning x=0")
            return
        end
      
        % y
        % diff'
        % ratios
        % leaving_idx
        % entering_idx
        
        pivot_value = y(leaving_idx, entering_idx);
        v(leaving_idx) = entering_idx;
        xb(leaving_idx) = xb(leaving_idx)/pivot_value;
        y(leaving_idx, :) = y(leaving_idx, :)/pivot_value;

        for j = 1:m
            if j == leaving_idx
                continue
            end
            factor = y(j, entering_idx) / y(leaving_idx, entering_idx);
            y(j, :) = y(j, :) - y(leaving_idx, :) * factor;
            xb(j) = xb(j) - xb(leaving_idx) * factor;
        end 
    end
    for i = 1:m
        if v(i) <= sz
            x(v(i)) = xb(i);
        end
    end
end