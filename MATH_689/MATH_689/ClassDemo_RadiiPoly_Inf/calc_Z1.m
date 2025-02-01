function Z1 = calc_Z1(p, nu, N)
    % Extract a0 from the input array a
    a0 = p(1);
    
    % Initialize the sum
    sum_term = 0;
    
    % Loop through n from 1 to N
    for n = 1:N
        % Access the corresponding element from the input array a
        an = p(n+1);
        
        % Add the current term to the sum
        sum_term = sum_term + an * nu^n;
    end
    
    % Multiply the sum by 1/|a0| to get the final result
    Z1 = sum_term * (1/abs(a0));
end