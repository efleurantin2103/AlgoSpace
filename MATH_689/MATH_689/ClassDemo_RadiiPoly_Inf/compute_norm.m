function norm_value = compute_norm(B, nu, max_iter)
    
    % Initialize the sum
    norm_sum = 0;
    
    % Iterate from n=0 to max_iter
    for n = 0:max_iter
        
      
        % Compute the absolute value of [AF(ā)]_n
        abs_B_a_n = abs(B(n+1));
        
        % Compute the term |[AF(ā)]_n|ν^n
        term = abs_B_a_n * nu^n;
        
        % Add the term to the sum
        norm_sum = norm_sum + term;
    end
    
    % Assign the final sum to the norm value
    norm_value = norm_sum;
end
