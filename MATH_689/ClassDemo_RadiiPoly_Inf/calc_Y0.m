function Y0 = calc_Y0(p, AF, nu, N)

a0=p(1);

% Initialize the sum
    norm_sum1 = 0;
    
    % Iterate from n=0 to max_iter
    for i = 0:N
      
        % Compute the absolute value of [AF(ā)]_n
        abs_AF_a_n = abs(AF(i+1));
        
        % Compute the term |[AF(ā)]_n|ν^n
        term = abs_AF_a_n * nu^i;
        
        % Add the term to the sum
        norm_sum1 = norm_sum1 + term;
    end

    n=N+1;

% Second and third summation terms combined
norm_sum23 = 0;
for n = N+1:2*N
        for j = 0:2*N-n
            norm_sum23 = norm_sum23 + nchoosek(N,N-j)*nchoosek(n,n-N+j)*a0^n*nu^n;
        end
end
norm_sum23 = norm_sum23 * (1/(2*abs(a0)));

% Final output
Y0 = norm_sum1 + norm_sum23;

end
