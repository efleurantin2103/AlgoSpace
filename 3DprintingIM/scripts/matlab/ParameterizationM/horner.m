%Function Horner’s Method
function Px = horner(x, A, order)

if (order == 0)
    p_i = A(1);
    else
    p_i = A(order+1)*x + A(order);
        for i = 1:(order-1)
            p_i = p_i*x + A(order-i);
        end
end

Px = p_i;