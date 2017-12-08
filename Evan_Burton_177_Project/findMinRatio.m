function index = findMinRatio( b, a )

    n = size(a, 1);

    index = -1;
    min_val = inf;

    for i = 1:n
        
        if a(i) > 0 && b(i) >= 0
            
            temp = b(i)/a(i);
            
            if temp < min_val
                min_val = temp;
                index = i;
            end
            
        end
        
    end
    
end

