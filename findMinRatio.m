function index = findMinRatio( b, a )

    %VECTORS_TESTING = [b a]
    %FIND_MIN_RATIO = [b./a]
    %pause

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
    
    %FIND_MIN_RATIO_RESULT = FIND_MIN_RATIO(index)

end

