function CARRY = pivot( piv_row, CARRY, col )

m = size(col, 1);

%BEGIN_PIVOT = [CARRY, col]

piv_row = piv_row+1;

% Scale pivot row
CARRY(piv_row, :) = CARRY(piv_row, :)/col(piv_row);
col(piv_row) = 1;

% Apply elimination only to CARRY
for i = 1:m
    if i ~= piv_row
        factor = -col(i);
        CARRY(i, :) = CARRY(i, :) + factor*CARRY(piv_row, :);
        col(i) = 0;
       % pause
    end
end

%PIVOT_RESULT = [CARRY, col]

end

