function CARRY = pivot( piv_row, CARRY, col )

m = size(col, 1);

% This is to account for CARRY and col having their obj. coefficients in
% the first row.
piv_row = piv_row+1;

% Scale pivot row
CARRY(piv_row, :) = CARRY(piv_row, :)/col(piv_row);
% Set pivot value to 1 in case of round-off error
col(piv_row) = 1;

% Apply elimination
for i = 1:m
    if i ~= piv_row
        factor = -col(i);
        CARRY(i, :) = CARRY(i, :) + factor*CARRY(piv_row, :);
        col(i) = 0;
    end
end

end

