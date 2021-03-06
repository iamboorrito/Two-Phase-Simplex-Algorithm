function [ x_opt, z_opt, CARRY, basis ] = rsimplex(z0, c, A, b, basis)

[m, n] = size(A);

if isempty(basis)
    B = eye(m);
    basis = n+1:n+m;
else
    B = A(:, basis);
end

% Construct CARRY
CARRY = [z0, zeros(1,m); b B];

simplex_exit = 1;

while simplex_exit > 0
    
    % Compute current soln
    x_opt = zeros(n+m, 1);
    x_opt(basis) = CARRY(2:m+1, 1);
    z_opt = -CARRY(1, 1);
    %disp(z_opt)
      
    % Check if any of the objective coefficients will be negative
    pi_t = CARRY(1, 2:m+1);
    
    c_j_new = 0;
    %tableau = [ CARRY(2:m+1, 2:m+1)*A, CARRY(2:m+1, 2:m+1), CARRY(2:m+1, 2:m+1)*b; (c + pi_t*A), pi_t, -z_opt]
    %basis
    
    %pause
    
    % Generate the c_j values until you find a negative one
    for j = 1:n
        
        % Only consider pivoting on columns corresponding to vars not
        % in the basis
        if ismember(j, basis)
            continue;
        end
        
        % This is c_j* in the book, its cB * B^(-1) * A, but
        % cB * B^(-1) is already computed as pi_t
        %j
        c_j_new = c(j) + pi_t*A(:, j);
        
        % Round towards to zero
        if abs( c_j_new ) < 10^(-14)
            c_j_new = 0;
        end
        
        %pause
        
        % Found pivot column corresponding to var not in basis
        if c_j_new < 0 %&& ~ismember(j, basis)
            piv_col = j;
            %pause
            break;
            
        end
        
    end
    
    % Min w = 0 and artificial vars remain, but we don't have a
    % negative obj. fn coefficient, so use pivot op to remove them.
    % Possible TODO: use partial pivoting
    if c_j_new > 0 && abs(z_opt) < 10^(-12) && max(basis) > n
        
        % Get the artificial variable's row from basis
        [~, piv_row] = max(basis);
        
        % Save CARRY in var this time because we need to access its entries
        B_inv = CARRY(2:m+1, 2:m+1);
        
        % Look through the original variables to find one to replace
        % artificial var in basis by generating and searching the row
        % for a nonzero entry we can pivot on
        for i = 1:n
            
            % calculate the row entries to find a nonzero entry
            val = B_inv(piv_row, :)*A(:, i);
            
            %pause
            
            % found nonzero entry, pivot on this one
            if abs(val) > 10^(-8)
                
                % Select the column as entering var
                piv_col = i;
                
                % Update basis before pivoting
                %basis = updateBasis(basis, piv_col, piv_row);
                basis(piv_row) = piv_col;
                
                % Generate the next obj. coeff and column of A* then
                % perform pivot op.
                c_j_new = c(piv_col) + pi_t*A(:, piv_col);
                A_j_new = B_inv*A(:, piv_col);
                CARRY = pivot(piv_row, CARRY, [c_j_new ; A_j_new]);
                
                %pause
                
            end
            
        end
        
        % Skip the rest of the alg. steps since we had to do a special
        % pivot operation
        continue;
        
    end
            
    % No negative coeffs found, soln must be optimal
    if c_j_new >= 0
        x_opt = x_opt(1:n);
        z_opt = -CARRY(1, 1);
        return
    end
    
    % Find pivot row
    A_j_new = CARRY(2:m+1, 2:m+1)*A(:, piv_col);
    piv_row = findMinRatio(CARRY(2:m+1, 1), A_j_new);
    
    % No pivot row found, must be unbounded
    if piv_row == -1
        MSG = sprintf('%s\n', 'Problem unbounded, nonpositive column found')
        x_opt = NaN;
        z_opt = -inf;
        return
    end
    
    % Update basic variables
    basis(piv_row) = piv_col;
    
    % Pivot
    CARRY = pivot(piv_row, CARRY, [c_j_new ; A_j_new]);
    
end

end

