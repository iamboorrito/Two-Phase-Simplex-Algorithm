'''
Created on Oct 10, 2017

@author: Evan Burton
'''
import numpy as np
from builtins import range

class Linear_Program():
    def __init__(self, m, obj, constraints):
        con_rows = constraints.shape[0]
        con_cols = constraints.shape[1]
        
        inequalities = 0
        initialized = False
        
        cols_to_add = np.empty((con_rows, 1))
        self.basic_vars = np.zeros((con_rows))
        
        # count inequalities
        for i in range(con_rows):
            
            val = constraints[i, con_cols-2];
            
            if val != 0:
                
                inequalities += 1
                ei = np.zeros((con_rows, 1))
                                
                if val == 1:
                    
                    ei[i] = -1
                    
                elif val == -1:
                    
                    ei[i] = 1
                    self.basic_vars[i] = con_cols+inequalities-3
                
                if initialized:
                    cols_to_add = np.append(cols_to_add, ei, axis=1)
                else:
                    cols_to_add[:,:] = ei
                    initialized = True
                
        num_art_vars = 0
        aux_row_needed = 0
        
        ####### Add artificial vars to basic var list (negative) #######
        for i in range(con_rows):
            if self.basic_vars[i] == 0:
                self.basic_vars[i] = -(con_cols+inequalities-2+i);
                num_art_vars += 1
                if aux_row_needed == 0:
                    aux_row_needed = 1
                    
        ####### Adjust row size and columns accordingly #######
        self.table = np.zeros((con_rows+1+aux_row_needed, con_cols+inequalities+num_art_vars-1))
        
        self.rows = self.table.shape[0]
        self.cols = self.table.shape[1]
        
        # insert LHS constraint coeffs
        self.table[0:con_rows, 0:con_cols-2] = constraints[:, 0:con_cols-2]

        # insert slack/surplus vars
        self.table[0:con_rows, con_cols-2:con_cols+inequalities-2] = cols_to_add
    
        ####### insert artificial variables #######
        for i in range(con_rows):
            index = self.basic_vars[i]
            if index < 0:
                self.table[i, -int(index)] = 1
                self.table[self.rows-1, -int(index)] = 1
        
        # insert RHS
        self.table[0:con_rows, con_cols+inequalities+num_art_vars-2] = constraints[:, con_cols-1]
    
        # insert objective
        if m.lower() == 'min':
            self.table[con_rows, 0:obj.shape[0]] = obj
        else:
            self.table[con_rows, 0:obj.shape[0]+2] = -obj
        
        ####### insert auxillary objective if not in canonical form #######
    
        for i in range(con_rows):
            if self.basic_vars[i] < 0:
                self.table[self.rows-1, :] -= self.table[i, :]

    def pivot(self, row, col):
        val = self.table[row, col]
        
        # update basic variable
        self.basic_vars[row] = col
        
        self.table[row, :] /= val
        
        # perform elimination
        for i in range(self.table.shape[0]):
            if i != row:
                self.table[i, :] -= self.table[i, col]*self.table[row, :]
    
    def solve(self):
        
        #simplex_alg(self)
                
        last_col = self.table.shape[1]-1
        
        soln = np.zeros((last_col))
        
        for i in range(len(self.basic_vars)):
            soln[int(self.basic_vars[i])] = self.table[i, last_col]
            
        return soln
        
    # TODO: 
    def select_pivot(self):
        return [0, 0]

    def simplex(self):
        
        while min(self.table[self.rows-1, 0:self.cols-1]) < 0:
            
            pivot = self.select_pivot()
            
            if pivot is None:
                print("Problem is unbounded")
                break
            
            self.pivot(pivot[0], pivot[1])
            
            
program = Linear_Program('min', np.array([-6, 2]), np.array([[1, 1, -1, 10], [1, -2, -1, 5]]))

print(program.basic_vars)
print(program.table)
print(program.solve())







