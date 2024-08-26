# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 21:18:32 2024

@author: rravu
"""

class Matrix():
    
    def __init__(self, matrix):
        self.matrix = matrix;
        self.NumRows = len(self.matrix)
        self.NumCols = len(self.matrix[0])
        
    class RowOperations: # contains various static methods for matrix row operations
        '''
        When solving a linear system of equations modeled by a matrix, various row operations can be used.
        For computing the row echelon form or reduced row echelon form (REF/RREF) of a matrix, the following
        3 row operations can be performed:
            
            - Row Swapping: Swapping the positions of two given rows in the matrix
            
            - Scalar-Row Multipication: Multiplying all the entries in a given row by a nonzero scalar
            
            - Current Row -/+ Constant*Pivot Row: Adding the entries of the pivot row to the corresponding entires
                                                  of the current row. A nonzero scalar can be multiplied to the
                                                  entries of the pivot row as well.
        '''

        def AddMultipleOfRow(CurrRow, PivotRow, scalar):
            '''
            This method adds the entries of the given pivot row to the corresponding entries of the current row.
            A nonzero scalar can be used to augment the entries of the pivot row during the addition of the two rows
            for optimal row reduction.
            
            Three arguments are passed: CurrRow, which is the current row of the matrix we wish to operate on, 
            pivot row, which is the row containing the current pivot of the matrix (a pivot is the first nonzero
            entry in a row), and a scalar, containing a nonzero scalar.
            ''' 
            
            for i in range(len(CurrRow)):
                CurrRow[i] -= PivotRow[i]*scalar
        
        def ScalarMultiplyRow(row, scalar):
            '''
            This method multiplies all the entries of a given row by a specified nonzero constant
             
            Two arguments are passed: row, containg the row we wish to operate on, and scalar, which is a nonzero
            constant
            '''
            for i in range(len(row)):
                row[i] *= scalar
        
        def RowSwap(Matrix, row1Index, row2Index): 
            '''
            This method swaps the positions of two rows within a given matrix. For this calculator,
            matrices are represented as 2-dimensional lists, with matrix rows being represented as 1-dimensional lists.
            Hence, the swapping of two rows entails reassigning the indices of the rows in the matrix.
            
            Three arguments are taken: the matrix, the index of the first row, and the index of the second row
            '''
            Matrix[row1Index], Matrix[row2Index] = Matrix[row2Index], Matrix[row1Index]
        
        def FirstNonZeroPivotAfter(matrix, rowindex, pivotIndex):
            '''
            Although not strictly a row operation, this method acts as a helper function for locating the
            first row that does NOT contain a 0 in the current pivot column. The pivot column is the column location
            of the pivot entry. The pivot of a row in a matrix is the first nonzero entry in that row. We use pivots
            to keep track of which rows we must reduce in order to get an efficient and clean solution to a linear
            system.
             
            Three arguments are passed: the matrix, the index of the row, and the current pivot index, which holds
            either the row or column index of the current pivot.
            '''
            
            for i in range(rowindex+1, len(matrix)):
                if matrix[i][pivotIndex] != 0:
                    return i
            return -1 # if no row with a pivot column entry != 0 was found, return 0
        
        def LastNonZeroRowIndex(matrix, zero_row): 
            '''
            This method is a helper method that locates the last occurrence of a nonzero row. A zero row is a 
            row of the matrix that contains all zeros. This method traverses the matrix backwards to locate the
            last nonzero row. This method is useful for row swapping during the process of converting a matrix to
            echelon form.
            
            Two arguments are passed: the matrix, as well as a default zero row, which is a list the maximum possible
            number of zeros for a given row in the matrix.
            '''
            for j in range(len(matrix)-1, -1, -1): # reverse iteration
                if matrix[j] != zero_row:
                    return j
            return None # if no nonzero row was found, return None
        
    def EchelonForm(self, choice): # Computes either the row echelon form or reduced row echelon form of a matrix
        '''
        Converting a matrix to "Echelon Form" gives us information relating to the solution set of the linear system
        modeled by the matrix. You can convert a matrix either to standard row echelon form (REF) or reduced row
        echelon form (RREF). 
        
        REF (Row Echelon Form) Requirements:
            
            - All zero rows (rows that contain all zeros) must be at the BOTTOM of the matrix
            
            - The pivot entries of the matrix (the first nonzero entry in each row) must be to the right of the pivot entry of the row above it

        RREF (Reduced Row Echelon Form) Requirements:
            
            - The matrix must be in REF (Standard Reduced Row Echelon Form)
            
            - The leading entries/pivots of each row must be 1's 
            
            - The leading 1's must contain 0's everywhere else in their column.
            
        
        It should be noted that although standard row echelon form (REF) does not require the leading entries
        of each row to be 1's, unlike RREF, we still convert the leading entries of each row to be 1 just for
        simplicity.
            
        '''
          
        NumPivots = self.NumCols if self.NumCols < self.NumRows else self.NumRows
        
        zero_row = [0 for num in range(self.NumCols)]
        
        RowOps = Matrix.RowOperations
        
        for i in range(NumPivots):
            
            if self.matrix[i]==zero_row:
                LastNonZero = RowOps.LastNonZeroRowIndex(self.matrix, zero_row)
                if LastNonZero != None:
                    RowOps.RowSwap(self.matrix, i, LastNonZero)
                
                else:
                    return "Your matrix is a null matrix"
                
            if self.matrix[i][i]==0:
                FirstNonZero = RowOps.FirstNonZeroPivotAfter(self.matrix, i, i)
                if FirstNonZero != -1:
                    RowOps.RowSwap(self.matrix, i, FirstNonZero)
                
            if self.matrix[i][i] != 1:
                RowOps.ScalarMultiplyRow(self.matrix[i], 1/self.matrix[i][i])
                
            if choice == 'j':  # If standard row-echelon form is specified
                for j in range(i+1, len(self.matrix)):
                    if self.matrix[j][i] != 0:
                        RowOps.AddMultipleOfRow(self.matrix[j], self.matrix[i], self.matrix[j][i]/self.matrix[i][i])
                        
            elif choice == 'k': # if reduced row-echelon form is specified
                 for j in range(0, len(self.matrix)):
                     if j==i:
                         continue
                     if self.matrix[j][i] != 0:
                         RowOps.AddMultipleOfRow(self.matrix[j], self.matrix[i], self.matrix[j][i]/self.matrix[i][i])
                
        return self
             
    def __mul__(self, OtherMatrix): # Matrix Multiplication 
        '''
        This method multiples two given matrices. A result matrix is allocated
        for storing the elements resulting from multiplying the two matrices.
        The result matrix is then returned.
        
        
        Matrix Multiplication Restrictions:
            
            - Matrix multiplication can only be done when the number of COLUMNS of the first matrix
              is equal to the number of ROWS of the second matrix
            
            - As a result, matrix multiplication is NOT commutative 
            (i.e. You can multiply a 4x2 matrix and 2x3 matrix in that specific order, but NOT the other way around)
        '''
        
        result = [[0 for elem in OtherMatrix.matrix[0]] for row in self.matrix]
        
        for i in range(len(self.matrix)):
            for j in range(len(OtherMatrix.matrix)):
                for v in range(len(OtherMatrix.matrix[0])):
                    result[i][v] += (self.matrix[i][j]*OtherMatrix.matrix[j][v])
                    
        result = Matrix(result)
        
        
        return result
                   
    def Identity(self):
        '''
        This method acts as a helper function by returning an identity matrix of same dimension as the
        given matrix. The identity matrix is a matrix containing all zeros except for the main diagonal, which contains
        all 1's. The identity matrix is commonly used in matrix algebra for operations such as the power and inverse of
        a matrix.
        '''
        return Matrix([[1 if i==j else 0 for i in range(len(self.matrix[0]))] for j in range(len(self.matrix))])
    
    def __pow__(self, Exponent):
        '''
        The power of a matrix entails repeatedly multiplying a matrix by itself k times, where k is some scalar.
        When the specified exponent is positive, traditional matrix multiplication is performed repeatedly, k times.
        
        When a matrix is raised to the zero power, an identity matrix of equivalent dimension is returned.
        
        When a matrix is raised to the negative power, the inverse of a matrix must be computed, as a negative exponent
        is defined as the multiplicative inverse of the base matrix.
        
        RESTRICTIONS: only SQUARE matrices can be taken to a given power.
        '''
        
        
        if Exponent==0:
            return self.Identity()
        
        if Exponent > 0:
            
            result = self.Identity()
            for num in range(Exponent):
                result *= self
            
            return result
        
        if Exponent < 0:
            
            if self.Determinant()==0:
                return f"Your matrix has no inverse because the determinant is 0. Hence, we cannot compute your matrix to power {Exponent}"
                
            else:
                
                result = (self.Inverse())**(-1*Exponent)
                
                return result
    
    def __add__(self, OtherMatrix): 
        '''
        This function adds two matrices together.
        An additional matrix (2D List) is allocated for storing the
        resulting elements from adding the two input matrices.
        
        The resulting matrix is then returned.
        '''
        result_matrix = [[] for row in self.matrix]
        
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[0])):
                result_matrix[i].append(self.matrix[i][j] + OtherMatrix.matrix[i][j])
            
        result_matrix = Matrix(result_matrix)
        
        return result_matrix
    
    def __sub__(self, OtherMatrix):
        '''
        This function subtracts two given matrices.
        An additional matrix (2D List) is allocated for storing the
        resulting elements from subtracting the two input matrices.
        
        The resulting matrix is then returned.
        '''
        result_matrix = [[] for row in self.matrix]
        
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[0])):
                result_matrix[i].append(self.matrix[i][j] - OtherMatrix.matrix[i][j])
            
        result_matrix = Matrix(result_matrix)
        
        return result_matrix
    
    def Scalar_Multiplication(self, Scalar): # Scalar Multiplication 
        '''
        This method multiplies all the entries of the matrix by a specified nonzero scalar quantity.
        '''
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[0])):
                self.matrix[i][j] *= Scalar
        
        return self
    
    def Transpose(self):
        '''
        The transpose of a matrix entails swapping the rows and columns of the matrix. This method
        swaps the rows and columns of the matrix by allocating space for an additional matrix, used for storing
        the entries of the transpose of the matrix. The resulting matrix is then returned.
        '''
        result = [[] for column in self.matrix[0]]
        
        for i in range(len(self.matrix[0])): 
            for j in range(len(self.matrix)): 
                result[i].append(self.matrix[j][i]) 
                
        result = Matrix(result)
        return result
    
    def Rank(self):
        '''
        The rank of a matrix provides information regarding the dimension of the column space of the matrix. The column
        space of a given matrix is the set of ALL possible linear combinations of the column vectors of the matrix
        (linear combinations are, like the name implies, just combinations of vectors formed by adding or substracting
        the vectors from each other, with possible nonzero scalars).
        
        The rank of a matrix can be calculated by examining the number of nonzero rows of the matrix in either its
        row echelon form or reduced row echelon form.
        '''
        
        self.EchelonForm("j") # specifying choice "j" reduces the matrix to standard row echelon form (REF)
                
        zero_row = [0 for num in self.matrix[0]] # zero rows are rows containing all zeros
        
        return len([row for row in self.matrix if row != zero_row]) # counting the number of non-zero rows
    
    def Determinant(self):
        '''
        The determminant of a matrix is a special single numeric value that can provide information relating to the
        linear map or area of a region resulting after a linear transformation. There are numerous methods to compute
        the determinant, but this method will use laplace/cofactor expansion.
        
        The determinant of a 2x2 matrix can be found by taking the difference of the product of the its diagonal entries
        
        example:
            
            A = [[1, 2],
                 [3, 4]]
            
            det(A) = 1(4)-3(2) = -2
            
        Laplace expansion computes higher dimensional matrices by constantly creating minor matrices of entries along a certain
        row/column, and computing the weighted sum of the determinant of all minor matrices.
        
        For this method, we implement a recursive laplace expansion solution, using the 2x2 matrix as a base case.
        
        Restrictions: The determinant can only be computed for SQUARE matrices
        '''
        
        # This is a recursive solution implementing laplace expansion
    
        def Minor_Matrix(matrix, index):  # creates minor matrices for each entries along the first row of the matrix (and sub matrices)
            return Matrix([[row[i] for i in range(len(row)) if i != index] for row in matrix.matrix[1:]])
            
        if len(self.matrix)==2: # The base case where the given matrix is 2x2 in dimension. Simply return the difference between
                                # product of diagonal entries
            return (self.matrix[0][0]*self.matrix[1][1])-(self.matrix[0][1]*self.matrix[1][0])
        
        else:
            # If input matrix is larger than 2x2, this method will recursively call itself, creating more minor matrices for 
            # sub matrices of the matrix. The weighted total is then taken to calculate the overall determminant of the matrix
            CurrMatrix = self
            Total_Det = 0
            for i in range(len(CurrMatrix.matrix[0])):
                Minor = Minor_Matrix(CurrMatrix, i)
                sign = (-1)**i
                Det_Minor = CurrMatrix.matrix[0][i]*sign*Matrix.Determinant(Minor) # recursively calling Matrix.Determinant
                Total_Det += Det_Minor # totalling the determinant across all minor matrices
                
            return Total_Det # returning the determinant value
    
    def Inverse(self):
        '''
        A matrix X is deemed "Invertible" if there exists another matrix, Y, such that 
        X*Y = Y*X = I (where I is the identity matrix). If such a matrix exists, then matrix Y is deemed 
        the "Inverse" of matrix X.
        
        Two prominently used methods to compute the inverse of a matrix (if it exists) are the adjugate formula, and 
        gauss-jordan elimination. In order to showcase both methods, while not completely sacrificing computational complexity,
        the adjugate formula is used in this method for computing matrices that are between 2x2 and 4x4 in dimension (inclusive). 
        For matrices larger than 4x4, gauss-jordan elimination is used for computing the inverse of the matrix.
        
        Restrictions:
            
            - The inverse can only be computed for SQUARE matrices
            
            - A matrix is not invertible if its determinant is equivalent to 0
        '''
        
        if self.Determinant()==0: 
            return "Your matrix does not have an inverse due to the determinant being 0"
        
        if len(self.matrix)==2:
            
            DeterminantInverse = 1/self.Determinant()
            self.matrix[0][0], self.matrix[1][1] = self.matrix[1][1], self.matrix[0][0]
            self.matrix[0][1] *= -1
            self.matrix[1][0] *= -1
            return self.Scalar_Multiplication(DeterminantInverse)
        
        if len(self.matrix) > 2 and len(self.matrix) <= 4: 
            '''
            if the matrix is between 2x2 and 4x4 in dimension (inclusive), 
            then the adjugate formula is used for computing the inverse
            '''                                              
            
            def Cofactor(matrix, column_num, row_num):
               Cofactor = [[matrix[j][i] for i in range(len(matrix[j])) if i != column_num] for j in range(len(matrix)) if j != row_num]
               return Matrix(Cofactor)
           
            
            Cofac_Matrix = Matrix([row[:] for row in self.matrix]) 
            # setting Cofac_Matrix directly to self causes errors, as the following operations then update both
            # matrices in memory simultaneously. Hence, we must set the cofactor matrix to the
            # inherent 2D list of the original matrix, and create a new instance of class Matrix
           
            for i in range(len(self.matrix)):
                for j in range(len(self.matrix[0])):
                    Cofac_Matrix.matrix[i][j] = ((-1)**((i+j)))*Matrix.Determinant(Cofactor(self.matrix, j, i))
                    
            InverseDeterminant = 1/self.Determinant()
            
            Cofac_Matrix = Matrix.Transpose(Cofac_Matrix)
            
            inverse = Matrix.Scalar_Multiplication(Cofac_Matrix, InverseDeterminant)
    
            return inverse
        
        else:
            '''
            The adjoint/adjugate method for computing the inverse of a matrix becomes computationally
            inefficient for matrices of larger size. Hence, for matrices larger than 4x4, the gauss-jordan
            elimination method for inverse is used.
            
            Computing the inverse of a matrix by gauss-jordan elimination involves reducing the matrix to reduced-row echelon form,
            and then performing all the operations used for reduction on the equivalent identity matrix.
            
            This can be proven by the fact that if a matrix is invertible, then the product of its invertible counterpart is equivalent
            to the identity matrix:
                
                        assume matrix X is invertible (x^-1 is the inverse of X):
                        
                        X*(X^-1) = (X^-1)*X = I (the identity matrix)
                        
                        The identity matrix is essentially the matrix equivalent of 1 for matrices:
                            
                        identity property: X*I = X (I is the identity matrix)
                        
                        We can reduce a given matrix to its identity form through row operations (gauss-jordan elimination):
                            
                            [X|I]
                            
                        Then, if we perform the SAME steps we used to reduce matrix X to its identity form, we end up with the inverse:
                            
                            [X|I] -> Row Operations -> [I|X^-1]
    
            Thus, if we reduce a given matrix X to reduced row echelon form, and then perform those vary same steps, in the very same order,
            on an identity matrix of equivalent dimension, we wound up with the inverse of matrix .
            '''
            
            Identity_Matrix = self.Identity()
            for i in range(len(self.matrix)):
                self.matrix[i].extend(Identity_Matrix.matrix[i])
                
            self.EchelonForm("k") # here, we specify choice "k" so that we put our matrix in reduced row echelon form.
            
            Inverse = Matrix([row[len(self.matrix[0])//2:] for row in self.matrix]) # returning the modified identity matrix
            
            return Inverse
            
    def __str__(self):           
        '''
        Provides a custom string representation of the matrix. 
        '''
        string = ""
        for row in self.matrix:
            string+=f"{row}\n"
        return string

def Create():
    '''
    This function is in charge of initializing and creating a user-defined matrix, and calculating the matrix operation
    specified by the user. Two nested functions comprise this function:
        
        - User_Matrix(): in charge of creating the user's matrix,
                         as well as verifying that it is a valid matrix for the specified user operation. For matrix operations that
                         require two matrices (i.e. matrix multiplication, addition, subtraction, etc) this function gets two arguments
                         passed: a tuple, consisting of the FIRST matrix's dimensions (Number of rows, Number of columns), and a
                         string variable, named choice, which is the user's specified choice. Exception handling is used for
                         verifying if the operation is valid and permitted between the two matrices. When creating the first matrix for
                         operations requiring two matrices, None is passed for both arguments, but when creating the second matrix,
                         both arguments are passed to the function.
                         
                         If the matrix operation only requires one matrix (i.e. scalar multiplication, echelon form, etc) then
                         None is passed for the tuple dimensions argument. 
                         
                         
        - Operation_Input(): Provides the user with the various different available matrix operations, and makes the user select
                             their desired operation. This function calls User_Matrix() to create the user's matrix, and then
                             calls the Matrix class to perform the specified matrix operation.
    '''
    
    from math import pi, e
    
    def User_Matrix(FirstMatrixDim: (int, int), choice): # Creates and validates the matrix given by user
        '''
        In charge of creating and initializing the user's matrices, as well as verifying if the operation specified is valid
        for the matrices (or matrix).
        '''
      
        Valid = False
        
        while not Valid:
            try: 
               num_rows = input("Please enter the number of rows of your matrix: ").strip()
               
               if not num_rows.isdigit():
                   raise TypeError("You did not provide a valid numeric entry for the number of rows. Try again.")
                   
               num_rows = int(num_rows)
               
               if num_rows < 1:
                   raise ValueError("You did not provide a valid numeric entry for the number of rows. Try again.")
               
               if (choice=='a') and (FirstMatrixDim[1] != num_rows):
                   
                   raise ValueError("The # of rows of your second matrix does not match the # of columns of your first matrix, hence, "
                                     + f"multiplication cannot be done. Please specify a matrix with {FirstMatrixDim[1]} rows to proceed")
              
               num_columns = input("Please enter the number of columns of you matrix: ").strip()
               
               if not num_columns.isdigit():
                   
                   raise TypeError("You did not provide a valid numeric entry for the number of columns. Try again.")
                   
                   
               num_columns = int(num_columns)
               
               if num_columns < 1:
                   raise ValueError("You did not provide a valid numeric entry for the number of columns. Try again.")
                   
               if (choice=='b' or choice=='c') and (num_columns!=num_rows):
                   
                   raise ValueError('''Matrix Addition/Subtraction can only be performed if the two matrices are of equal dimension.
                                   The dimension of your first and second matrices differ. Please respecify the number of rows and columns
                                   of your second matrix to continue''')
                                   
               if (choice=='i' or choice=='l' or choice=='f') and (num_columns!=num_rows):
                   
                  operation = None
                  
                  match choice:
                      case 'i':
                          operation = 'determinant'
                      case 'l':
                          operation = 'inverse'
                      case 'f':
                          operation = 'power'
                          
                          
                  print(f"The dimensions of your matrix is {num_rows}x{num_columns}")
                  raise ValueError(f"The {operation} can only be computed for square matrices."
                                   + " Please respecify the dimensions of your matrix to continue")
                                       
              
               print(f"The dimensions of your matrix are {num_rows} rows x {num_columns} columns. Is this correct?")
               
               confirmation = input("Enter yes or no: ").strip().lower()
               
               while confirmation != "yes" and confirmation != "no":
                   print("You did not provide a valid response. Try again")
                   confirmation = input("Enter yes or no: ").strip().lower()
                   
               Valid = confirmation=='yes'
               
               if Valid==False:
                   print("Sorry for the mixup! Let's re-analyze your matrix.")
              
            except (ValueError, TypeError) as e:
                print(e)
                Valid = False
        
        matrix = [[] for i in range(num_rows)]
        
        
        for i in range(len(matrix)): # accessing rows of the matrix
            matrix[i] = input(f"Please enter the entries of row {i+1} of your matrix seperated by spaces: ")
            matrix[i] = [eval(num) for num in matrix[i].split(" ") if num!=" " and num!=""]
            while len(matrix[i]) != num_columns:
                print("The number of entries does not match your specified number of columns.")
                matrix[i] = input(f"Please re-enter the entries of row {i+1} of your matrix seperated by spaces: ")
                matrix[i] = [eval(num) for num in matrix[i].split(" ")]
                
        matrix = Matrix(matrix)
        
        return matrix

    def Operation_Input(): # Performs the matrix operation chosen by the user.by calling the Matrix class
        '''
        Calls User_Matrix() to create the user's matrices (or matrix), and performs the specified matrix operation by calling the
        Matrix class.
        '''
        
        from string import ascii_lowercase
        
        print("Welcome to my matrix operations calculator!")
    
        valid = False
        
        while not valid:
            try:
                valid = True
                print("My calculator offers the following operations: ")
                print("(a): Matrix Multiplication\n(b): Matrix Addition\n(c): Matrix Subtraction\n(d): Scalar Multiplication")
                print("(e): Transpose\n(f): Power\n(g): Rank\n(h): Diagonal\n(i): Determinant\n(j): REF\n(k): RREF\n(l): Inverse")
                choice = input("Please select your desired operation: ").lower()
                if choice not in ascii_lowercase[0:12]:
                    raise ValueError("You did not select a valid operation. Try again")
            
            except ValueError as v:
                valid = False
                print(v)
                
        
        if choice=='a' or choice=='b' or choice=='c':
            
            print("let's analyze your 1st matrix")
            
            matrix1 = User_Matrix(None, None)
            
            print("let's analyze your second matrix")
            
            Matrix1Dimensions = (len(matrix1.matrix), len(matrix1.matrix[0]))
            
            matrix2 = User_Matrix(Matrix1Dimensions, choice)
            
            
            match choice:
                
                case "a": print(matrix1*matrix2)
                
                case "b": print(matrix1+matrix2)
                
                case "c": print(matrix1-matrix2)    
                
        else: 
             
            print("Let's analyze your matrix")
            matrix = User_Matrix(None, choice)
            
            if choice=='d':
                
                valid = False
                
                while not valid:
                    try:
                        valid = True
                        scalar = eval(input("Please specify your scalar quantity: "))
                    
                    except (SyntaxError, NameError):
                        valid = False
                        print("You did not enter a valid scalar quantity. Please try again")
                
                print(matrix.Scalar_Multiplication(scalar))
                
            if choice=='e':
        
                print("Your matrix transposed is:")
                print(matrix.Transpose())
                
            if choice=='f':
               
                valid = False
                
                while not valid:
                    try:
                        valid = True
                        power = int(input("Please specify an integer power: "))
                        
                    except (ValueError):
                        valid = False
                        print("You did not specify a valid integer exponent. Try again")
                    
                print("Result: ")
                print(matrix**power)
        
            if choice=='g':
                
                print("Result: ")
                print(matrix.Rank())
                
            if choice=='h':
                
                print("Result: ")
                print(matrix.Diagonal())
                
            if choice=='i':
                
                print("Result: ")
                print(matrix.Determinant())
                
            if choice=='j' or choice=='k':
                
                operation = "Reduced Row Echelon Form" if choice=='k' else "Row Echelon Form"
                print(f"Your Matrix in {operation} is: ")
                print(matrix.EchelonForm(choice))
                
            if choice=='l':
                
                print("The Inverse of your Matrix is: ")
                
                print(matrix.Inverse())
            
    Operation_Input()
            
def main():
    
    Create()

main()












    
