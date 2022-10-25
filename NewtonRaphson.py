import numpy as np
from math import *
from random import *

class Solver():
    def solve(self, tries):
        pass

# ======== Macros ======== #
MAX_ITERATIONS = 10000000
EPSILON = 0.000001

# ======== Trigonometric Functions ======== #
class TrigonometricFunction:
    def __init__(self, cos_coef, sin_coef, cos_pos, sin_pos, cons):
        # For
        # 4sin(x) + sin(3x) -cos(x) + 2cos(2x) + 2cos(3x) - 3 
        # We have :
        self.cos_coef = cos_coef    # np.array([-1,2,2])
        self.cos_pos = cos_pos      # np.array([1,2,3])
        self.sin_coef = sin_coef    # np.array([4,1])
        self.sin_pos = sin_pos      # np.array([1,3])
        self.cons = cons            # -3

    def Calculate_TrigonometricFunction(self, point):
        value = 0.0
        point = float(point)
        for i in range(len(self.cos_coef)):
            # a[1]*cos(1*point) + a[2]*cos(2*point) + ... + a[i]*cos(i*point)
            value += self.cos_coef[i] * cos(self.cos_pos[i]*point)
            
        for i in range(len(self.sin_coef)):
            # b[1]*sin(1*point) + b[2]*sin(2*point) + ... + b[i]*sin(i*point)
            value += self.sin_coef[i] * sin(self.sin_pos[i]*point)
            
        value += self.cons
        return value

    def DerivativeCalculator(self):
            new_sin_coef = []
            new_cos_coef = []

            for i in range(len(self.cos_coef)):
                # a[1]* 1 * (-sin(1*x)) + ... + a[i]* i * (-sin(i*x))
                new_sin_coef.append(-self.cos_coef[i]*self.cos_pos[i])
                
            for i in range(len(self.sin_coef)):
                # a[1]* 1 * cos(1*x) + ... + a[i]* i * cos(i*x)
                new_cos_coef.append(self.sin_coef[i]*self.sin_pos[i])

            # the derivative
            return TrigonometricFunction( np.array(new_cos_coef), np.array(new_sin_coef), self.sin_pos, self.cos_coef,0)

    def solve(self, tries):
        X1 = uniform(-1000,1000)
        derivative = self.DerivativeCalculator()
        Precision = 1.
        starting_point = X1
        X_next = 0.
        iter_nb = 0
        while ( Precision > EPSILON):
            if derivative.Calculate_TrigonometricFunction(X1) != 0 and iter_nb < MAX_ITERATIONS:
                X_next = X1 - (self.Calculate_TrigonometricFunction(X1) /
                            derivative.Calculate_TrigonometricFunction(X1))
                Precision = abs(X1 - X_next)
                X1 = X_next
                iter_nb += 1
            else:
                if tries == 10:
                    return "Newton-Raphson cannot find a solution. This equation might not have any solutions."
                else:
                    if iter_nb >= MAX_ITERATIONS:
                        print("Cannot find a solution with", "{:.3f}".format(starting_point), "as a starting point.\n Trying a new starting point.")
                    else:
                        print("Cannot apply the Newton-Raphson Method (Divison by zero).\n Trying a new starting point.")
                    return self.solve(tries+1)
        return X1

# ======== Polynomial Functions ======== #
class PolynomialFunction:
    def __init__(self, coef, pow):
        # For
        # 3x^3 + 2x + 4
        # We have :
        self.coeficients = coef #np.array([3,2,4])
        self.power = pow        #np.array([3,1,0])

    def Calculate_PolynomialFunction(self, point):
        value = 0.0
        for i in range(len(self.coeficients)):
            value += self.coeficients[i]*(point)**self.power[i]
        return value
    
    def DerivativeCalculator(self):
        new_coef = []
        new_pow = []
        for i in range(len(self.coeficients)):
            # a[1] * 1 * x^(1-1) + a[2] * 2 * x^(2-1) + ... + a[i] * 1 * x^(i-1)
            new_coef.append(self.coeficients[i]*self.power[i])
            if self.power[i] != 0:
                new_pow.append(self.power[i]-1)
            else:
                new_pow.append(self.power[i])
        
        return PolynomialFunction(np.array(new_coef), 
                                        np.array(new_pow))  # the derivative

    def solve(self, tries):
        X1 = uniform(-1000,1000)
        derivative = self.DerivativeCalculator()  # create derivative
        Precision = 1.
        starting_point = X1
        X_next = 0.
        iter_nb = 0
        while (Precision > EPSILON):
            if derivative.Calculate_PolynomialFunction(X1) != 0 and iter_nb < MAX_ITERATIONS:
                X_next = X1 - (self.Calculate_PolynomialFunction(X1) /
                            derivative.Calculate_PolynomialFunction(X1))  # Use Newton-Raphson formula
                Precision = abs(X1-X_next)
                X1 = X_next
                iter_nb += 1
            else:
                if tries == 10:
                    return "Newton-Raphson cannot find a solution. This equation might not have any solutions."
                else:
                    if iter_nb >= MAX_ITERATIONS:
                        print("Cannot find a solution with", "{:.3f}".format(starting_point), "as a starting point.\n Trying a new starting point.")
                    else:
                        print("Cannot apply the Newton-Raphson Method (Divison by zero).\n Trying a new starting point.")
                    return self.solve(tries+1)

        return X_next

# ======== System of Equations ======== #
class System():
    def __init__(self, Equation_Nb, equations, cons):
        self.Equation_Nb = Equation_Nb  # Ex for 2 equations : 2
        
        self.equations = equations      # Ex for y^2 + 3xy + 5x ; 2yx^2 + x + xy :
                                        #   [[[1,0,2],[3,1,1],[5,1,0]], 
                                        #    [[2,2,1],[1,1,0],[1,1,1]]]
        
        # General form => [ [[constant[1], power of var[1], power of var[2], ..., power of var[n]]
        #                    [constant[2],power of var[1], power of var[2], ..., power of var[n]]], <-- Equation 1
        #                   [[cons[1],power[1], power[2], ...], [.    .    .], .................], <-- Equation 2
        #                    ...,
        #                  ]
        
        # Ex for the system y^2 + 3xy + 5x = 3 ; 2yx^2 + x + xy = 2, we have : np.array([3,2])
        self.cons = cons


    def Calcultate_PartialDerivative(self, Eq, var_nb):
        derivative = []
        for i in range(len(Eq)):

            derivative.append(Eq[i].copy())
            derivative[i][0] *= derivative[i][var_nb]
            if(derivative[i][var_nb] != 0):
                derivative[i][var_nb] -= 1

        return derivative

    def Calculate_MultiVariablesEquation_Values(self, eq, vector):
        tot_val = 0
        for i in range(len(eq)): #nb_terme 

            value = eq[i][0] #1

            for j in range(len(eq[i])-1): #0=>2
                value *= vector[j]**eq[i][j+1]

            tot_val += value

        return tot_val

    def Calculate_FxK(self, vector):
        func = []
        for i in range(self.Equation_Nb):
            func.append(self.Calculate_MultiVariablesEquation_Values(self.equations[i], vector))
        func = np.array(func)
        func += self.cons
        return func

    def Jacobian_Matrix(self, vector):
        matrix = []                        # Create an empty matrix

        for i in range(self.Equation_Nb):  # Loops until every equation has been used

            results = []  # Creates an empty list that will receive the eventual results for each equation

            # Loops in order to calculate the partial derivative for each variable
            for j in range(self.Equation_Nb):
                # Calculates the partial derivative for each variable and stocks the intended value in the list "results"
                results.append(self.Calculate_MultiVariablesEquation_Values(self.Calcultate_PartialDerivative(self.equations[i], j+1), vector))
            # Adds the results list which contains the values of the pariatl derivatives for equation[i] in the matrix
            matrix.append(results)

        # returns a 2D matrix(=array) which contains lists of partial derivatives for every equation (in other words, the Jacobian Matrix)
        return np.array(matrix)
    
    def solve(self, tries):
        X1 = [uniform(-1000,1000) for _ in range(self.Equation_Nb)]
        Precision = 1.
        while (Precision > EPSILON):
            if np.linalg.det(self.Jacobian_Matrix(X1)) != 0:
                val = np.linalg.solve(self.Jacobian_Matrix(X1), self.Calculate_FxK(X1))
                X_Next = X1 - val
                Precision = abs(np.mean(X1 - X_Next))
                X1 = X_Next
            else:
                if tries != 10:
                    print("Cannot apply the Newton-Raphson Method (Jacobian Matrix is singular).\n Trying a new starting vector.")
                    return self.Newton_Raphson_System(tries+1)
                else:
                    print("Newton-Raphson cannot find a solution. This equation might not have any solutions.")
        return X1

""""Old tests
poly_func1 = PolynomialFunction(np.array([1, 1, -1]), 
                                np.array([2, 1, 0]))  # x^2 + x - 1 ==> 2x + 1
print("Newton Raphson Polynomial:", poly_func1.Newton_Raphson_Poly(-1,0, 0.000001))

trig_func1 = TrigonometricFunction(np.array([4]), np.array([4, -4]), 
                                   np.array([1]), np.array([1, 2]), 
                                   5) # 4cos(x) + 4sin(x) -4sin(2x) - 5
print("Newton Raphson Trigonometric:",trig_func1.Newton_Raphson_Trigo(-1, 0, 0.000001))

system1 = System(2, [[[1,1,3], [1,2,1] ], [[1,2,0], [1,0,2] ]], [-1,4])

print("Newton Raphson System:", system1.Newton_Raphson_System([0,0],0, 0.000001))
"""