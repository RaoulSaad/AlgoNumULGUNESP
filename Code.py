from re import X
import numpy as np
from math import *

# ======== Trigonometric and Polynomial Functions ======== #
def Create_TrigonometricFunction(cos_coef, sin_coef, cons):
    TrigonometricFunctions = {
        "Cos_Coeficients" : np.array([]), #Ex: np.array([-1,0,2,1]) -> -cos(x) + 2cos(3x) + cos(4x)
        "Sin_Coeficients" : np.array([]), #Ex: np.array([4,0,2]) -> 4sin(x) + 2cos(2x)
        "Constant" : None # Ex: 3 -> -cos(x) + 2cos(3x) + cos(4x) + 4sin(x) + 2cos(2x) = 3
    }
    TrigonometricFunctions["Cos_Coeficients"] = cos_coef
    TrigonometricFunctions["Sin_Coeficients"] = sin_coef
    TrigonometricFunctions["Constant"] = cons
    
    return TrigonometricFunctions

def Calculate_TrigonometricFunction(function, point):
    value = 0
    for i in range(len(function["Cos_Coeficients"])):
        value += function["Cos_Coeficients"][i]*cos((i+1)*point) # a[1]*cos(1*point) + a[2]*cos(2*point) + ... + a[i]*cos(i*point)
    for i in range(len(function["Sin_Coeficients"])):
        value += function["Sin_Coeficients"][i]*sin((i+1)*point) # b[1]*sin(1*point) + b[2]*sin(2*point) + ... + b[i]*sin(i*point)
    value -= function["Constant"]
    return value

def Create_PolynomialFunction(coef, pow):
    PolynomialFunctions = {
        "Coeficients" : np.array([]), #Ex: np.array([1,2,3,4,-5])
        "Power" : np.array([]) # Ex: np.array([0,1,3,5,7]) -> 1 + 2x + 3x^3 + 4x^5 - 5x^7
    }
    PolynomialFunctions["Coeficients"] = coef
    PolynomialFunctions["Power"] = pow
    return PolynomialFunctions

def Calculate_PolynomialFunction(function, point):
    value = 0
    for i in range(len(function["Coeficients"])):
        value += function["Coeficients"][i]*(point)**function["Power"][i]
    return value

def DerivativeCalculator(Function, point, Polynomial = True):
    if Polynomial: #if the function is a polynomial
        new_coef = []
        new_pow = []
        for i in range(1,len(Function["Coeficients"])):
            new_coef.append(Function["Coeficients"][i]*Function["Power"][i]) # a[1] * 1 * x^(1-1) + a[2] * 2 * x^(2-1) + ... + a[i] * 1 * x^(i-1)
            new_pow.append(Function["Power"][i]-1)
        derivative = Create_PolynomialFunction(np.array(new_coef), np.array(new_pow)) #the derivative
        return Calculate_PolynomialFunction(derivative, point) #the value of the derivative on a specific point given
    
    else: #if the function is trigonometric
        new_sin_coef = []
        new_cos_coef = []
        
        for i in range(len(Function["Cos_Coeficients"])):
            new_sin_coef.append(-Function["Cos_Coeficients"][i]*(i+1)) # a[1]* 1 * (-sin(1*x)) + ... + a[i]* i * (-sin(i*x))
        for i in range(len(Function["Sin_Coeficients"])):
            new_cos_coef.append(Function["Sin_Coeficients"][i]*(i+1)) # a[1]* 1 * cos(1*x) + ... + a[i]* i * cos(i*x)
        
        derivative = Create_TrigonometricFunction(np.array(new_cos_coef), np.array(new_sin_coef), 0) #the derivative
        return Calculate_TrigonometricFunction(derivative, point) #the value of the derivative on a specific point given

# ======== System of Equations ======== #
def Create_EqSys(Equation_Nb, equations, cons):
    System = {
        "NumberOfEq" : None, #Ex: 2 -> 2 equations
        "Equations" : [], #Ex: [ [[1,0,2],[3,1,1],[5,1,0]], [[2,2,1],[1,1,0],[1,1,1]] ] -> Ex: y^2 + 3xy + 5x ; 2yx^2 + x + xy
                                    #General form => [ [[constant[1], power of var[1], power of var[2], ..., power of var[n]]
                                    #                   [constant[2],power of var[1], power of var[2], ..., power of var[n]]], <-- Equation 1
                                    #                  [[cons[1],power[1], power[2], ...], [.    .    .], .................], <-- Equation 2
                                    #                    .........,
                                    #                    [.....]]
        "Constant" : np.array([]) #Ex: np.array([3,2]) -> y^2 + 3xy + 5x = 3 ; 2yx^2 + x + xy = 2
    }
    System["NumberOfEq"] = Equation_Nb
    System["Equations"] = equations
    System["Constant"] = cons
    return System

def PartialDerivative_Calculator(Eq, var_nb):
    derivative = []
    
    for i in range(len(Eq)): # [(1,0,2),(3,1,1), (5,1,0)] => [(0,-1,2), (3,0,1), (5,0,0)]
    
        derivative.append(Eq[i].copy())
        derivative[i][0] *= derivative[i][var_nb]
        derivative[i][var_nb] -= 1
    
    return derivative

def Calculate_PartialDerivative(partial_deriv, vector): #calculates the value of a partial derivative 
                                                        #after replacing each variable with the intended number
    tot_val = 0
    
    for i in range(len(partial_deriv)):
    
        value = partial_deriv[i][0]
    
        for j in range(len(partial_deriv[i])-1):
            value*= vector[j]**partial_deriv[i][j+1]
    
        tot_val += value
    
    return tot_val

def Jacobian_Matrix(EqSys, vector):
    matrix = []                              # Create an empty matrix
    
    for i in range(EqSys["NumberOfEq"]):     #Loops until every equation was used
    
        results = []                         #Creates an empty list that will receive the eventual results for each equation
        
        for j in range(EqSys["NumberOfEq"]): #Loops in order to calculate the partial derivative for each variable
            results.append(Calculate_PartialDerivative(PartialDerivative_Calculator(EqSys["Equations"][i], j+1), vector)) #Calculates the partial derivative for each variable and stocks the intended value in the list "results"
        matrix.append(results)               #Adds the results list which contains the values of the pariatl derivatives for equation[i] in the matrix
    
    return np.array(matrix)                  #returns a 2D matrix(=array) which contains lists of partial derivatives for every equation (in other words, the Jacobian Matrix)




# ===== Some tests ===== #

poly_func1 = Create_PolynomialFunction(np.array([1,2,2]), np.array([0,1,3])) # 1 + 2x + 2x^3
print(poly_func1["Coeficients"])
print(poly_func1["Power"])
poly_func1_prime = DerivativeCalculator(poly_func1, 2) # 2 + 6*4 = 26
print(poly_func1_prime)

trig_func1 = Create_TrigonometricFunction(np.array([4]), np.array([4, -4]), 5)
print("Coef of cos:", trig_func1["Cos_Coeficients"])
print("Coef of sin:", trig_func1["Sin_Coeficients"])
trig_func1_prime = DerivativeCalculator(trig_func1, 1, False)
print(trig_func1_prime)

system1 = Create_EqSys(2, [ [[1,0,2], [2,2,0], [-3,1,1]], [[2,2,2],[3,1,0],[-2,0,1],[1,1,1]] ], np.array([3,4]))
print('\n'.join(['\t'.join([str(cell) for cell in row]) for row in system1["Equations"]]))
print(Jacobian_Matrix(system1, [2,2]))