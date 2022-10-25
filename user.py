
import re
import numpy as np
from NewtonRaphson import TrigonometricFunction, PolynomialFunction, System

#
# Load a trignometric function from a string and return a TrigonometricFunction object
# Example : 5 cos(10 x) -8 sin(x) -25.2
#
def interpretTrigonometricFunction(functionString: str) -> TrigonometricFunction:
    functionString = functionString.replace(" ", "")

    regex_cos_coeffs='([+-]{0,1}\d*\.{0,1}\d+)cos'
    cos_coeffs = []
    cos_coeffs=[float(x) for x in re.findall(regex_cos_coeffs, functionString)]

    regex_cos_pos='cos\(([-]{0,1}\d*\d+)x\)'
    cos_pos = []
    cos_pos=[float(x) for x in re.findall(regex_cos_pos, functionString)]

    regex_sin_coeffs='([+-]{0,1}\d*\.{0,1}\d+)sin'
    sin_coeffs = []
    sin_coeffs=[float(x) for x in re.findall(regex_sin_coeffs, functionString)]

    regex_sin_pos='sin\(([-]{0,1}\d*\d+)x\)'
    sin_pos = []
    sin_pos=[float(x) for x in re.findall(regex_sin_pos, functionString)]

    regex_constant='(?<![\(])([+-]{0,1}\d+\.{0,1}\d*)(?![a-zA-Z])' #cos(4x) + 2sin(x) + 4sin(2x) - 2
    constant=[float(x) for x in re.findall(regex_constant, functionString)]

    if constant is not None and len(constant) > 0:
        if len(constant) > 1:
            constant = np.sum(np.array(constant))
        else:
            constant = float(constant[0])
    else:
        constant = 0

    print("cos_coeffs :", cos_coeffs)
    print("cos_pos :", cos_pos)
    print("sin_coeffs :", sin_coeffs)
    print("sin_coeffs :", sin_pos)
    print("const", constant)

    return TrigonometricFunction(cos_coeffs, sin_coeffs, cos_pos, sin_pos, constant)

#
# Load a polynomial function form a string and return a PolynomialFunction object.
# Note : every x term must have a coefficient and a power -> "-1x^1" for -x and "[+]1x^1" for [+]x ([+] means the + is optionnal)
#        spaces arent standard (any space can be added any where)
# 
# Format : a_1 x^b_1 + a_2 x^b_2 + a1 x^b_3 + ... +a_n x^b_n + k
# 
# Example : 1x^120 + 25 + 18 + 3x^2 +2 x^3 +1 x^4 -2.2 x^5 + 10
#
def interpretPolynomialFunction(functionString: str) -> PolynomialFunction:
    functionString = functionString.replace(" ", "")
    
    # get coefficients
    regex_coeff='([+-]{0,1}\d*\.{0,1}\d+)x'
    coeffs=[float(x) for x in re.findall(regex_coeff, functionString)]

    # get exponents
    regex_expo='x\^(\d+)'
    exponents=[int(x) for x in re.findall(regex_expo, functionString)]

    # make the retrieve of constants easier
    operatorsRegex = '([\+\-])'
    splittedEquation = re.split(operatorsRegex, functionString)

    # get all the constants 
    regex_constants="(?<![\^])([+-]{0,1}\d+\.{0,1}\d*)(?![x])"
    constants = [float(x) for x in re.findall(regex_constants, functionString)]

    # sum the constants and add them as a coef of an x^0 term
    if len(constants) > 0:
        if len(constants) > 1:
            coeffs.append(np.sum(np.array(constants)))
        else :
            coeffs += constants
        exponents.append(0)

    return PolynomialFunction(coeffs, exponents)


#
# Create a custom regex for a given variable in the function
#
def variableSelected(variableNb):
    return str(variableNb) + ']\^(\d+)'

#
# Load a system from an array of strings representing the functions of the systems and return a System object
#
def interpretSystem(EquationsString: str) -> System:
    regex_coeff = '([+-]{0,1}\d*\.{0,1}\d+)x'
    regex_constants="(?<![\^\[])([+-]{0,1}\d+\.{0,1}\d*)(?![x\]])" 
    
    system = []
    cons = []
    for i in range(len(EquationsString)):
        func = EquationsString[i] #2x[1]^3*x[2]^2 + 2x[1]^0*x[2]^1
        func = func.replace(" ", "")
        #gets the coeficients of all terms
        coeffs=[float(x) for x in re.findall(regex_coeff, func)]
        
        exp = []
        for j in range(len(EquationsString)):
            var_exp = [int(x) for x in re.findall(variableSelected(j+1), func)]
            exp.append(var_exp)

        eq = []
        for k in range(len(coeffs)):
            term = []
            term.append(coeffs[k])
            for exponant in exp:
                term.append(exponant[k])
            eq.append(term)
        system.append(eq)

        constants = [float(x) for x in re.findall(regex_constants, func)]

        if len(constants) > 0:
            if len(constants) > 1:
                cons.append(np.sum(np.array(constants)))
            else :
                cons += constants

    return System(len(EquationsString), system, cons)

#
# Main function
#
def main():
    print("**********************************")
    print("   📈 Zero of function finder 📈   ")
    print()
    print("    ULiege - Unesp - 2022-2023    ")
    print("**********************************")
    print()

    exit = False
    while(not exit):
        loadChoice = -1         #1 / 2      (resp. in terminal, in a file)
        typeOfFunction = -1     #1 / 2 / 3  (resp. trigonometric function, polynomial function, system)
        toSolve = None          #The function object created from code.py (TrigonometricFuntion, PolynomialFunction, System)
        fileName = ""

        #Choose the loading method (terminal or file)
        while(loadChoice not in [1,2]):
            loadChoice = int(input("Enter a number for the function loading method :\n\t⏩ 1 to load the function/system from the terminal\n\t⏩ 2 to load from a file\n"))
        if loadChoice == 2:
            fileName = input("Enter a file name : ")

        #Choose the type of function
        while(typeOfFunction not in [1,2,3]):
            typeOfFunction = int(input("Enter a number for the function form :\n\t⏩ 1 for a trigonometric function\n\t⏩ 2 for a polynomial function\n\t⏩ 3 for system of functions\n"))

        #Create the function object for the given function type accordingly to the chosen loading method
        if(typeOfFunction == 1):
            if(loadChoice == 1):
                toSolve = interpretTrigonometricFunction(input("Function : "))
            else:
                exit()
                file = open(fileName, "r")
                function = file.readline().strip()
                file.close()
                toSolve = interpretTrigonometricFunction(function)
            
        elif(typeOfFunction == 2):
            if(loadChoice == 1):
                toSolve = interpretPolynomialFunction(input("Function: "))
            else:
                file = open(fileName, "r")
                function = file.readline().strip()
                file.close()
                toSolve = interpretPolynomialFunction(function)

        elif(typeOfFunction == 3):
            if(loadChoice == 1):
                functionsNb = int(input("Number of functions: "))
                functions = []
                for i in range(functionsNb):
                    functions.append(input("Function n•"+str(i+1)+" (e.g. ax[1]^2*x[2]^1...*x[n]^e + ...): "))
                toSolve = interpretSystem(functions)
            else:
                file = open(fileName, "r")
                functions = []
                for function in file.readlines():
                    functions.append(function.strip())
                file.close()
                toSolve = interpretSystem(functions)
        else:
            exit()
        
        #Solve the function
        print("x = ",toSolve.solve(0)," ✅")

        #Ask the user if he want to input a new function
        newFunction = -1
        while(newFunction not in ["Y", "N"]):
            newFunction = input("Do you want to load a new function ? [Y]es or [N]o: ")
            if(newFunction == "N"):
                exit = True

# If the file main.py is the file executed (i.e. not imported is other file which is then executed), 
# we call main()
if __name__ == "__main__":
    main()
