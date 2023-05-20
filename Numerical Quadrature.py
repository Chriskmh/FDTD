from cProfile import label
import numpy as np
import matplotlib.pyplot as plt
def p0(x):
    return 1

def p1(x):
    return x

def p2(x):
    return 0.5*(3*(np.power(x,2))-1)

def p3(x):
    return 0.5*(5*(np.power(x,3))-3*x)

def p4(x):
    return 1/8*(35*np.power(x,4)-30*(x*x)+3)

# anti_derivative
def  ad_p0(x):
    return x+1

def  ad_p1(x):
    return 0.5*x*x-0.5

def  ad_p2(x):
    return 0.5*x*(x*x-1)

def  ad_p3(x):
    return 0.5*x*x*((5/4)*x*x-(3/2))+1/8

def  ad_p4(x):
    return x/8*(7*np.power(x,4)-10*x*x+3)

def analytical_anti_derivative_p(w,x):
    lenx = len(x)
    y = np.zeros(lenx)
    for i in range(lenx):
        ad_p = np.array([ad_p0(x[i]),ad_p1(x[i]),ad_p2(x[i]),ad_p3(x[i]),ad_p4(x[i])]) 
        y[i] = w@ad_p

    return y

def eval_legendre(w,x):
    lenx = len(x)
    y = np.zeros(lenx)
    for i in range(lenx):
        p = np.array([p0(x[i]),p1(x[i]),p2(x[i]),p3(x[i]),p4(x[i])])
        y[i] = w@p
    
    return y

def left_riemann_sum(y,x):
    deltax = x[1]-x[0]
    sum_y = np.zeros(len(y)-1)
    for i in range(len(sum_y)):
        sum_y[i] = np.sum(y[0:i+1])

    sum_y = deltax*sum_y
    return sum_y

def right_riemann_sum(y,x):
    deltax = x[1]-x[0]
    sum_y = np.zeros(len(y)-1)
    for i in range(len(sum_y)):
        sum_y[i] = np.sum(y[1:i])

    sum_y = deltax*sum_y
    return sum_y

def  trapezoid_rule(y,x):
    deltax = x[1]-x[0]
    temp = np.zeros(len(y)-1)
    sum_y = np.zeros(len(y)-1)
    for i in range(len(temp)):
        temp[i] = y[i]+y[i+1]

    for i in range(len(sum_y)):
        sum_y[i] = np.sum(temp[0:i])

    sum_y = sum_y*deltax/2
    return sum_y


def problem_a():
    x = np.arange(-1,1,0.05)

    plt.figure()
    plt.plot(x,np.ones(len(x)),label = 'p0')
    plt.plot(x,p1(x),label = 'p1')
    plt.plot(x,p2(x),label = 'p2')
    plt.plot(x,p3(x),label = 'p3')
    plt.plot(x,p4(x),label = 'p4')

    plt.legend()
    plt.show()


# problem_a()


def problem_b_and_c_and_d():
    
    # generate random vectors w0-w4
    w = np.random.randn(25).reshape(5,5)

    # calculate eval_legendre with w0-w4
    x = np.arange(-1,1,0.05)
    y = []
    for i in range(5):
        wi = w[i]
        yi = eval_legendre(wi,x)
        y.append(yi)

    y = np.array(y)

    # Plot item c)
    plt.figure()
    plt.title('item c): eval_legendre')
    for i in range(5):
        plt.plot(x,y[i],label = 'y'+str(i))

    plt.legend()
    
    # Calculate and plot numerical integral of function eval_legendre with w0-w4
    plt.figure()
    plt.title('item d) : Numerical integral_y')
    integral_y = []
    for i in range(5):
        integral_y.append(trapezoid_rule(y[i],x))
        plt.plot(x[0:len(x)-1],integral_y[i],label = 'integral_y'+str(i))
    plt.legend()

    # Calculate and plot analytical integral Legendre polynomials
    plt.figure()
    plt.title('item d) :analytical integral_y ')
    analytical_integral_y = []
    for i in range(5):
        wi = w[i]
        temp = analytical_anti_derivative_p(wi,x)
        analytical_integral_y.append(temp)

    analytical_integral_y = np.array(analytical_integral_y)
    for i in range(5):
        plt.plot(x,analytical_integral_y[i],label = 'analytical_integral_y'+str(i))

    plt.legend()
    plt.show()


# problem_b_and_c_and_d()



def problem_e():
# Define relative_error
    def relative_error(Iquad, Iexact):
        n = len(Iquad)
        error = 0
        for i in range(n):
            if i == 0:
                continue
            error = error + np.abs((Iquad[i]-Iexact[i])/Iexact[i])
        return error/n

    error_left_sum = []
    error_right_sum = []
    error_trapezoid = []

    N = [10,20,50,100,200]
    for n in N:
        # Generate 1000r andom  vector w
        w = np.random.randn(5000).reshape(1000,5)
        x = np.linspace(-1,1,n+1)
        lenx = len(x)

# Calculate eval_legendre function
        y = []
        for i in range(1000):
            wi = w[i]
            yi = eval_legendre(wi,x)
            y.append(yi)

        y = np.array(y)

# Claculate numerical integral
        left_riemann_integral = []
        right_riemann_integral = []
        trapezoid_rule_integral = []

        for i in range(1000):
            left_riemann_integral.append(left_riemann_sum(y[i],x))
            right_riemann_integral.append(right_riemann_sum(y[i],x))
            trapezoid_rule_integral.append(trapezoid_rule(y[i],x))

        left_riemann_integral = np.array(left_riemann_integral)
        right_riemann_integral = np.array(right_riemann_integral)
        trapezoid_rule_integral = np.array(trapezoid_rule_integral)


        # Calculate analytical_integral

        analytical_integral = []
        for i in range(1000):
            wi = w[i]
            temp = analytical_anti_derivative_p(wi,x)
            analytical_integral.append(temp)
        
        analytical_integral = np.array(analytical_integral)
        analytical_integral = analytical_integral[:,:lenx-1]

        left_riemann_error = 0
        right_riemann_error = 0
        trapezoid_rule_error = 0
        # midpoint_rule_error = 0

# Calculate error
        for i in range(lenx):
            left_riemann_error = left_riemann_error + relative_error(left_riemann_integral[i],analytical_integral[i])
            right_riemann_error = right_riemann_error + relative_error(right_riemann_integral[i] , analytical_integral[i])
            trapezoid_rule_error = trapezoid_rule_error + relative_error(trapezoid_rule_integral[i],analytical_integral[i])
            # midpoint_rule_error = midpoint_rule_error + relative_error(midpoint_rule_integral[i],analytical_integral[i])

        left_riemann_error = np.log10(left_riemann_error)
        right_riemann_error = np.log10(right_riemann_error)
        trapezoid_rule_error = np.log10(trapezoid_rule_error)
        # midpoint_rule_error = np.log10(midpoint_rule_error)

        # print(left_riemann_error)
        # print(right_riemann_error)
        # print(trapezoid_rule_error)

        # print('__________')

        error_left_sum.append(left_riemann_error)
        error_right_sum.append(right_riemann_error)
        error_trapezoid.append(trapezoid_rule_error)

# Plot error vs N on log log plot 

    N = np.array(N)
    N = np.log10(N)

    plt.figure()
    plt.scatter(N,error_right_sum,label = 'error_right_sum')
    plt.scatter(N,error_left_sum,label = 'error_left_sum')
    plt.scatter(N,error_trapezoid,label = 'error_trapezoid')

    logn = np.arange(1,2.5,0.05)

# linear regression

    z1 = np.polyfit(N, error_right_sum, 1)  
    p1 = np.poly1d(z1)
    log1 = p1(logn)
    plt.plot(logn,log1,label = ' linear regression error_right_sum')

    z2 = np.polyfit(N, error_left_sum, 1)  
    p2 = np.poly1d(z2)
    log2 = p2(logn)
    plt.plot(logn,log2,label = ' linear regression error_left_sum')

    z3 = np.polyfit(N, error_trapezoid, 1)  
    p3 = np.poly1d(z3)
    log3 = p3(logn)
    plt.plot(logn,log3,label = ' linear regression error_trapezoid')

    plt.legend()
    plt.show()

problem_e()
