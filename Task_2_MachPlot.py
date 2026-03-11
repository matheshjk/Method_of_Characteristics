import numpy as np
from sympy import Eq, solve, symbols
import matplotlib.pyplot as plt
import matplotlib.tri as tri

def prandtl_meyer_function(M, gamma): #find nu from M
    nu = (np.sqrt((gamma+1)/(gamma-1)))*(np.arctan(np.sqrt(((gamma-1)/(gamma+1))*(M**2 - 1)))) - np.arctan(np.sqrt(M**2 - 1))
    return nu

def isentropic_mach(p1, p2, M1, gamma): #find M2 from isentropic pressure ratio and M1
    M2 = np.sqrt((((1+((gamma-1)/2)*M1**2)*(p1/p2)**((gamma-1)/gamma))-1)*(2/(gamma-1)))
    return M2

def isentropic_pressure_from_M(p1, M1, M2, gamma):
    p2 = p1*((1+((gamma - 1)*0.5*M1**2))/(1+((gamma - 1)*0.5*M2**2)))**(gamma/(gamma -1))
    return p2

def Mach_angle(M): #find mu from M
    mu = np.arcsin(1/M)
    return mu
def intersection(x1,y1,a1,x2,y2,a2):
    #if a1*a2 > 0:
        ##exit()
    #else:
    x,y = symbols('x y')
    line1 = Eq(y-y1, a1*(x-x1))
    line2 = Eq(y-y2, a2*(x-x2))
    intersection = solve((line1,line2),(x,y))
    return intersection[x], intersection[y] 
def M_from_nu(nu, gamma):
    a = 1
    b = 100
    while ((b-a) >= 0.01):
        c = (a+b)/2
        eqn_a = prandtl_meyer_function(a, gamma) - nu
        eqn_b = prandtl_meyer_function(b,gamma) - nu
        eqn_c = prandtl_meyer_function(c,gamma) -nu
    
        if (eqn_a*eqn_c < 0):
            b = c
        elif(eqn_b*eqn_c < 0):
            a = c
    return c


def gamma_minus_reflections_initial(nu_arr_m,phi_arr_m):    #initial expansion fan from nozzle exit
    
    #initialize arrays
    nu_mat = np.zeros((n_chr,n_chr))
    mu_mat = np.zeros((n_chr,n_chr))
    phi_mat = np.zeros((n_chr,n_chr))
    M_mat = np.zeros((n_chr,n_chr))
    P_mat = np.zeros((n_chr,n_chr))
    x_mat = np.zeros((n_chr,n_chr))
    y_mat = np.zeros((n_chr,n_chr))
    slope_m_mat = np.zeros((n_chr,n_chr))
    slope_p_mat = np.zeros((n_chr,n_chr))
    
    
    #main loop
    for j in range(0,n_chr):
        for i in range(0,n_chr):
           
            if j == 0:
                if i == 0:  #(0,0)
                    nu_mat[i][j] = nu_arr_m[i] 
                    phi_mat[i][j] = phi_arr_m[i]
                    M_mat[i][j] = M_from_nu(nu_mat[i][j], gamma)
                    mu_mat[i][j] = Mach_angle(M_mat[i][j])
                    P_mat[i][j] = isentropic_pressure_from_M(Pe, Me, M_mat[i][j], gamma)
                    
                    slope_m_mat[i][j] = np.tan(phi_mat[i][j]- mu_mat[i][j])
                    slope_p_mat[i][j] = np.tan(phi_mat[i][j]+ mu_mat[i][j])
                    
                    coordinates_1 = intersection(x0, y0, slope_m_mat[i][j], 0 , 0 , 0 )
                    
                    x_mat[i][j] = coordinates_1[0]
                    y_mat[i][j] = coordinates_1[1]
                    x_values = [x0, x_mat[i][j]]
                    y_values = [y0, y_mat[i][j]]
                    
                    plt.plot(x_values, y_values, c = 'black')
                else:   #(i,0)
                    phi_mat[i][j] = phi_arr_m[i]
                    nu_mat[i][j] = nu_arr_m[i]
                    M_mat[i][j] = M_from_nu(nu_mat[i][j], gamma)
                    P_mat[i][j] = isentropic_pressure_from_M(Pe, Me, M_mat[i][j], gamma)
                    mu_mat[i][j] = Mach_angle(M_mat[i][j])
                    
                    slope_m_mat[i][j] = np.tan(phi_mat[i][j] - mu_mat[i][j])
                    slope_p_mat[i][j] = np.tan(phi_mat[i][j] + mu_mat[i][j])
                    
                    coordinates_2 = intersection(x0, y0,slope_m_mat[i][j],x_mat[i-1][j], y_mat[i-1][j], slope_p_mat[i-1][j]) 
                    x_mat[i][j] = coordinates_2[0]
                    y_mat[i][j] = coordinates_2[1]
                    
                    x_values = [x0, x_mat[i][j]]
                    x2_values = [x_mat[i-1][j],x_mat[i][j]]
                    y_values = [y0, y_mat[i][j]]
                    y2_values = [y_mat[i-1][j],y_mat[i][j]]
                    
                    plt.plot(x_values, y_values, color = 'black')
                    plt.plot(x2_values,y2_values,color = 'black')
                    
            if  j !=0 and i != 0:
                if i == j:  #(i,i)
                    phi_mat[i][j] = phi_e
                    nu_mat[i][j] = nu_mat[i][j-1] + phi_mat[i][j-1]
                    M_mat[i][j] = M_from_nu(nu_mat[i][j], gamma)
                    P_mat[i][j] = isentropic_pressure_from_M(Pe, Me, M_mat[i][j], gamma)
                    mu_mat[i][j] = Mach_angle(M_mat[i][j])

                    slope_m_mat[i][j] = np.tan(phi_mat[i][j] - mu_mat[i][j])
                    slope_p_mat[i][j] = np.tan(phi_mat[i][j] + mu_mat[i][j])
                    
                    coordinates_2 = intersection(0, 0, 0,x_mat[i][j-1], y_mat[i][j-1], slope_m_mat[i][j-1]) 
                    x_mat[i][j] = coordinates_2[0]
                    y_mat[i][j] = coordinates_2[1]
                    
                  
                    x2_values = [x_mat[i][j-1],x_mat[i][j]]
                    
                    y2_values = [y_mat[i][j-1],y_mat[i][j]]
                    
                    
                    plt.plot(x2_values,y2_values,color = 'black')
                    
                elif i > j:     #rest points
                    phi_mat[i][j] = 0.5*(-nu_mat[i-1][j] + phi_mat[i-1][j] + nu_mat[i][j-1] + phi_mat[i][j-1])
                    nu_mat[i][j] = 0.5*(nu_mat[i-1][j] - phi_mat[i-1][j] + nu_mat[i][j-1] + phi_mat[i][j-1])
                    M_mat[i][j] = M_from_nu(nu_mat[i][j], gamma)
                    P_mat[i][j] = isentropic_pressure_from_M(Pe, Me, M_mat[i][j], gamma)
                    mu_mat[i][j] = Mach_angle(M_mat[i][j])
                    
                    slope_m_mat[i][j] = np.tan(phi_mat[i][j] - mu_mat[i][j])
                    slope_p_mat[i][j] = np.tan(phi_mat[i][j] + mu_mat[i][j])
                    
                    coordinates_2 = intersection(x_mat[i][j-1],y_mat[i][j-1],slope_m_mat[i][j],x_mat[i-1][j], y_mat[i-1][j], slope_p_mat[i-1][j]) 
                    x_mat[i][j] = coordinates_2[0]
                    y_mat[i][j] = coordinates_2[1]
                    
                    x_values = [x_mat[i][j-1], x_mat[i][j]]
                    x2_values = [x_mat[i-1][j],x_mat[i][j]]
                    y_values = [y_mat[i][j-1], y_mat[i][j]]
                    y2_values = [y_mat[i-1][j],y_mat[i][j]]
                    
                    plt.plot(x_values, y_values, color = 'black')
                    plt.plot(x2_values,y2_values,color = 'black')
                    
                else: #skip instances where j>i
                    continue   
    
    #initialise arrays for returning and plotting
    nu_end_arr = []
    phi_end_arr = []
    x_end_arr = []
    y_end_arr = []
    M_end_arr =[]
    M_plot = []
    x_plot = []
    y_plot = []
    x_start = [0]
    y_start = [1]
    M_start = [M1]
               
    for j in range(0,n_chr):    #loop to find array values to return and plot
        nu_end_arr = np.append(nu_end_arr, nu_mat[n_chr-1][j])
        phi_end_arr = np.append(phi_end_arr, phi_mat[n_chr-1][j])
        x_end_arr = np.append(x_end_arr, x_mat[n_chr-1][j])
        y_end_arr = np.append(y_end_arr, y_mat[n_chr-1][j])
        M_end_arr = np.append(M_end_arr, M_mat[n_chr-1][j])
        for i in range(0,n_chr):
            if i >= j:
                M_plot = np.append(M_plot, M_mat[i][j])  
                x_plot = np.append(x_plot, x_mat[i][j])
                y_plot = np.append(y_plot, y_mat[i][j])     
    for i in range (0,n_chr): #loop for plotting M/P in the initial expansion fan
        x_start = np.append(x_start, x_mat[i][0])
        y_start = np.append(y_start, y_mat[i][0])
        M_start = np.append(M_start, M_mat[i][0])
    
    
    triang = tri.Triangulation(x_plot,y_plot)   #plot values in reflection region
    contour = plt.tricontourf(triang, M_plot, levels=50, cmap='jet',vmin = 2, vmax = 3)  
    plt.clim(2, 3)
    
    triang_Machinitial = tri.Triangulation(x_start, y_start) #plot values in the expansion fan ahead
    plt.tricontourf(triang_Machinitial, M_start, levels = 50, cmap = 'jet', vmin = 2, vmax = 3)
    # Add a colorbar
    cbar = plt.colorbar(contour, extend = 'both')   #legend colorbar
    
    
    return nu_end_arr, phi_end_arr ,x_end_arr, y_end_arr, M_end_arr
    
def gamma_plus_reflections(nu_arr_p,phi_arr_p, x_arr, y_arr,M_arr):
    
    #initialise the arrays
    nu_mat = np.zeros((n_chr,n_chr))
    mu_mat = np.zeros((n_chr,n_chr))
    phi_mat = np.zeros((n_chr,n_chr))
    M_mat = np.zeros((n_chr,n_chr))
    P_mat = np.zeros((n_chr,n_chr))
    x_mat = np.zeros((n_chr,n_chr))
    y_mat = np.zeros((n_chr,n_chr))
    slope_m_mat = np.zeros((n_chr,n_chr))
    slope_p_mat = np.zeros((n_chr,n_chr))
   
    #main loop
    for j in range(0,n_chr):
        for i in range(0,n_chr):
           
            if i == 0 :
                if j == 0:  #(0,0)
                    nu_mat[i][j] = nu_1
                    phi_mat[i][j] = phi_arr_p[i]
                    M_mat[i][j] = M_from_nu(nu_mat[i][j], gamma)
                    mu_mat[i][j] = Mach_angle(M_mat[i][j])
                    P_mat[i][j] = isentropic_pressure_from_M(Pe, Me, M_mat[i][j], gamma)
                    
                    slope_m_mat[i][j] = np.tan(phi_mat[i][j]- mu_mat[i][j])
                    slope_p_mat[i][j] = np.tan(phi_mat[i][j]+ mu_mat[i][j])
                    
                    coordinates_1 = intersection(x_arr[j], y_arr[j], slope_p_mat[i][j], x0 , y0 , np.tan(phi_mat[i][j]) )
                    
                    x_mat[i][j] = coordinates_1[0]
                    y_mat[i][j] = coordinates_1[1]
                    x_values = [x_arr[j], x_mat[i][j]]
                    y_values = [y_arr[j], y_mat[i][j]]
                    x_boundary = [x0,x_mat[i][j] ]
                    y_boundary = [y0,y_mat[i][j] ]
                    
                    plt.plot(x_values, y_values, c = 'black')
                    plt.plot(x_boundary, y_boundary, c = 'red')
                else:   #(0,j)
                    phi_mat[i][j] = 0.5*(nu_mat[i][j-1]+ phi_mat[i][j-1]- nu_arr_p[j] + phi_arr_p[j])
                    nu_mat[i][j] = 0.5*(nu_mat[i][j-1]+ phi_mat[i][j-1]+ nu_arr_p[j] - phi_arr_p[j])
                    M_mat[i][j] = M_from_nu(nu_mat[i][j], gamma)
                    P_mat[i][j] = isentropic_pressure_from_M(Pe, Me, M_mat[i][j], gamma)
                    mu_mat[i][j] = Mach_angle(M_mat[i][j])
                    
                    slope_m_mat[i][j] = np.tan(phi_mat[i][j] - mu_mat[i][j])
                    slope_p_mat[i][j] = np.tan(phi_mat[i][j] + mu_mat[i][j])
                    
                    coordinates_2 = intersection(x_arr[j], y_arr[j],slope_p_mat[i][j],x_mat[i][j-1], y_mat[i][j-1], slope_m_mat[i][j]) 
                    x_mat[i][j] = coordinates_2[0]
                    y_mat[i][j] = coordinates_2[1]
                    
                    x_values = [x_arr[j], x_mat[i][j]]
                    x2_values = [x_mat[i][j-1],x_mat[i][j]]
                    y_values = [y_arr[j], y_mat[i][j]]
                    y2_values = [y_mat[i][j-1],y_mat[i][j]]
                    
                    
                    
                    plt.plot(x_values, y_values, color = 'black')
                    plt.plot(x2_values,y2_values,color = 'black')
                    
            if  i != 0 and j !=0:
                if j > i: #(j>i) is before j==i because j==i points use data from j>i points
                    phi_mat[i][j] = 0.5*(-nu_mat[i-1][j] + phi_mat[i-1][j] + nu_mat[i][j-1] + phi_mat[i][j-1])
                    nu_mat[i][j] = 0.5*(nu_mat[i-1][j] - phi_mat[i-1][j] + nu_mat[i][j-1] + phi_mat[i][j-1])
                    M_mat[i][j] = M_from_nu(nu_mat[i][j], gamma)
                    P_mat[i][j] = isentropic_pressure_from_M(Pe, Me, M_mat[i][j], gamma)
                    mu_mat[i][j] = Mach_angle(M_mat[i][j])
                    
                    slope_m_mat[i][j] = np.tan(phi_mat[i][j] - mu_mat[i][j])
                    slope_p_mat[i][j] = np.tan(phi_mat[i][j] + mu_mat[i][j])
                    
                    coordinates_2 = intersection(x_mat[i][j-1],y_mat[i][j-1],slope_m_mat[i][j],x_mat[i-1][j], y_mat[i-1][j], slope_p_mat[i-1][j]) 
                    x_mat[i][j] = coordinates_2[0]
                    y_mat[i][j] = coordinates_2[1]
                    
                    x_values = [x_mat[i][j-1], x_mat[i][j]]
                    x2_values = [x_mat[i-1][j],x_mat[i][j]]
                    y_values = [y_mat[i][j-1], y_mat[i][j]]
                    y2_values = [y_mat[i-1][j],y_mat[i][j]]
                    
                    plt.plot(x_values, y_values, color = 'black')
                    plt.plot(x2_values,y2_values,color = 'black')
                    
                elif i ==j: #(i,i)
                    nu_mat[i][j] = nu_1 
                    phi_mat[i][j] = nu_mat[i][j] - nu_mat[i-1][j] + phi_mat[i-1][j]
                    
                    M_mat[i][j] = M_from_nu(nu_mat[i][j], gamma)
                    P_mat[i][j] = isentropic_pressure_from_M(Pe, Me, M_mat[i][j], gamma)
                    mu_mat[i][j] = Mach_angle(M_mat[i][j])

                    slope_m_mat[i][j] = np.tan(phi_mat[i][j] - mu_mat[i][j])
                    slope_p_mat[i][j] = np.tan(phi_mat[i][j] + mu_mat[i][j])
                    
                    coordinates_2 = intersection(x_mat[i-1][j-1], y_mat[i-1][j-1], np.tan(phi_mat[i][j]),x_mat[i-1][j], y_mat[i-1][j], slope_p_mat[i][j]) 
                    x_mat[i][j] = coordinates_2[0]
                    y_mat[i][j] = coordinates_2[1]
                    
                    x_boundary = [x_mat[i-1][j-1], x_mat[i][j]]
                    x_values = [x_mat[i-1][j],x_mat[i][j]]
                    y_boundary = [y_mat[i-1][j-1], y_mat[i][j]]
                    y_values = [y_mat[i-1][j],y_mat[i][j]]
                    
                    plt.plot(x_boundary, y_boundary, color = 'red')
                    plt.plot(x_values,y_values,color = 'black')
                else:
                    continue    
    #initialise array for plotting            
    nu_end_arr = []
    phi_end_arr = []
    nu_end_arr = []
    phi_end_arr = []
    x_end_arr = []
    y_end_arr = []
    M_end_arr =[]
    M_plot = []
    x_plot = []
    y_plot = []
    x_start = x_arr
    y_start = y_arr
    M_start = M_arr    
    
    
    
    
    
    for i in range(0,n_chr):
        nu_end_arr = np.append(nu_end_arr, nu_mat[i][n_chr-1])
        phi_end_arr = np.append(phi_end_arr, phi_mat[i][n_chr-1])
        x_end_arr = np.append(x_end_arr, x_mat[i][n_chr-1])
        y_end_arr = np.append(y_end_arr, y_mat[i][n_chr-1]) 
        M_end_arr = np.append(M_end_arr, M_mat[i][n_chr-1]) 
        for j in range(0,n_chr):
            if j >= i:
                M_plot = np.append(M_plot, M_mat[i][j])  
                x_plot = np.append(x_plot, x_mat[i][j])
                y_plot = np.append(y_plot, y_mat[i][j])   
                
    for j in range (0,n_chr):
        x_start = np.append(x_start, x_mat[0][j])
        y_start = np.append(y_start, y_mat[0][j])
        M_start = np.append(M_start, M_mat[0][j])  
    
    #for plotting values for region with P = P_atm
    x_1region = [0, x_arr[0], x_start[n_chr]]
    y_1region = [1, y_arr[0], y_start[n_chr]]
    M_1region = [M1, M_arr[0], M_start[n_chr]]
    
    triang = tri.Triangulation(x_plot,y_plot)
    contour = plt.tricontourf(triang, M_plot, levels=50, cmap='jet',vmin = 2, vmax = 3)  # Adjust levels and cmap as needed
    plt.clim(2, 3)
    
    triang_Machinitial = tri.Triangulation(x_start, y_start)
    plt.tricontourf(triang_Machinitial, M_start, levels = 50, cmap = 'jet', vmin = 2, vmax = 3)
    
    triang_constant = tri.Triangulation(x_1region, y_1region)   #region P_atm
    plt.tricontourf(triang_constant, M_1region, levels = 50, cmap = 'jet', vmin = 2, vmax = 3)
    
    return nu_end_arr, phi_end_arr ,x_end_arr, y_end_arr, M_end_arr    

def gamma_minus_reflections(nu_arr_m,phi_arr_m, x_arr, y_arr,M_arr,x_end,y_end, M_end):
    
    nu_mat = np.zeros((n_chr,n_chr))
    mu_mat = np.zeros((n_chr,n_chr))
    phi_mat = np.zeros((n_chr,n_chr))
    M_mat = np.zeros((n_chr,n_chr))
    P_mat = np.zeros((n_chr,n_chr))
    x_mat = np.zeros((n_chr,n_chr))
    y_mat = np.zeros((n_chr,n_chr))
    slope_m_mat = np.zeros((n_chr,n_chr))
    slope_p_mat = np.zeros((n_chr,n_chr))
    
    
    for j in range(0,n_chr):
        for i in range(0,n_chr):
            
            if j == 0:
                if i == 0:
                    nu_mat[i][j] = nu_arr_m[i] 
                    phi_mat[i][j] = phi_arr_m[i]
                    M_mat[i][j] = M_from_nu(nu_mat[i][j], gamma)
                    mu_mat[i][j] = Mach_angle(M_mat[i][j])
                    P_mat[i][j] = isentropic_pressure_from_M(Pe, Me, M_mat[i][j], gamma)
                    
                    slope_m_mat[i][j] = np.tan(phi_mat[i][j]- mu_mat[i][j])
                    slope_p_mat[i][j] = np.tan(phi_mat[i][j]+ mu_mat[i][j])
                    
                    coordinates_1 = intersection(x_arr[i], y_arr[i], slope_m_mat[i][j], 0 , 0 , 0 )
                    
                    x_mat[i][j] = coordinates_1[0]
                    y_mat[i][j] = coordinates_1[1]
                    x_values = [x_arr[i], x_mat[i][j]]
                    y_values = [y_arr[i], y_mat[i][j]]
                    
                    plt.plot(x_values, y_values, c = 'black')
                else:
                    phi_mat[i][j] = phi_arr_m[i]
                    nu_mat[i][j] = nu_arr_m[i]
                    M_mat[i][j] = M_from_nu(nu_mat[i][j], gamma)
                    P_mat[i][j] = isentropic_pressure_from_M(Pe, Me, M_mat[i][j], gamma)
                    mu_mat[i][j] = Mach_angle(M_mat[i][j])
                    
                    slope_m_mat[i][j] = np.tan(phi_mat[i][j] - mu_mat[i][j])
                    slope_p_mat[i][j] = np.tan(phi_mat[i][j] + mu_mat[i][j])
                    
                    coordinates_2 = intersection(x_arr[i], y_arr[i],slope_m_mat[i][j],x_mat[i-1][j], y_mat[i-1][j], slope_p_mat[i-1][j]) 
                    x_mat[i][j] = coordinates_2[0]
                    y_mat[i][j] = coordinates_2[1]
                    
                    x_values = [x_arr[i], x_mat[i][j]]
                    x2_values = [x_mat[i-1][j],x_mat[i][j]]
                    y_values = [y_arr[i], y_mat[i][j]]
                    y2_values = [y_mat[i-1][j],y_mat[i][j]]
                    
                    plt.plot(x_values, y_values, color = 'black')
                    plt.plot(x2_values,y2_values,color = 'black')
                    
            if  j !=0 and i != 0:
                if i == j:
                    phi_mat[i][j] = phi_e
                    nu_mat[i][j] = nu_mat[i][j-1] + phi_mat[i][j-1]
                    M_mat[i][j] = M_from_nu(nu_mat[i][j], gamma)
                    P_mat[i][j] = isentropic_pressure_from_M(Pe, Me, M_mat[i][j], gamma)
                    mu_mat[i][j] = Mach_angle(M_mat[i][j])

                    slope_m_mat[i][j] = np.tan(phi_mat[i][j] - mu_mat[i][j])
                    slope_p_mat[i][j] = np.tan(phi_mat[i][j] + mu_mat[i][j])
                    
                    coordinates_2 = intersection(0, 0, 0,x_mat[i][j-1], y_mat[i][j-1], slope_m_mat[i][j-1]) 
                    x_mat[i][j] = coordinates_2[0]
                    y_mat[i][j] = coordinates_2[1]
                    
                  
                    x2_values = [x_mat[i][j-1],x_mat[i][j]]
                    
                    y2_values = [y_mat[i][j-1],y_mat[i][j]]
                    
                    #plt.plot(x_values, y_values, color = 'blue')
                    plt.plot(x2_values,y2_values,color = 'black')
                    
                elif i > j:
                    phi_mat[i][j] = 0.5*(-nu_mat[i-1][j] + phi_mat[i-1][j] + nu_mat[i][j-1] + phi_mat[i][j-1])
                    nu_mat[i][j] = 0.5*(nu_mat[i-1][j] - phi_mat[i-1][j] + nu_mat[i][j-1] + phi_mat[i][j-1])
                    M_mat[i][j] = M_from_nu(nu_mat[i][j], gamma)
                    P_mat[i][j] = isentropic_pressure_from_M(Pe, Me, M_mat[i][j], gamma)
                    mu_mat[i][j] = Mach_angle(M_mat[i][j])
                    
                    slope_m_mat[i][j] = np.tan(phi_mat[i][j] - mu_mat[i][j])
                    slope_p_mat[i][j] = np.tan(phi_mat[i][j] + mu_mat[i][j])
                    
                    coordinates_2 = intersection(x_mat[i][j-1],y_mat[i][j-1],slope_m_mat[i][j],x_mat[i-1][j], y_mat[i-1][j], slope_p_mat[i-1][j]) 
                    x_mat[i][j] = coordinates_2[0]
                    y_mat[i][j] = coordinates_2[1]
                    
                    x_values = [x_mat[i][j-1], x_mat[i][j]]
                    x2_values = [x_mat[i-1][j],x_mat[i][j]]
                    y_values = [y_mat[i][j-1], y_mat[i][j]]
                    y2_values = [y_mat[i-1][j],y_mat[i][j]]
                    
                    plt.plot(x_values, y_values, color = 'black')
                    plt.plot(x2_values,y2_values,color = 'black')
                    
                else:
                    continue   
                
                
    nu_end_arr = []
    phi_end_arr = []
    x_end_arr = []
    y_end_arr = []
    M_plot =[]
    x_plot =[]
    y_plot =[]  
    x_start= x_arr
    y_start = y_arr
    M_start =M_arr
    
                  
    for j in range(0,n_chr):
        nu_end_arr = np.append(nu_end_arr, nu_mat[n_chr-1][j])
        phi_end_arr = np.append(phi_end_arr, phi_mat[n_chr-1][j])
        x_end_arr = np.append(x_end_arr, x_mat[n_chr-1][j])
        y_end_arr = np.append(y_end_arr, y_mat[n_chr-1][j])    
        for i in range(0,n_chr):
            if i >= j:
                M_plot = np.append(M_plot, M_mat[i][j])  
                x_plot = np.append(x_plot, x_mat[i][j])
                y_plot = np.append(y_plot, y_mat[i][j])
                
    for i in range (0,n_chr):
        x_start = np.append(x_start, x_mat[i][0])
        y_start = np.append(y_start, y_mat[i][0])
        M_start = np.append(M_start, M_mat[i][0])   
    
    slope_end = np.tan(phi_mat[n_chr-10][0] + nu_mat[n_chr-1][0])    
    x_1region = [x_end,x_arr[0], x_start[n_chr]]  
    y_1region = [y_end,y_arr[0], y_start[n_chr]]  
    M_1region = [M_end,x_arr[0], M_start[n_chr]]  
    
      
    triang = tri.Triangulation(x_plot,y_plot)
    contour = plt.tricontourf(triang, M_plot, levels=50, cmap='jet',vmin = 2, vmax = 3)  # Adjust levels and cmap as needed
    plt.clim(2, 3)
        
    triang_Machinitial = tri.Triangulation(x_start, y_start)
    plt.tricontourf(triang_Machinitial, M_start, levels = 50, cmap = 'jet', vmin = 2, vmax = 3)
    
    triang_constant = tri.Triangulation(x_1region, y_1region)
    plt.tricontourf(triang_constant, M_1region, levels = 50, cmap = 'jet', vmin = 2, vmax = 3)
    
    return nu_end_arr, phi_end_arr ,x_end_arr, y_end_arr,slope_end      
          


#nozzle exit point and initial fan properties
x0 = 0
y0 = 1
Me = 2
phi_e = 0
P_atm = 1
Pe = 1.2*P_atm
gamma = 1.4

M1 = isentropic_mach(Pe,P_atm, Me, gamma)
nu_e = prandtl_meyer_function(Me, gamma)
nu_1 = prandtl_meyer_function(M1, gamma)
mu_e = Mach_angle(Me)
mu_1 = Mach_angle(M1)

n_chr = int(input('Number of characteristics:'))    #input number of characteristics to generate

phi_1 = nu_1 - nu_e
nu_arr = []
M_arr = []

phi_arr = [phi_e]
del_phi = (phi_1 - phi_e)/(n_chr - 1)   #incremental change in characteristics

#list of Phi, mu and M values at the end points
for i in range(1, n_chr):
    phi_arr = np.append(phi_arr, phi_e +(del_phi*i))
for i in range(0, n_chr):
    nu_arr = np.append(nu_arr, nu_e + phi_arr[i])
#for i in range(0, n_chr):
   # M_arr = np.append(M_arr, M_from_nu(nu_arr[i], gamma))

new_arr = gamma_minus_reflections_initial(nu_arr, phi_arr)

nu_arr = new_arr[0]
phi_arr = new_arr[1]
x_arr = new_arr[2]
y_arr = new_arr[3]
M_arr = new_arr[4]

#for triangular region with P < P_atm
x_end = x_arr[n_chr-1]
y_end = y_arr[n_chr-1]
M_end = M_arr[n_chr-1]


new_arr = gamma_plus_reflections(nu_arr, phi_arr, x_arr, y_arr, M_arr)
nu_arr = new_arr[0]
phi_arr = new_arr[1]
x_arr = new_arr[2]
y_arr = new_arr[3]
M_arr = new_arr[4]


new_arr = gamma_minus_reflections(nu_arr, phi_arr, x_arr, y_arr,M_arr,x_end,y_end, M_end)
slope_final = new_arr[4]
X_final = new_arr[2]
Y_final = new_arr[3]

Phi_extension = intersection(x_arr[n_chr-1], y_arr[n_chr-1], np.tan(phi_arr[n_chr-1]), X_final[0], Y_final[0],np.tan(slope_final))

x_extension = Phi_extension[0]
y_extension = Phi_extension[1]

xvalues = [x_extension,x_arr[n_chr-1] ]
yvalues = [y_extension,y_arr[n_chr-1] ]
#Plot all the computed contours and lines

plt.ylim(0,5)

plt.plot(xvalues, yvalues, c = 'red')
plt.xlabel('X-coordinate')
plt.ylabel('Y-coordinate')
plt.title('Mach number Contour for M_exit =2')

plt.show()

