
#Newton-Raphson method

import numpy as np
import matplotlib.pyplot as plt

R_E = 6378
u = 398600

#Newton-Raphson method 
def Kepler(M,e,tol):
    
    E = M
    delta_E = 1
    no_iter =  0
    while abs(delta_E) > tol:

        f = E - e*np.sin(E) - M
        delta_f = 1 - e*np.cos(E)
        delta_E = -(f/delta_f)
        E = E + delta_E
        no_iter += 1
    
    return E,no_iter

def true_anomaly(E):
    
    theta = 2*np.arctan(np.tan(E/2)*np.sqrt((1+e)/(1-e)))
    
    return theta
    
def radial_distance(theta,a):
    
    r = a * (1-e**2) / (1+e*np.cos(theta))

    return r

def altitude(p,n,tol,e,a):
    alt_series = []
    time_series = []
    theta_series = []
    for t in range(0,37005,15):
        M = n*t
        E, no_iter = Kepler(M,e,tol)
        theta = true_anomaly(E)
        r = radial_distance(theta,a)
        alt = r - R_E
        alt_series.append(alt)
        time_series.append(t/p)
        theta_series.append(theta)
    print(no_iter)
    return alt_series,time_series, no_iter, theta_series

def orbit(theta_series):
    series_x = []
    series_y = []
    for i in range(len(alt_series)):
        altitude = alt_series[i]
        x = (altitude + R_E) * np.cos(theta_series[i])
        y = (altitude + R_E) * np.sin(theta_series[i])
        series_x.append(x)
        series_y.append(y)
    
    return series_x, series_y


if __name__ == "__main__":
    # M = 21*(np.pi/180)
    # #M = 1.792 
    # e = 0.25
    # tol = 1e-12
    # a = 24000
    # E, no_iter = Kepler(M,e,tol)
    # theta = true_anomaly(E)
    # r = radial_distance(theta,a)
    # print(E, no_iter,theta,r)

    u = 398600
    a = 24000
    e = 0.72
    tol = 1e-12
    n= np.sqrt(u/a**3)
    p= (2*np.pi)/n
    
    alt_series, time_series,no_iter, theta_series = altitude(p,n,tol,e,a)

    series_x,series_y = orbit(theta_series)

    plt.plot(time_series,alt_series)
    plt.xlabel("Time (% of orbital period)")
    plt.ylabel("Altitude (km)")
    plt.show()
    
    plt.scatter(series_x,series_y)
    plt.show

    

