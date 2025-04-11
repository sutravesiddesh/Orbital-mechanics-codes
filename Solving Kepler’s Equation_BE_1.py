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

def altitude(tp,n,tol,e,a):
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
        time_series.append(t/tp)
        theta_series.append(theta)
    print(no_iter)
    return alt_series,time_series, no_iter, theta_series

def orbit(theta_series,alt_series):
    series_x = []
    series_y = []
    for i in range(len(alt_series)):
        altitude = alt_series[i]
        x = (altitude + R_E) * np.cos(theta_series[i])
        y = (altitude + R_E) * np.sin(theta_series[i])
        series_x.append(x)
        series_y.append(y)
    
    return series_x, series_y

def theta(e,p,a,R_E):
    
    alpha = R_E**2 * e**2 + p**2
    beta = 2 * R_E**2 * e
    gamma = R_E**2 - p**2

    discriminant = beta**2 - 4*alpha*gamma
    if discriminant < 0:
        return None
    
    cos_theta1 = (-beta + np.sqrt(discriminant)) / (2*alpha)
    cos_theta2 = (-beta - np.sqrt(discriminant)) / (2*alpha)

    theta1 = np.arccos(cos_theta1)
    theta2 = np.arccos(cos_theta2)

    angles = [theta1, -theta1, theta2, -theta2]
    markers = ['*', '*', 'o', 'o']
    colors = ['red', 'red', 'red', 'red']

    for angle, marker, color in zip(angles, markers, colors):
        r = radial_distance(angle, a)
        x = (r) * np.cos(angle)
        y = (r) * np.sin(angle)
        plt.scatter(x, y, color=color, marker=marker,s=100, label=f"θ={angle:.2f} rad")
    
    return angles, theta1, theta2

def eclipse_time(theta1,theta2,e,n,tp):
    E1 = 2*np.arctan(np.tan(theta1/2)*np.sqrt((1-e)/(1+e)))
    E2 = 2*np.arctan(np.tan(-theta1/2)*np.sqrt((1-e)/(1+e)))
    E3 = 2*np.arctan(np.tan(theta2/2)*np.sqrt((1-e)/(1+e)))
    E4 = 2*np.arctan(np.tan(-theta2/2)*np.sqrt((1-e)/(1+e)))
    M1 = E1 - e*np.sin(E1)
    M2 = E2 - e*np.sin(E2)
    M3 = E3 - e*np.sin(E3)
    M4 = E4 - e*np.sin(E4)
    periapsis_eclipse_time = (abs(M2 - M1) / n)/60
    apoapsis_eclipse_time = (tp - M3/n + M4/n)/60 
    return periapsis_eclipse_time, apoapsis_eclipse_time

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
    tp= (2*np.pi)/n
    p = a*(1-e**2)
    x_center, y_center = 0, 0
    
    alt_series, time_series,no_iter, theta_series = altitude(tp,n,tol,e,a)

    series_x,series_y = orbit(theta_series,alt_series)

    plt.plot(time_series,alt_series)
    plt.xlabel("Time (% of orbital period)")
    plt.ylabel("Altitude (km)")
    plt.show()
    
    plt.scatter(series_x,series_y)
    circle = plt.Circle((x_center, y_center), R_E, color='green', fill=True)
    ax = plt.gca()  # Get the current axes
    ax.add_patch(circle)
    ax.set_xlim(-R_E-40000, R_E+5000)
    ax.set_ylim(-R_E-15000, R_E+15000)
    plt.gca().set_aspect('equal', adjustable='box')
    angles, theta1, theta2 = theta(e, p, a, R_E)
    plt.show()
  
    eclipse_duration = eclipse_time(theta1, theta2, e, n,tp)
    print(f"Periapsis eclipse duration: {eclipse_duration[0]} minutes, Apoapsis eclipse duration: {eclipse_duration[1]} minutes")
    print(f"Periapsis angle: {theta1}, Apoapsis angle: {theta2}")
