import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


R_E = 6378
u = 398600

def two_body(t, X):

    r = np.sqrt(X[0]**2 + X[1]**2 + X[2]**2)
    X_dot = np.array([X[3], X[4], X[5], -u*X[0]/r**3, -u*X[1]/r**3, -u*X[2]/r**3])
    
    return X_dot

if __name__ == "__main__":
    theta = 100 * np.pi / 180 
    altitude = 404
    r_iss = R_E + altitude
    x_iss = r_iss * np.cos(theta)
    y_iss = r_iss * np.sin(theta)
    z_iss = 0
    vx_iss = -np.sqrt(u/r_iss) * np.sin(theta)
    vy_iss = np.sqrt(u/r_iss) * np.cos(theta)
    vz_iss = 0
    X_iss = np.array([x_iss, y_iss, z_iss, vx_iss, vy_iss, vz_iss])
    a_iss = R_E + altitude
    tp_iss = (2 * np.pi * np.sqrt(a_iss**3 / u) )
    x_center, y_center = 0, 0

    print(f"Period of ISS: {tp_iss/3600:.2f}  hours")

    t_eval = np.arange(0, 15*tp_iss, 10)
    two_body(t_eval, X_iss)
    sol_iss = solve_ivp(two_body, [0, 15*tp_iss], X_iss, t_eval=t_eval, method='RK45', rtol=1e-12, atol=1e-12)
    x_iss = sol_iss.y[0]
    y_iss = sol_iss.y[1]
    z_iss = sol_iss.y[2]

    N_rev = 12
    a_chaser = a_iss * (1 - (theta / (2 * np.pi * N_rev)))**(2/3)
    e_chaser = a_iss / a_chaser - 1
    
    v_chaser = np.sqrt(u*((2/r_iss) - (1/a_chaser)))
    print(f"Semi-major axis of chaser: {a_chaser:.2f} km")
    print(f"Eccentricity of chaser: {e_chaser:.2f}")

    x_chaser = r_iss
    y_chaser = 0
    z_chaser = 0
    vx_chaser = 0
    vy_chaser = v_chaser
    vz_chaser = 0
    tp_chaser = (2 * np.pi * np.sqrt(a_chaser**3 / u) )
    X_chaser = np.array([x_chaser, y_chaser, z_chaser, vx_chaser, vy_chaser, vz_chaser])
    print(X_chaser)
    sol_chaser = solve_ivp(two_body, [0, 15*tp_iss], X_chaser, t_eval=t_eval, method='RK45', rtol=1e-12, atol=1e-12)
    x_chaser, y_chaser, z_chaser = sol_chaser.y[0], sol_chaser.y[1], sol_chaser.y[2]

    distance = np.sqrt((x_iss - x_chaser)**2 + (y_iss - y_chaser)**2 + (z_iss - z_chaser)**2)
 
    plt.plot(x_iss, y_iss, label="ISS Orbit")
    plt.plot(x_chaser, y_chaser, label="Chaser Orbit", linestyle='--')
    plt.scatter(x_iss[0], y_iss[0], color='red', label="Initial Point ISS", zorder=5, marker='o') 
    plt.scatter(x_chaser[0], y_chaser[0], color='red', label="Initial Point Chaser", zorder=5, marker='*')  
    circle = plt.Circle((x_center, y_center), R_E, color='green', fill=True)
    ax = plt.gca()  
    ax.add_patch(circle)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel('x (km)')
    plt.ylabel('y (km)')
    plt.title('ISS and Chaser Orbits')
    plt.legend()
    plt.legend(loc="upper right") 
    plt.grid()
    plt.show() 
 
 
 
    plt.plot(t_eval/tp_chaser, distance)
    plt.xlabel('Number of Orbit')
    plt.ylabel('Distance (km)')
    plt.title('Distance Between ISS and Chaser')
    plt.grid()


    zero_distance_index = np.argmin(distance)  
    plt.scatter(t_eval[zero_distance_index]/tp_chaser, distance[zero_distance_index], color='red', label='Zero Distance', zorder=5, marker='o')

    plt.legend()
    plt.show()

    v_iss = np.sqrt(u/r_iss)
    N_rev_range = np.arange(2, 31)
    delta_v_values = []

    for N_rev in N_rev_range:
        a_chaser = a_iss * (1 - (theta / (2 * np.pi * N_rev)))**(2/3)
        e_chaser = a_iss / a_chaser - 1
        r_apoapsis = a_chaser * (1 + e_chaser)
        v_apoapsis_chaser = np.sqrt(u*((2/r_apoapsis) - (1/a_chaser)))
        delta_v = abs(v_apoapsis_chaser - v_iss)
        delta_v_values.append(delta_v)
    
    plt.plot(N_rev_range, delta_v_values, marker='o', color='blue', label=r'$\Delta V$')
    plt.xlabel('Number of Revolutions (N_rev)')
    plt.ylabel(r'$\Delta V$ (km/s)')
    plt.title(r'$\Delta V$ Required to Match Apogee Speed to ISS')
    plt.grid(True)
    plt.legend()
    plt.show()

    v_circular = np.sqrt(u/r_iss)
    perigee_altitudes = np.arange(60, 211, 10)  # Perigee altitudes from 60 km to 210 km
    delta_v_values = []

    for perigee_altitude in perigee_altitudes:
        r_perigee = R_E + perigee_altitude
        a_new = (r_iss + r_perigee) / 2
        v_apogee_new = np.sqrt(u*((2/r_iss) - (1/a_new)))
        delta_v = abs(v_apogee_new - v_circular)
        delta_v_values.append(delta_v)
    
    plt.plot(perigee_altitudes, delta_v_values, marker='o', color='blue', label=r'$\Delta V$')
    plt.xlabel('Perigee Altitude (km)')
    plt.ylabel(r'$\Delta V$ (km/s)')
    plt.title(r'$\Delta V$ Required to Lower Perigee')
    plt.grid(True)
    plt.legend()
    plt.show()
    
