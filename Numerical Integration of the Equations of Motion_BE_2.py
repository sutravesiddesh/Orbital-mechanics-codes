import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


R_E = 6378
u = 398600

# Two body equations of motion in 3D

def two_body(t, X):

    r = np.sqrt(X[0]**2 + X[1]**2 + X[2]**2)
    # v = r_dot = np.array([vx, vy, vz])
    X_dot = np.array([X[3], X[4], X[5], -u*X[0]/r**3, -u*X[1]/r**3, -u*X[2]/r**3])
    
    return X_dot
    
# def specific_energy(X):
#     r = np.sqrt(X[0]**2 + X[1]**2 + X[2]**2)
#     v = np.sqrt(X[3]**2 + X[4]**2 + X[5]**2)
#     K_E = (v**2)/2
#     U_E = -u/r
#     Total_E = K_E + U_E
#     return Total_E,K_E,U_E

if __name__ == "__main__":
    X0 = np.array([7115.804, 3391.696,3492.221,-3.762, 4.063, 4.184])
    t_val = np.arange(0, 86400, 10)
    two_body(t_val, X0)
    sol = solve_ivp(two_body, [0, 86400], X0, t_eval=t_val, method='RK45', rtol=1e-12, atol=1e-12)
    time = sol.t / 3600

    # Plotting Magnitude of position vs time
    radius = np.sqrt(sol.y[0]**2 + sol.y[1]**2 + sol.y[2]**2)
    plt.plot(time, radius, label="Radius (km)", color="blue")
    plt.xlabel("Time (hours)")
    plt.ylabel("Radius (km)")
    plt.title("Magnitude of Position Vector vs Time")
    plt.grid(True)
    plt.legend()
    plt.show()

    # Plotting Magnitude of velocity vector vs time
    velocity = np.sqrt(sol.y[3]**2 + sol.y[4]**2 + sol.y[5]**2)
    plt.plot(time, velocity, label="Velocity (km/s)", color="red")
    plt.xlabel("Time (hours)")
    plt.ylabel("Velocity (km/s)")
    plt.title("Magnitude of Velocity Vector vs Time")
    plt.grid(True)
    plt.legend()    
    plt.show()

    # Plotting orbit in 3D
    x = sol.y[0]
    y = sol.y[1]
    z = sol.y[2]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot3D(x, y, z, label="Orbit", color="green")
    # Earth sphere for visualization
    phi = np.linspace(0, np.pi, 100)
    theta = np.linspace(0, 2 * np.pi, 100)
    phi, theta = np.meshgrid(phi, theta)
    x_Earth = R_E * np.sin(phi) * np.cos(theta)
    y_Earth = R_E * np.sin(phi) * np.sin(theta)
    z_Earth = R_E * np.cos(phi)
    ax.plot_surface(x_Earth, y_Earth, z_Earth, color='blue', alpha=0.3, label="Earth")
    ax.set_xlabel("X (km)")
    ax.set_ylabel("Y (km)")
    ax.set_zlabel("Z (km)")
    ax.set_title("3D Orbit with Earth")
    ax.legend()
    plt.show()

    # Plotting specific energies vs time
    vx = sol.y[3]
    vy = sol.y[4]
    vz = sol.y[5]
    K_E = (velocity**2)/2
    U_E = -u/radius
    Total_E = K_E + U_E
    plt.plot(time, K_E, label="Kinetic Energy", color="blue")
    plt.plot(time, U_E, label="Potential Energy", color="red")
    plt.plot(time, Total_E, label="Total Energy", color="green")
    plt.xlabel("Time (hours)")
    plt.ylabel("Specific Energy (km²/s²)")
    plt.title("Specific Energies vs. Time")
    plt.grid(True)
    plt.legend()
    plt.show()

    # Specific angular momentum
    hx = y * vz - z * vy
    hy = z * vx - x * vz
    hz = x * vy - y * vx
    h = np.sqrt(hx**2 + hy**2 + hz**2)
    # Plotting specific angular momentum vs time
    plt.plot(time, h, label="Specific Angular Momentum", color="purple",linestyle="-")
    plt.xlabel("Time (hours)")
    plt.ylabel("Specific Angular Momentum (km²/s)")
    plt.title("Specific Angular Momentum vs. Time")
    plt.grid(True)
    plt.legend()
    plt.show() 

    X1 = np.array([0,0,8550,0,-7,0])
    sol1 = solve_ivp(two_body, [0, 86400], X1, t_eval=t_val, method='RK45', rtol=1e-12, atol=1e-12)
    time1 = sol1.t / 3600
    x1, y1, z1 = sol1.y[0], sol1.y[1], sol1.y[2]
    vx1, vy1, vz1 = sol1.y[3], sol1.y[4], sol1.y[5]
    radius1 = np.sqrt(x1**2 + y1**2 + z1**2)
    velocity1 = np.sqrt(vx1**2 + vy1**2 + vz1**2)
    K_E1 = (velocity1**2)/2
    U_E1 = -u/radius1
    Total_E1 = K_E1 + U_E1
    hx1 = y1 * vz1 - z1 * vy1
    hy1 = z1 * vx1 - x1 * vz1
    hz1 = x1 * vy1 - y1 * vx1
    h1 = np.sqrt(hx1**2 + hy1**2 + hz1**2)
    i = np.arccos(hz1.mean()/h1.mean())*180/np.pi
    a = u / (2 * abs(Total_E1.mean()))
    time_period = 2 * np.pi * np.sqrt(a**3 / u)
    print("Semi-major axis (a):", a)
    print("Inclination (i):", i)
    print("Time period (T):", time_period/60, "minutes")

    # Plot radial distance vs. time
    plt.plot(time1, radius1, label="Radial Distance (km)", color="blue")
    plt.xlabel("Time (hours)")
    plt.ylabel("Radial Distance (km)")
    plt.title("Radial Distance vs. Time")
    plt.grid(True)
    plt.legend()
    plt.show()

    # Plot velocity magnitude vs. time
    plt.plot(time1, velocity1, label="Velocity (km/s)", color="red")
    plt.xlabel("Time (hours)")
    plt.ylabel("Velocity (km/s)")
    plt.title("Velocity Magnitude vs. Time")
    plt.grid(True)
    plt.legend()
    plt.show()
