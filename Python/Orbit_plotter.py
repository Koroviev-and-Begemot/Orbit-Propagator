from utils import*

# Constants and initialization
mu_earth = 398600
rad_earth = 6371
mu_moon = 4.904E3
rad_moon = 1737

delta_t = 10
sim_time = 31 * 24 * 60 ** 2
t = 0
k = int(np.ceil(sim_time / delta_t)) + 1

moon_state = np.zeros((6, k))
moon_state_dot = np.zeros((6, k))
moon_state_bar = np.zeros((6, k))
moon_state_dot_bar = np.zeros((6, k))

sat_state = np.zeros((6, k))
sat_state_dot = np.zeros((6, k))
sat_state_bar = np.zeros((6, k))
sat_state_dot_bar = np.zeros((6, k))
a_b = np.zeros((3, k))

m_sat = 1000
T_sat = 1

t_b = [4 * 60 ** 2, 10 * 60 ** 2]
tau_b = [500, 5000]
dir = [1, 2]

t_b.append(sim_time)
tau_b.append(0)
dir.append(1)

moon_state[:, 0] = peri2ECI(384314.45, 0.026, 18.28, 0, 0, 0, mu_earth)
moon_state_dot[:, 0] = np.concatenate((moon_state[3:6, 0], -mu_earth * moon_state[0:3, 0] / np.linalg.norm(moon_state[0:3, 0]) ** 3))

sat_state[:, 0] = peri2ECI(35000, 0.866, 80, 270, 0, 0, mu_moon) + moon_state[:, 0]
sat_state_dot[:, 0] = np.concatenate((
    sat_state[3:6, 0],
    -mu_earth * sat_state[0:3, 0] / np.linalg.norm(sat_state[0:3, 0]) ** 3 -
    mu_moon * (sat_state[0:3, 0] - moon_state[0:3, 0]) / np.linalg.norm(sat_state[0:3, 0] - moon_state[0:3, 0]) ** 3
))

# Calculation loop
i = 1
p = 1

for j in range(len(t_b)):
    while t < t_b[j]:
        moon_state_bar[:, i] = moon_state[:, i - 1] + delta_t * moon_state_dot[:, i - 1]
        moon_state_dot_bar[:, i] = np.concatenate((moon_state_bar[3:6, i], -mu_earth * moon_state_bar[0:3, i] / np.linalg.norm(moon_state_bar[0:3, i]) ** 3))
        moon_state[:, i] = moon_state[:, i - 1] + delta_t / 2 * (moon_state_dot[:, i - 1] + moon_state_dot_bar[:, i])
        moon_state_dot[:, i] = np.concatenate((moon_state[3:6, i], -mu_earth * moon_state[0:3, i] / np.linalg.norm(moon_state[0:3, i]) ** 3))

        sat_state_bar[:, i] = sat_state[:, i - 1] + delta_t * sat_state_dot[:, i - 1]
        sat_state_dot_bar[:, i] = np.concatenate((
            sat_state_bar[3:6, i],
            -mu_earth * sat_state_bar[0:3, i] / np.linalg.norm(sat_state_bar[0:3, i]) ** 3 -
            mu_moon * (sat_state_bar[0:3, i] - moon_state_bar[0:3, i]) / np.linalg.norm(sat_state_bar[0:3, i] - moon_state_bar[0:3, i]) ** 3
        ))
        sat_state[:, i] = sat_state[:, i - 1] + delta_t / 2 * (sat_state_dot[:, i - 1] + sat_state_dot_bar[:, i])
        sat_state_dot[:, i] = np.concatenate((
            sat_state[3:6, i],
            -mu_earth * sat_state[0:3, i] / np.linalg.norm(sat_state[0:3, i]) ** 3 -
            mu_moon * (sat_state[0:3, i] - moon_state[0:3, i]) / np.linalg.norm(sat_state[0:3, i] - moon_state[0:3, i]) ** 3
        ))

        t += delta_t
        i += 1

    while t < t_b[j] + tau_b[j]:
        if tau_b[j] > 0:
            a_b[:, i] = T_sat / m_sat / 1000 * sat_state[3:6, i - 1] / np.linalg.norm(sat_state[3:6, i - 1])
            if dir[j] == 2:
                a_b[:, i] = -a_b[:, i]
        else:
            a_b[:, i] = 0

        moon_state_bar[:, i] = moon_state[:, i - 1] + delta_t * moon_state_dot[:, i - 1]
        moon_state_dot_bar[:, i] = np.concatenate((moon_state_bar[3:6, i], -mu_earth * moon_state_bar[0:3, i] / np.linalg.norm(moon_state_bar[0:3, i]) ** 3))
        moon_state[:, i] = moon_state[:, i - 1] + delta_t / 2 * (moon_state_dot[:, i - 1] + moon_state_dot_bar[:, i])
        moon_state_dot[:, i] = np.concatenate((moon_state[3:6, i], -mu_earth * moon_state[0:3, i] / np.linalg.norm(moon_state[0:3, i]) ** 3))

        sat_state_bar[:, i] = sat_state[:, i - 1] + delta_t * sat_state_dot[:, i - 1]
        sat_state_dot_bar[:, i] = np.concatenate((
            sat_state_bar[3:6, i],
            -mu_earth * sat_state_bar[0:3, i] / np.linalg.norm(sat_state_bar[0:3, i]) ** 3 -
            mu_moon * (sat_state_bar[0:3, i] - moon_state_bar[0:3, i]) / np.linalg.norm(sat_state_bar[0:3, i] - moon_state_bar[0:3, i]) ** 3 +
            a_b[:, i]
        ))
        sat_state[:, i] = sat_state[:, i - 1] + delta_t / 2 * (sat_state_dot[:, i - 1] + sat_state_dot_bar[:, i])
        sat_state_dot[:, i] = np.concatenate((
            sat_state[3:6, i],
            -mu_earth * sat_state[0:3, i] / np.linalg.norm(sat_state[0:3, i]) ** 3 -
            mu_moon * (sat_state[0:3, i] - moon_state[0:3, i]) / np.linalg.norm(sat_state[0:3, i] - moon_state[0:3, i]) ** 3 +
            a_b[:, i]
        ))

        t += delta_t
        i += 1
    if j == len(t_b)-1:
        print('Calculation Complete')
    else:
        print(f'burn {j+1} occurred at {t}s and lasted {tau_b[j]}')


# Plots

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111, projection='3d')
ax1.grid(True)
q = 4E5
ax1.set_xlim([-q, q])
ax1.set_ylim([-q, q])
ax1.set_zlim([-q, q])

drawsphere(0, 0, 0, rad_earth, ax1)

ax1.plot3D(moon_state[0, :], moon_state[1, :], moon_state[2, :], linewidth=2, color='r')
ax1.plot3D(sat_state[0, :], sat_state[1, :], sat_state[2, :], linewidth=2, color='k')

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111, projection='3d')
ax2.grid(True)
drawsphere(0, 0, 0, rad_moon, ax2)
q = 10E3
ax2.set_xlim([-q, q])
ax2.set_ylim([-q, q])
ax2.set_zlim([-q, q])

ax2.plot3D(sat_state[0, :] - moon_state[0, :], sat_state[1, :] - moon_state[1, :], sat_state[2, :] - moon_state[2, :], linewidth=2, color='r')

for i in range(len(t_b) - 1):
    xx = range(int(t_b[i] / delta_t), int((t_b[i] + tau_b[i]) / delta_t))
    ax2.plot3D(sat_state[0, xx] - moon_state[0, xx], sat_state[1, xx] - moon_state[1, xx], sat_state[2, xx] - moon_state[2, xx], linewidth=2, color='b')

plt.show()



