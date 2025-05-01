# AE508 Course Project: Intercept and Deflect an Asteroid
Design an optimal trajectory for a spacecraft on a geostationary orbit
to intercept an asteroid that is on a collision course with the Earth. 


## Team

* petern4@illinois.edu
* pc46@illinois.edu
* davisr2@illinois.edu

## Contents

* `final_report.m`  - plots the converged solutions. The converged solutions are in `reports/  tradeoff_studies_results_final.txt`.

* `deflection_report.m` - iterates over the converged solutions, and for each solution computes a change in the velocity of the asteroid caused by the impactor using the following equation (Wang et al):
```
dv(tf) = beta*(m(tf)/m_ast(tf))*(v(tf) - v_ast(tf))
```
Where `dv(tf),v(tf),v_ast(tf)` are vectors corresponding to the change in the asteroid velocity, the velocity of the impactor spacecraft at the time of impact and the velocity of the asteroid before the collision. The velocities are in ECI, and the masses are in kg. `beta=1` is a scaling parameter, and can be used to vary the amount of material ejected due to the collision. 
The script uses the estimated change in the asteroid velocity to propagate the asteroid trajectory (two-body EOMs, no perturbing forces, ECI reference frame, the Earth is the primary central body, so this may not be correct) to find a set of solutions that would push the asteroid further away from the Earth. Finally, the script plots four best solutions, the estimated change in the velocity of the asteroid, the time of flight and the asteroid mass.

* `eom.m` - computes the state and costate dynamics.

* `hamiltonian.m` - computes the Hamiltonian.

* `initial_values.m` - sets up the the location of the impactor on a departure orbit, the location of the asteroid, and other parameters.

* `lambert.m` - the Lambert's problem solver (not used).

* `max_momentum.m` - evaluates the cost function.

* `plots.m` - a set of utility functions to plot the trajectory, states, the control and costates.

* `tradeoff_studies.m` - solves the asteroid intercept optimization problem by sweeping over thrust, longitude (true anomaly) and time of flight. Converged solutions are stored in `reports/tradeoff_studies_results.txt`.

* `two_body.m` - two-body EOMs, no perturbing forces.

* `classical2posvel.m` - converts the classical orbital elements to the position/velocity vectors.

* `Asteroid_Intercept_AE508_Spring2025.m` - solves the optimal trajectory problem, given the time of flight is 17 days, the engine's thrust magnitude is 0.5 Newton, and the longitude (the true anomaly) on the departure orbit is 340 deg. 

`reports` - contains the plots to be included into the final report.

## Initial conditions
- The impactor satellite is in a location correspoding to a GEO orbit. It is also assumed that the impactor has already done a burn to put in a parabolic trajectory.
- The satellite's longitude on the departure orbit varies
- The reference frame is Earth' Inertial (J2000)
- The position of the hypothetical asteroid is set to JD=2462239.715277778, UT = A.D. 2029-Apr-13 05:10:00.0000
```
   position_asteroid = [-2.956932341462374E+05, -2.043495023251738E+05, -9.847039118095844E+04] km
   velocity_asteroid = [4.292596028662857E+00,  3.878022383303825E+00,  1.677451136476388E+00] km/s 
```
- The mass of the hypothetical asteroid used to determine the change in the asteroid velocity after colliding with the impactor spacecraft is set to be a range between 2.2e6 to 2.2e8 kg, where 2.2e8 is the mass of 2024 YR4 asteroid
- The impactor's wet mass is 500 kg.
- The impactor's engine thrust magnitude and the specific impulse are varied, as shown:

| Thrust Magnitude, N | Isp (Specific Impulse), sec|
|---------------------|---------------------------|
| 0.05 | 4190.0 |
| 0.1 | 4000.0 |
| 0.25 | 3000.0 |
| 0.5 | 1000.0 |
| 1.0 | 800.0 |
| 1.5 | 650.0 |
| 2 |  400.0 |
| 2.5 |  250.0 |

## Methodology
See the notes for the stationary conditions, Euler-Lagrange equations and the final boundary conditions.

## How to run 
To find an optimal trajectory for the impactor satellite departing on a parabolic orbit to intercept a hypothetical asteroid, open `Asteroid_Intercept_AE508_Spring2025.m` and click `Run`. The time of flight is fixed and is set to 10 days, the maximum thrust magnitude is 1 N.

## How to run tradeoff studies
In MATLAB, open `tradeoff_studies.m` and click `Run`. The script will sweep over the time of flight from 4 to 10 days, thrust magnitude between 0.05 and 2 N and the true anomaly angle of the impactor satellite on the departure orbit from 0 to 360 deg with 10 deg step. The script will generate the results and store them in `reports/tradeoff_studies_results.txt`

## Remarks
Getting the solution to converge has been a challenge, especially when using low-thrust engine with the thrust magnitude less than 1 N and the time of flight > 10-15 days. The MATLAB environment frequently reported memory allocation errors when the time of flight exceeded 10 or 15 days. To mitigate the large memory requirement issues, use `linspace(t0,tf,1000)` to explicitly specify how the time interval is discretized.

## How the results compare to DART Mission
According to the DART mission report (https://dart.jhuapl.edu/Mission/Impactor-Spacecraft.php#:~:text=DART%20navigated%20to%20crash%20itself,(580%20kilograms)%20at%20impact) the DART crashed into Dimorphos at a speed of approximately 6.1 km/s, with its total mass 580 kg at the time 
of impact (the mission used hydrazine propellant (50 kg) for spacecraft maneuvers and altitude control, and xenon (60 kg) to operate the Ion propulsion).



Source: https://dart.jhuapl.edu/Mission/Impactor-Spacecraft.php#:~:text=DART%20navigated%20to%20crash%20itself,(580%20kilograms)%20at%20impact.

