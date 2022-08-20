# DELTA PARALLEL ROBOT

![delta robot](https://www.linearmotiontips.com/wp-content/uploads/2016/01/Delta-Robot-Diagram.jpg)

this will be a step by step demonstration of how to control trajectory of a Delta robot End-Effector
check out my [telegram channel](https://t.me/engineering_stuff_69)

# 1 - ROBOTICS (INVERSE AND FORWARD KINEAMTICS)
with trajectory planning there are two fundemental questions we need an answer for:
1. given a specific location in the real world, what values should my robot's joint be set to in order to get the End-Effector there? (inverse kinematics)
2. given the setting of my joints, where is my EE in real world coordinates? (forward kinematics)

### Theory
![delta robot kinemtaics figure 1](https://i.ibb.co/cc29GYf/Delta-robot-kin.png)
![delta robot kinemtaics figure 2](https://i.ibb.co/VVVQfkF/Delta-robot-kin-2.png)

The Delta robot is a 3-DOF robot that consists of two parallel platforms. One of them, which we call the EE platform, is capable of moving and the other one, which we call base platform, is not capable of that. These two platforms are connected with three arms. Each arm has one pin and two universal joints (as shown in the figure) that connect two solid rods. The rod which is connected to base platform by the pin, is called the active rod, and the other one is called the passive rod. The center of the base and EE platforms are marked as $O$ and $O'$ respectively.  The pin is denoted as A, the universal joint connecting the passive and active rods, is denoted as B. The universal joint connecting the EE platform and the passive rod is denoted as C. 

$\overrightarrow{(O O')} + \overrightarrow{(O' C_i)} = \overrightarrow{(O A_i)} + \overrightarrow{(A_i B_i)} + \overrightarrow{(B_i C_i)}$ (Equation 1)

- IK: so you'll need to solve this equation and find $\theta_{ij}$ with respect to the other variables. this will solve the inverse kinematics problem. $\theta_{1j}$ are the angles of the actuator joints.
- FK: for forward kinematics it, the equation number 1 is solved for $p_x, p_y, p_z$ (position of the EE) with respect to other variables (which is a much easier problem that IK) 

### References 
- here's a [good playlist](https://www.youtube.com/playlist?list=PLjx2FAhpTe3FGbcjBbxlhf56qVR0XbVNO) for learning FK and IK

- [Delta Robot Inverse Kinematics](https://sites.google.com/site/deltarobotberkeley/how-it-works)

- [Delta Robot Inverse Kinematics method 1](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/theory/Inverse%20Kinematics%20(Delta%20Robot).pdf)

- [Delta Robot Inverse Kinematics method 2](https://www.researchgate.net/publication/242073945_Delta_robot_Inverse_direct_and_intermediate_Jacobians)

### Python implementation
you can see my [python implementation of IK](https://github.com/ArthasMenethil-A/Delta_Robot/tree/main/inverse%20and%20forward%20kinematics):
- [method 1](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/inverse%20and%20forward%20kinematics/IK_method_1.py)
- [method 2](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/inverse%20and%20forward%20kinematics/IK_method_2.py)

# 2 - POINT TO POINT MOVEMENT
this section is dedicated to answer how should you go about writing a code for point to point movement (moving the EE from point 1 to point 2 in 3d space )

### theory
3-4-5 polynomial and 4-5-6-7 polynomial point to point movement, you can learn about this in the book ["Fundamentals of Robotic 
Mechanical Systems, theory, methods, and Algorithms, Fourth Edition by Jorge Angeles - chapter 6"](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/theory/Angles_a3hfp_Fundamentals_of_Robotic.pdf)

### python implementation
- POINT TO POINT MOVEMENT (3-4-5 polynomial):
  moveing from point to point in the direction desired --> [trajectory_planning_345.py](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/point%20to%20point%20movement%20(python)/trajectory_planning_345.py)

![345 point-to-point](https://i.ibb.co/HGQvTwT/345.png)
- POINT TO POINT MOVEMENT (4-5-6-7 polynomial):
  we repeat what we've done for sub-step 2 but with a 7th order polynomial --> [trajectory_planning_4567.py](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/point%20to%20point%20movement%20(python)/trajectory_planning_4567.py)

# 3 - TRAJECTORY PLANNING (CUBIC SPLINE AND SIMILAR ALGORITHMS)
this section is dedicated to planning out a specific trajectory for the robot to go through

### theory
cubic spline (book Trajectory Planning for Automatix Machines and Robots by Luigi Biagiotti and Claudio Melchiorri)
1. cubic spline with assigned initial and final velocities (part 4.4.1)
2. cubic spline with assigned intial and final velocities and acceleration (part 4.4.4)
3. smoothing cubic spline (part 4.4.5)

### python implementation
- CIRCLE MOVEMENT:
  cubic spline with assigned initial and final velocities --> [trajectory_planing_cubic_spline.py](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/trajectory%20planning%20-%20cubic%20splin%20(python)/trajectory_planning_cubic_spline.py)

- CIRCLE MOVEMENT:
  cubic spline with assigned initial and final velocities and acceleration --> [trajectory_planing_cubic_spline_4.4.4.py](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/trajectory%20planning%20-%20cubic%20splin%20(python)/trajectory_planning_cubic_spline_4.4.4.py)

# 2.5 - JACOBIAN
The jacobian matrix relates velocity of EE to the velocity of actuator joints with the relation: $\vec{v} = J \dot{\vec{\theta}}$

### theory
for the IK we have the equations as followed 
- $p_x \cos(\Phi_i) - p_y \sin(\Phi_i) = R - r + a \cos(\theta_{1i}) + b \sin(\theta_{3i}) \cos(\theta_{2i} + \theta_{1i})$       (equation 1)
- $p_x \sin(\Phi_i) + p_y \cos(\Phi_i) = b \cos(\theta_{3i})$       (equation 2)
- $p_z = a \sin(\theta_{2i}) + b \sin(\theta_{3i}) \sin(\theta_{2i} + \theta_{1i})$       (equation 3)

from these equation we solve for $\theta_{ij}$ (solution of IK) and solve for $p_x, p_y, p_z$ (solution of FK)

1. [what is Jacobian matrix and how it is calculated](https://www.sciencedirect.com/science/article/pii/S1877050918309876)
2. [how jacobian matrix is used in delta robot trajectory palanning](http://jai.front-sci.com/index.php/jai/article/view/505)
3. [detailed calculation and analysis of Jacobian matrix in delta robot](http://jai.front-sci.com/index.php/jai/article/view/505)

### python implementation

- [jacobian file: IK, FK, jacobian matrix calculation](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/inverse%20and%20forward%20kinematics/IK_method_2.py)

# 3 - TRAPEZOIDAL TRAJECTORY PLANNING
  
