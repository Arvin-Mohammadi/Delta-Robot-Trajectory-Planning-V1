# Delta_Robot

![delta robot](https://www.linearmotiontips.com/wp-content/uploads/2016/01/Delta-Robot-Diagram.jpg)

this will be a step by step demonstration of how to control trajectory of a Delta robot End-Effector
check out my [telegram channel](https://t.me/engineering_stuff_69)

# 1. ROBOTICS (INVERSE AND FORWARD KINEAMTICS)
with trajectory planning there are two fundemental questions we need an answer for:
1. given a specific location in the real world, what values should my robot's joint be set to in order to get the End-Effector there? (inverse kinematics)
2. given the setting of my joints, where is my EE in real world coordinates? (forward kinematics)

### theory
![delta robot kinemtaics figure 1](https://i.ibb.co/cc29GYf/Delta-robot-kin.png)
![delta robot kinemtaics figure 2](https://i.ibb.co/VVVQfkF/Delta-robot-kin-2.png)

The Delta robot is a 3-DOF robot that consists of two parallel platforms. One of them, which we call the EE platform, is capable of moving and the other one, which we call base platform, is not capable of that. These two platforms are connected with three arms. Each arm has one pin and two universal joints (as shown in the figure) that connect two solid rods. The rod which is connected to base platform by the pin, is called the active rod, and the other one is called the passive rod. The center of the base and EE platforms are marked as $O$ and $O'$ respectively.  The pin is denoted as A, the universal joint connecting the passive and active rods, is denoted as B. The universal joint connecting the EE platform and the passive rod is denoted as C. 

$\overrightarrow{(O O')} + \overrightarrow{(O' C_i)} = \overrightarrow{(O A_i)} + \overrightarrow{(A_i B_i)} + \overrightarrow{(B_i C_i)}$

so you'll need to solve this equation and find $\theta_ij$ with respect to the other variables. this will solve the inverse kinematics problem. $\theta_1j$ are the angles of the actuator joints.

here's a [good playlist](https://www.youtube.com/playlist?list=PLjx2FAhpTe3FGbcjBbxlhf56qVR0XbVNO) for learning FK and IK <br />
<br />
[here's how you would compute Delta Robot IK](https://sites.google.com/site/deltarobotberkeley/how-it-works) <br />
<br />
[also this pdf explains the theory and implementation pretty easily](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/theory/Inverse%20Kinematics%20(Delta%20Robot).pdf) <br />
<br />

### python implementation
you can see my python implementation of IK in the file [trajectory_planning_345.py](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/point%20to%20point%20movement%20(python)/trajectory_planning_345.py)

# 2. POINT TO POINT MOVEMENT
this section is dedicated to answer how should you go about writing a code for point to point movement (moving the EE from point 1 to point 2 in 3d space )

### theory
3-4-5 polynomial and 4-5-6-7 polynomial point to point movement, you can learn about this in the book ["Fundamentals of Robotic 
Mechanical Systems, theory, methods, and Algorithms, Fourth Edition by Jorge Angeles - chapter 6"](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/theory/Angles_a3hfp_Fundamentals_of_Robotic.pdf)

### python implementation
- POINT TO POINT MOVEMENT (3-4-5 polynomial):
  moveing from point to point in the direction desired --> [trajectory_planning_345.py](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/point%20to%20point%20movement%20(python)/trajectory_planning_345.py)

- POINT TO POINT MOVEMENT (4-5-6-7 polynomial):
  we repeat what we've done for sub-step 2 but with a 7th order polynomial --> [trajectory_planning_4567.py](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/point%20to%20point%20movement%20(python)/trajectory_planning_4567.py)

# 3. TRAJECTORY PLANNING (CUBIC SPLINE AND SIMILAR ALGORITHMS)
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

# 3. TRAJECTORY PLANNING (TRYING OTHER ALGORITHMS) 

## 1) Jacobian 

### theory
1. [what is Jacobian matrix and how it is calculated](https://www.sciencedirect.com/science/article/pii/S1877050918309876)
2. [how jacobian matrix is used in delta robot trajectory palanning](http://jai.front-sci.com/index.php/jai/article/view/505)
3. [detailed calculation and analysis of Jacobian matrix in delta robot](http://jai.front-sci.com/index.php/jai/article/view/505)

### python implementation

