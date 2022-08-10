# Delta_Robot

![delta robot](https://www.linearmotiontips.com/wp-content/uploads/2016/01/Delta-Robot-Diagram.jpg)

this will be a step by step demonstration of how to control trajectory of a Delta robot End-Effector
check out my [telegram channel](https://t.me/engineering_stuff_69)

# 1. ROBOTICS (INVERSE AND FORWARD KINEAMTICS)
with trajectory planning there are two fundemental questions we need an answer for:
1. given a specific location in the real world, what values should my robot's joint be set to in order to get the End-Effector there? (inverse kinematics)
2. given the setting of my joints, where is my EE in real world coordinates? (forward kinematics)

### theory
   here's a [good playlist](https://www.youtube.com/playlist?list=PLjx2FAhpTe3FGbcjBbxlhf56qVR0XbVNO) for learning FK and IK <br />
   <br />
   [here's how you would compute Delta Robot IK](https://sites.google.com/site/deltarobotberkeley/how-it-works) <br />
   <br />
   [also this pdf explains the theory and implementation pretty easily](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/theory/Inverse%20Kinematics%20(Delta%20Robot).pdf) <br />
   <br />

### python implementation
you can see my python implementation of IK in the file [trajectory_planning_345.py](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/python%20implementations/trajectory_planning_345.py)

# 2. POINT TO POINT MOVEMENT
this section is dedicated to answer how should you go about writing a code for point to point movement (moving the EE from point 1 to point 2 in 3d space )

### theory
3-4-5 polynomial and 4-5-6-7 polynomial point to point movement, you can learn about this in the book ["Fundamentals of Robotic 
Mechanical Systems, theory, methods, and Algorithms, Fourth Edition by Jorge Angeles - chapter 6"](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/theory/Angles_a3hfp_Fundamentals_of_Robotic.pdf)

### python implementation
- POINT TO POINT MOVEMENT (3-4-5 polynomial):
  moveing from point to point in the direction desired --> [trajectory_planning_345.py](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/python%20implementations/trajectory_planning_345.py)

- POINT TO POINT MOVEMENT (4-5-6-7 polynomial):
  we repeat what we've done for sub-step 2 but with a 7th order polynomial --> [trajectory_planning_4567.py](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/python%20implementations/trajectory_planning_4567.py)

# 3. TRAJECTORY PLANNING
this section is dedicated to planning out a specific trajectory for the robot to go through

### theory
cubic spline (book Trajectory Planning for Automatix Machines and Robots by Luigi Biagiotti and Claudio Melchiorri)
1. cubic spline with assigned initial and final velocities (part 4.4.1)
2. cubic spline with assigned intial and final velocities and acceleration (part 4.4.4)
3. smoothing cubic spline (part 4.4.5)

### python implementation
- CIRCLE MOVEMENT:
  cubic spline with assigned initial and final velocities --> [trajectory_planing_cubic_spline.py](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/python%20implementations/trajectory_planning_cubic_spline.py)

- CIRCLE MOVEMENT:
  cubic spline with assigned initial and final velocities and acceleration --> [trajectory_planing_cubic_spline_4.4.4.py](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/python%20implementations/trajectory_planning_cubic_spline_4.4.4.py)

# 3. TRAJECTORY PLANNING (TRYING OTHER ALGORITHMS) 

## 1) Jacobian 
### theory
i got the theory and algorithm from [this article](http://jai.front-sci.com/index.php/jai/article/view/505)

### python implementation

