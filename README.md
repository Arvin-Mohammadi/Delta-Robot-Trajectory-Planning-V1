# Delta_Robot

this will be a step by step demonstration of how to control trajectory of a Delta robot End-Effector
check out my [telegram channel](https://t.me/engineering_stuff_69)

### My Great Heading {#1} 









# step 1: THEORY

   1. first you'll need to know a bit of robotics(meaning what forward and inverse kinematics are)
      for that you can use this youtube play list:
      https://www.youtube.com/playlist?list=PLjx2FAhpTe3FGbcjBbxlhf56qVR0XbVNO
 
   2. after you get to know how FK and IK work, you should be able to derive these for the Delta pick-n-place robot, for IK you can 
      check this link out:
      https://sites.google.com/site/deltarobotberkeley/how-it-works

      also there is what i used: "theory\Inverse Kinematics (Delta Robot).pdf"

   3. 3-4-5 polynomial and 4-5-6-7 polynomial point to point movement, you can learn about this in the book "Fundamentals of Robotic 
       Mechanical Systems, theory, methods, and Algorithms, Fourth Edition by Jorge Angeles - chapter 6"

   4. cubic spline (book Trajectory Planning for Automatix Machines and Robots by Luigi Biagiotti and Claudio Melchiorri)
      1. cubic spline with assigned initial and final velocities (part 4.4.1)
      2. cubic spline with assigned intial and final velocities and acceleration (part 4.4.4)
      3. smoothing cubic spline (part 4.4.5)


# step 2: PYTHON IMPLEMENTATION 

   1. INVERSE KINEMATICS:
      finding the angles of the joints according to position of the end effector --> [trajectory_planning_345.py](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/python%20implementations/trajectory_planning_345.py)

      description:
      the whole idea is built around this concept of (where is the position of end-effector? so what configuration of joints can i use to get the end-effector to that position?)

   2. POINT TO POINT MOVEMENT (3-4-5 polynomial):
      moveing from point to point in the direction desired --> [trajectory_planning_345.py](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/python%20implementations/trajectory_planning_345.py)

      description:
      the idea behind this is, we get the end-effector's position at the start and the finish, plus we get how much of the velocity do we want to actually use. 
      from positions we find the starting and final theta config. and from maximum angular velocity we find T(the time needed for the operation)
      then with the help of boundary conditions we find a 5th order polynomial called "s" as a function of time. 
      with the help of "s" polynomial we find theta and theta_dot as a function of time 

   3. POINT TO POINT MOVEMENT (4-5-6-7 polynomial):
      we repeat what we've done for sub-step 2 but with a 7th order polynomial --> [trajectory_planning_4567.py](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/python%20implementations/trajectory_planning_4567.py)

   4. CIRCLE MOVEMENT:
      cubic spline with assigned initial and final velocities --> [trajectory_planing_cubic_spline.py](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/python%20implementations/trajectory_planning_cubic_spline.py)

      description: 
  5. CIRCLE MOVEMENT:
      cubic spline with assigned initial and final velocities and acceleration --> [trajectory_planing_cubic_spline_4.4.4.py](https://github.com/ArthasMenethil-A/Delta_Robot/blob/main/python%20implementations/trajectory_planning_cubic_spline_4.4.4.py)

      description: 

