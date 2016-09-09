

#Demonstration of an equivalent dynamics controller
#The idea is to design a non-linear controller for an inverted pendulum using the equation of a normal pendulum
#One only has to find what base acceleration one needs such that the total equation of motion is equal to the one of a normal pendulum 


import math
import random
import pylab
from rk import RK4


#Constants (all in SI units)
c   = 1.0  #Linear friction term
g   = 9.81 #Gravitational acceleration
L   = 0.6  #Height to center of gravity
a_c = g

dt      = 0.001
t_total = 10

degrees  = 80.0
degToRad = math.pi/180.0
radToDeg = 1.0/degToRad
rad_init = degrees*degToRad




#Full nonlinear damped pendulum equations
ydot = lambda t, x, y: -c*y - a_c/L*math.sin(x)
xdot = lambda t, x, y: -(ydot(t, x, y) + a_c/L*math.sin(x))/c

lv = RK4(xdot, ydot)
t_nl, y_nl = lv.solve([rad_init, 0], dt, t_total)


#Linearized damped pendulum equations
ydot = lambda t, x, y: -c*y - a_c/L*x
xdot = lambda t, x, y: -(ydot(t, x, y) + a_c/L*x)/c

lv = RK4(xdot, ydot)
t_l, y_l = lv.solve([rad_init, 0], dt, t_total)

#Damping ratio
omega = math.sqrt(a_c/L)
dr = c/2.0/omega


#Plot both damped pendulum equations
pylab.plot(t_nl,y_nl[0], '-b', label='Full eq.')
pylab.plot(t_l ,y_l[0] , '-r', label='Linearized eq.')
#Plot below makes sense only if underdamped (dr < 1)
if dr < 1:
	pylab.plot(t_l ,[rad_init*math.exp(-c/2.0*x) for x in t_l], '-g', label='Envelope function')
pylab.legend(loc='upper right')
pylab.ylim(-90.0*degToRad, 90.0*degToRad)
pylab.title('Full and linearized damped pendulums ($\zeta$ = %f)' % dr)
pylab.xlabel('Time in seconds')
pylab.ylabel('Angle from vertical and radians')
pylab.show()




#Inverted pendulum with no damping or control acceleration
rad_init = 1.0*degToRad
av_init  = 5.0*degToRad
t_total  = 1.0

def ydot(t, x,y):
	return g/L*math.sin(x)
def xdot(t, x,y):
	return y

lv = RK4(xdot, ydot)
t, y = lv.solve([rad_init, av_init], dt, t_total)

pylab.plot(t, y[0])
pylab.ylim(-90.0*degToRad, 90.0*degToRad)
pylab.title('Pendulum without damping')
pylab.xlabel('Time in seconds')
pylab.ylabel('Angle from vertical and radians')
pylab.show()




#Inverted pendulum with PD control
t_sampling = 0.1
num_loops  = 100
stdDev     = 0.5*degToRad

P = -17.0
D = -2.4

theta1 = 0.0
theta2 = 0.0

#Control law
def a_control():
	return P*theta2 + D*(theta2-theta1)/t_sampling

#Equations of motion
def ydot(t, x,y):
	return g/L*math.sin(x) + a_control()/L*math.cos(x)
def xdot(t, x,y):
	return y


#Solve equations of motion with control acceleration, sampling every t_sampling seconds
x_init = 50.0*degToRad
y_init = 10.0*degToRad

x = x_init
y = y_init 

totalTime  = []
totalTheta = []

for n in range(0,num_loops):
	#Define Runge-Kutta
	lv = RK4(xdot, ydot)

	#Solve motion for t_sampling seconds
	t, theta = lv.solve([x, y], dt, t_sampling)

	#Append data to total time and angle
	totalTime.extend([x+t_sampling*n for x in t])
	totalTheta.extend(theta[0])

	#Change initial angle and angular velocity for next loop
	x = theta[0][-1]
	y = theta[1][-1]

	#Update angle values for controller
	if n == 0:
		theta2 = x + random.normalvariate(0.0, stdDev)
	else:
		theta1 = theta2
		theta2 = x + random.normalvariate(0.0, stdDev)




#Inverted pendulum with equivalent dynamics controller
x = x_init
y = y_init 

theta1 = 0.0
theta2 = 0.0

totalTheta2 = []

a_c = 0.6*g
c   = 4.0

#Control law
def a_control2():
	return -(g+a_c)*math.tan(theta2) - L*c*(theta2-theta1)/t_sampling/math.cos(theta2)

#Equations of motion
def ydot2(t, x,y):
	return g/L*math.sin(x) + a_control2()/L*math.cos(x)
def xdot2(t, x,y):
	return y


for n in range(0,num_loops):
	#Define Runge-Kutta
	lv = RK4(xdot2, ydot2)

	#Solve motion for t_sampling seconds
	t, theta = lv.solve([x, y], dt, t_sampling)

	#Append data to total time and angle
	totalTheta2.extend(theta[0])

	#Change initial angle and angular velocity for next loop
	x = theta[0][-1]
	y = theta[1][-1]

	#Update angle values for controller
	if n == 0:
		theta2 = x + random.normalvariate(0.0, stdDev)
	else:
		theta1 = theta2
		theta2 = x + random.normalvariate(0.0, stdDev)


#Plot results
pylab.plot(totalTime,totalTheta , '-b', label='PD controller')
pylab.plot(totalTime,totalTheta2, '-r', label='Equivalent dynamics controller')
pylab.legend(loc='upper right')
pylab.ylim(-90.0*degToRad, 90.0*degToRad)
pylab.title('Pendulum with PD control using base acceleration')
pylab.xlabel('Time in seconds')
pylab.ylabel('Angle from vertical and radians')
pylab.show()


