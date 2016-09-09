# Equivalent-dynamics-controller
The python code demonstrates a simple non-linear control law for an inverted pendulum using an idea called 'Equivalent dynamics' (ED). The idea is the following; if one has the equations of motion for an [inverted pendulum with a control force] (https://en.wikipedia.org/wiki/Inverted_pendulum) (at the base in this case), and the equations of motion for the same mechanism with another set of forces that cause it to be stabe (i.e a normal [damped pendulum] (http://www2.phy.ilstu.edu/~rfm/380s13/chapters/CH2.4_Pendulum.pdf)), one can simply set them equal to each other and rearrange to get the needed control acceleration to make the inverted pendulum act as a damped pendulum.

I wrote the python code to demonstrate that not only the idea works, but it also stable for a wider set of initial disturbed conditions compared to a simple PD controller. Another feature of this control law is that for [small angles] (https://en.wikipedia.org/wiki/Small-angle_approximation) it is equivalent to a PD controller. The weights of the controller can be estimated directly from the equivalent linearized damped pendulum system, to give the required settling time.

The code uses a [costum Runge-Kutta solver] (https://github.com/sbillaudelle/runge-kutta) to solve the equations of motion in between the controller sampling times. After each iteration the pendulum angle is recorded with some normal noise and the control acceleration is updated. The pendulum's angular speed is estimated (for now) using a first order backward difference.

When the code is run it outputs the following graphs:
* A graph for a damped pendulum and a linearized damped pendulum with its envelope function, to show how fast its motion decays. One can see how the actual and linearized motion diverges over time. This graph can be used to determine what response one wants from the ED controller.
* A graph of the motion of an inverted pendulum without a controller, to get an idea of the time scales of the problem, etc.
* A graph of the same inverted pendulum, being controlled by a PD and a ED controller. The PD controller has very similar weights to the ED controller for small angles. One can increase the initial pendulum angle and see that the ED controller can function even with large initial pendulum angles.
