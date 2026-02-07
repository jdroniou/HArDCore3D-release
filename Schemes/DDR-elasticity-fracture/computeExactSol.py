#! python
import sympy as sym
#from sympy.printing import ccode
from sympy import *

x, y, z, pi, la = sym.symbols("x.x() x.y() x.z() PI lambda");

# Trigonometric function
# u1 = sym.sin(2 * pi * x) * sym.cos(2 * pi * y) * sym.cos(2 * pi * z);
# u2 = 0.5 * sym.cos(2 * pi * x) * sym.sin(2 * pi * y) * sym.cos(2 * pi * z);
# u3 = - sym.cos(2 * pi * x) * sym.cos(2 * pi * y) * sym.sin(2 * pi * z);

# u1 = -sym.sin(0.5*pi*x)*sym.cos(0.5*pi*y)*pow(z,2)*(1/la);
# u2 = pow(z,2)*(1/la)+3*pow(x,2)*sym.sin(z);
# u3 = pow(x,2)*pow(z,2)*(1/la);

# u1 = 2*cos((pi/2)*x)*pow(z,4)*(1/la);
# u2 = 8*cos((pi/2)*x)*pow(z,3)*(1/la)+3*pow(x,2)*sym.sin(z);
# u3 = -(16/pi)*sin((pi/2)*x)*pow(z,3)*(1/la);

# u1 = (-sym.sin(x)*sym.cos(y)*pow(z,2)-y)*(1/la);
# u2 = pow(z,2)*(1/la)+3*pow(x,2)*sym.sin(z);
# u3 = pow(x,2)*pow(z,2)*(1/la);

u1 = (-(0.25)*sym.sin(x)*sym.cos(y)*pow(z,2)-y)*(1/la);
u2 = (0.25)*pow(z,2)*(1/la)+3*pow(x,2)*sym.sin(z);
u3 = (0.25)*pow(x,2)*pow(z,2)*(1/la);


# Linear
#u1 = x;
#u2 = -y;
#u3 = z;

#Quadratic
#u1 = x*x;
#u2 = z*y;
#u3 = x*z+z*z;

# Computations
u = [u1, u2, u3];
grad_u1 = [sym.diff(u1,x), sym.diff(u1,y), sym.diff(u1,z)];
grad_u2 = [sym.diff(u2,x), sym.diff(u2,y), sym.diff(u2,z)];
grad_u3 = [sym.diff(u3,x), sym.diff(u3,y), sym.diff(u3,z)];
grad_u = [grad_u1, grad_u2, grad_u3];

lapl_u1 = sym.diff(grad_u1[0],x) + sym.diff(grad_u1[1],y) + sym.diff(grad_u1[2],z);
lapl_u2 = sym.diff(grad_u2[0],x) + sym.diff(grad_u2[1],y) + sym.diff(grad_u2[2],z);
lapl_u3 = sym.diff(grad_u3[0],x) + sym.diff(grad_u3[1],y) + sym.diff(grad_u3[2],z);

div_u = sym.diff(u1,x) + sym.diff(u2,y) + sym.diff(u3,z);
grad_div_u = [sym.diff(div_u,x),
                sym.diff(div_u,y),
                sym.diff(div_u,z)];

div_eps_u = [0.5*(lapl_u1 + grad_div_u[0]), 0.5*(lapl_u2 + grad_div_u[1]), 0.5*(lapl_u3 + grad_div_u[2])];


# Output
print("u=", ccode(u), "\n");
print("grad u1=", ccode(grad_u1), "\n");
print("grad u2=", ccode(grad_u2), "\n");
print("grad u3=", ccode(grad_u3), "\n");
print("div_eps_u=", ccode(div_eps_u), "\n");
print("grad_div_u=", ccode(grad_div_u), "\n");
#print("grad DIS=", ccode(grad_DIS), "\n");



