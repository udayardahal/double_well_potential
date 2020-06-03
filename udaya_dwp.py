"""
===============================================================================
This is double well potential used in constant pH molecular dynamics simulation.
This code helps you to interactively choose the value to meet your requirement. 
The barrier height, location of potential well and the shape of the well to the 
barrier is important to sample the protonation and deprotonation state and avoid 
unphysical state. The parameters of the double well potential are different for 
different chemical composition.
This potential is a function of U(l, k, d, b, a, w, m, s, r)
For the equation see the paper: J.  Chem.  Theory  Comput.2016, 12, 1040−1051
===============================================================================

Using the slider widget to control visual properties of your plot.

"""
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.55)

#k = np.arange(0,20,0.05)
#d = np.arange(0,20,0.05)
#b = np.arange(0,1,0.0005)
#a = np.arange(0,1,0.001)
#w = np.arange(0,500,1) 
#m = np.arange(0,1,0.01)
#s = np.arange(0,1,0.1)
#r = np.arange(0,50,0.5);
l = np.arange(-0.1,1.1,0.01)
# k = 4.417091;
# d = 3.5;
# b = 0.002957;
# a = 0.042082;
# w = 200.0; 
# m = 0.150672;
# s = 0.3;
# r = 16.457502;

k0 = 4.417091;
d0 = 3.5;
b0 = 0.002957;
a0 = 0.042082;
w0 = 200.0; 
m0 = 0.150672;
s0 = 0.3;
r0 = 16.457502;



def U(l,k = 4.417091,d = 3.5,b = 0.002957,a = 0.042082,w = 200.0,m = 0.150672,s = 0.3,r = 16
.457502):
        U_dwp_lambda = -k * ( np.exp( (- (l - 1 - b)**2) / (2 * a**2)) + np.exp( (- (l + b)*
*2) / (2 * a**2)) ) + d * ( np.exp( (- (l - 0.5)**2) / (2 * s**2)) ) + 0.5 * w * ( ( 1 - mat
h.erf( r * (l + m))) + (1 + math.erf( r * (l -1 - m))) )
        return U_dwp_lambda

U_dwp_lambda = np.array([U(mm) for mm in l])

ll, = plt.plot(l, U_dwp_lambda, lw=2)
plt.ylim([-10,50])
plt.xlim([-1,2])

delta_k = 0.05
delta_d = 0.05
delta_b = 0.0005
delta_a = 0.001
delta_w = 1
delta_m = 0.001
delta_s = 0.01
delta_r = 0.5
delta_l = 0.01
ax.margins(x=0)

axcolor = 'lightgoldenrodyellow'
axk = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
axd = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
axb = plt.axes([0.25, 0.20, 0.65, 0.03], facecolor=axcolor)
axa = plt.axes([0.25, 0.25, 0.65, 0.03], facecolor=axcolor)
axw = plt.axes([0.25, 0.30, 0.65, 0.03], facecolor=axcolor)
axm = plt.axes([0.25, 0.35, 0.65, 0.03], facecolor=axcolor)
axs = plt.axes([0.25, 0.40, 0.65, 0.03], facecolor=axcolor)
axr = plt.axes([0.25, 0.45, 0.65, 0.03], facecolor=axcolor)
#axl = plt.axes([0.25, 0.50, 0.65, 0.03], facecolor=axcolor)

sk = Slider(axk, 'k',    0,   50, valinit=k0, valstep=delta_k) 
sd = Slider(axd, 'd',    0,   60, valinit=d0, valstep=delta_d) 
sb = Slider(axb, 'b',    0,    1, valinit=b0, valstep=delta_b) 
sa = Slider(axa, 'a',    0,    1, valinit=a0, valstep=delta_a) 
sw = Slider(axw, 'w',    0,  500, valinit=w0, valstep=delta_w) 
sm = Slider(axm, 'm',    0,    1, valinit=m0, valstep=delta_m) 
ss = Slider(axs, 's',    0,    1, valinit=s0, valstep=delta_s) 
sr = Slider(axr, 'r',    0,   50, valinit=r0, valstep=delta_r) 
#sl = Slider(axl, 'l', 0.05, 1.05, valinit=l0, valstep=delta_l) 

def update(val):
   k = sk.val
   d = sd.val
   b = sb.val
   a = sa.val
   w = sw.val
   m = sm.val
   s = ss.val
   r = sr.val
#   l = sl.val
   ll.set_ydata(np.asarray([U(mm,k,d,b,a,w,m,s,r) for mm in l]))

   fig.canvas.draw_idle()

sk.on_changed(update)  
sd.on_changed(update) 
sb.on_changed(update) 
sa.on_changed(update) 
sw.on_changed(update) 
sm.on_changed(update) 
ss.on_changed(update) 
sr.on_changed(update) 
#sl.on_changed(update) 

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    sk.reset()  
    sd.reset()
    sb.reset()
    sa.reset()
    sw.reset()
    sm.reset()
    ss.reset()
    sr.reset()
#    sl.reset()

button.on_clicked(reset)

rax = plt.axes([0.025, 0.6, 0.15, 0.3], facecolor=axcolor)
radio = RadioButtons(rax, ('red', 'blue', 'green', 'black', 'orange', 'violet', 'magenta'), 
active = 0)


def colorfunc(label):
    ll.set_color(label)
    fig.canvas.draw_idle()
radio.on_clicked(colorfunc)

plt.show()
