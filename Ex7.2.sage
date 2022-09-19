################################################################################

R.<x> = QQ[]
X = HyperellipticCurve(R([5, 0, -4, 10, -8, 0, 1])) # LMFDB 3950.b.39500.1
f = X.hyperelliptic_polynomials()[0]
R, S = X(1,2), X(1,-2)
K = Qp(5,12)
XK = X.change_ring(K)
L.<a> = K.ext(x^4-5) # ramified extension
XL = X.change_ring(L)
P1 = XL.lift_x(a,all)[1] # ref pt in e_1
P2 = XL.lift_x(a,all)[0] # ref pt in e_2
P3 = XL.lift_x(a+2,all)[1] # ref pt in e_3
P4 = XL.lift_x(a+2,all)[0] # ref pt in e_4
P5 = XL.lift_x(a+3,all)[0] # ref pt in e_5
P6 = XL.lift_x(a+3,all)[1] # ref pt in e_6
def Log(z): return L(z).log(p_branch = 0, change_frac = True) # Iwasawa branch

################################################################################

fK = f.change_ring(K)
fac1 = list(fK.factor())[2][0] # fac1 modp = x^2
fac2 = list(fK.factor())[0][0] # fac2 modp = (x - 2)^2
fac3 = list(fK.factor())[1][0] # fac3 modp = (x - 3)^2
def coeff(f):
    B = f(0)
    A = f(1) - B - 1
    C = B - A^2/4
    return A, B, C
A1, B1, C1 = coeff(fac1)
A2, B2, C2 = coeff(fac2)
A3, B3, C3 = coeff(fac3)

################################################################################

N = 32 # truncation level
rang = 5
TT.<t1> = PowerSeriesRing(K, 't1', default_prec = N)
l1 = t1*(1+C1*t1^2)^(-1/2) # be careful with t1: t1 --> 1/(x + A1/2)
TT.<t2> = PowerSeriesRing(K, 't2', default_prec = N)
l2 = t2*(1+C2*t2^2)^(-1/2) # be careful with t2: t2 --> 1/(x + A2/2)
TT.<t3> = PowerSeriesRing(K, 't3', default_prec = N)
l3 = t3*(1+C3*t3^2)^(-1/2) # be careful with t3: t3 --> 1/(x + A3/2)

################################################################################

# 0) The components v_plus and v_minus
TT.<t1,t2,t3> = PowerSeriesRing(K, 't1,t2,t3', default_prec = 3*N)
l = TT(l1)*TT(l2)*TT(l3)

################################################################################

# first pole reduction : x-coordinate is -A1/2
TT.<t> = PowerSeriesRing(K, 't', default_prec = 2*N)
w = [l.truncate()(1/t,1/(t-A1/2+A2/2),1/(t-A1/2+A3/2))*(t-A1/2)^i/2 for i in range(rang)] # original forms around the first pole
d1 = [w[i].residue() for i in range(rang)] # coefficients of dx/(x+A1/2)
w = [w[i] - d1[i]*(1/t) for i in range(rang)]
ED = [0,0] + [-i/(t^(i+1)) for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
RR.<u> = QQ[]
F1 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F1[i] += FF[i][j]/(u^(N-2-j)) # be careful with the first coordinate: u --> x + A1/2

################################################################################

# test
[(l.truncate()(1/t,1/(t-A1/2+A2/2),1/(t-A1/2+A3/2))*(t-A1/2)^i/2 - F1[i](t).derivative() - d1[i]*(1/t)).valuation() for i in range(rang)]
### [0, 0, 0, 0, 0]

################################################################################

# second pole reduction : x-coordinate is -A2/2
w = [l.truncate()(1/(t-A2/2+A1/2),1/t,1/(t-A2/2+A3/2))*(t-A2/2)^i/2 - F1[i](t-A2/2+A1/2).derivative() - d1[i]*(1/(t-A2/2+A1/2)) for i in range(rang)]
d2 = [w[i].residue() for i in range(rang)] # coefficients of dx/(x+A2/2)
w = [w[i] - d2[i]*(1/t) for i in range(rang)]
ED = [0,0] + [-i/(t^(i+1)) for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
RR.<u> = QQ[]
F2 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F2[i] += FF[i][j]/(u^(N-2-j)) # be careful with the first coordinate: u --> x + A2/2

################################################################################

# test
[(l.truncate()(1/(t-A2/2+A1/2),1/t,1/(t-A2/2+A3/2))*(t-A2/2)^i/2 - F1[i](t-A2/2+A1/2).derivative() - d1[i]*(1/(t-A2/2+A1/2)) - F2[i](t).derivative() - d2[i]*(1/t)).valuation()  for i in range(rang)]
### [0, 0, 0, 0, 0]

################################################################################

# third pole reduction : x-coordinate is -A3/2
w = [l.truncate()(1/(t-A3/2+A1/2),1/(t-A3/2+A2/2),1/t)*(t-A3/2)^i/2 - F1[i](t-A3/2+A1/2).derivative() - d1[i]*(1/(t-A3/2+A1/2)) - F2[i](t-A3/2+A2/2).derivative() - d2[i]*(1/(t-A3/2+A2/2)) for i in range(rang)]
d3 = [w[i].residue() for i in range(rang)] # coefficients of dx/(x+A3/2)
w = [w[i] - d3[i]*(1/t) for i in range(rang)]
ED = [0,0] + [-i/(t^(i+1)) for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
RR.<u> = QQ[]
F3 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F3[i] += FF[i][j]/(u^(N-2-j)) # be careful with the first coordinate: u --> x + A3/2

################################################################################

# test
[(l.truncate()(1/(t-A3/2+A1/2),1/(t-A3/2+A2/2),1/t)*(t-A3/2)^i/2 - F1[i](t-A3/2+A1/2).derivative() - d1[i]*(1/(t-A3/2+A1/2)) - F2[i](t-A3/2+A2/2).derivative() - d2[i]*(1/(t-A3/2+A2/2)) - F3[i](t).derivative() - d3[i]*(1/t)).valuation() for i in range(rang)]
### [+Infinity, +Infinity, +Infinity, 0, 0]

################################################################################

dinf30 = (l.truncate()(1/(t-A3/2+A1/2),1/(t-A3/2+A2/2),1/t)*(t-A3/2)^3/2 - F1[3](t-A3/2+A1/2).derivative() - d1[3]*(1/(t-A3/2+A1/2)) - F2[3](t-A3/2+A2/2).derivative() - d2[3]*(1/(t-A3/2+A2/2)) - F3[3](t).derivative() - d3[3]*(1/t)).list()[0]
dinf41 = (l.truncate()(1/(t-A3/2+A1/2),1/(t-A3/2+A2/2),1/t)*(t-A3/2)^4/2 - F1[4](t-A3/2+A1/2).derivative() - d1[4]*(1/(t-A3/2+A1/2)) - F2[4](t-A3/2+A2/2).derivative() - d2[4]*(1/(t-A3/2+A2/2)) - F3[4](t).derivative() - d3[4]*(1/t)).list()[1]
dinf40 = (l.truncate()(1/(t-A3/2+A1/2),1/(t-A3/2+A2/2),1/t)*(t-A3/2)^4/2 - F1[4](t-A3/2+A1/2).derivative() - d1[4]*(1/(t-A3/2+A1/2)) - F2[4](t-A3/2+A2/2).derivative() - d2[4]*(1/(t-A3/2+A2/2)) - F3[4](t).derivative() - d3[4]*(1/t)).list()[0] + dinf41*A3/2

################################################################################

# test
[l.truncate()(1/(t-A3/2+A1/2),1/(t-A3/2+A2/2),1/t)*(t-A3/2)^i/2 == F1[i](t-A3/2+A1/2).derivative() + d1[i]*(1/(t-A3/2+A1/2)) + F2[i](t-A3/2+A2/2).derivative() + d2[i]*(1/(t-A3/2+A2/2)) + F3[i](t).derivative() + d3[i]*(1/t) for i in range(3)]
l.truncate()(1/(t-A3/2+A1/2),1/(t-A3/2+A2/2),1/t)*(t-A3/2)^3/2 == F1[3](t-A3/2+A1/2).derivative() + d1[3]*(1/(t-A3/2+A1/2)) + F2[3](t-A3/2+A2/2).derivative() + d2[3]*(1/(t-A3/2+A2/2)) + F3[3](t).derivative() + d3[3]*(1/t) + dinf30
l.truncate()(1/(t-A3/2+A1/2),1/(t-A3/2+A2/2),1/t)*(t-A3/2)^4/2 == F1[4](t-A3/2+A1/2).derivative() + d1[4]*(1/(t-A3/2+A1/2)) + F2[4](t-A3/2+A2/2).derivative() + d2[4]*(1/(t-A3/2+A2/2)) + F3[4](t).derivative() + d3[4]*(1/t) + dinf40 + dinf41*(t-A3/2)
### [True, True, True]
### True
### True

################################################################################

# RESULT:
# for i = 0,1,2, omega_i = dF1[i] + d1[i]*dx/(x+A1/2) + dF2[i] + d2[i]*dx/(x+A2/2) + dF3[i] + d3[i]*dx/(x+A3/2)
# omega_3 = dF1[3] + d1[3]*dx/(x+A1/2) + dF2[3] + d2[3]*dx/(x+A2/2) + dF3[3] + d3[3]*dx/(x+A3/2) + dinf30*dx
# omega_4 = dF1[4] + d1[4]*dx/(x+A1/2) + dF2[4] + d2[4]*dx/(x+A2/2) + dF3[4] + d3[4]*dx/(x+A3/2) + dinf40*dx + dinf41*x*dx

################################################################################

# 1) The component v_1
TT.<t2,t3> = PowerSeriesRing(K, 't2,t3', default_prec = 3*N)
l = TT(l2)*TT(l3)
g = fac1
gder = g.derivative()
TT.<t> = PowerSeriesRing(K, 't', default_prec = 2*N)

################################################################################

# first pole reduction : x-coordinate is -A2/2
x = t - A2/2 # x-local-coordinate at the pole
y = g(x).sqrt() # y-local-coordinate at the pole
w = [l.truncate()(1/t,1/(t-A2/2+A3/2))*(t-A2/2)^i*x.derivative()/(2*y) for i in range(rang)] # original forms around the first pole
d12 = [w[i].residue()/(1/t*x.derivative()/(2*y)).residue() for i in range(rang)] # coefficients of 1/(x+A2/2)*dx/2y
w = [w[i] - d12[i]*(1/t*x.derivative()/(2*y)) for i in range(rang)]
ED = [0,0] + [(t*gder(x)-2*i*g(x))/(t^(i+1))*x.derivative()/(2*y) for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
RR.<u,v> = QQ[]
F12 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F12[i] += FF[i][j]*v/(u^(N-2-j)) # be careful with the first coordinate: u --> x + A2/2

################################################################################

# test
[(l.truncate()(1/t,1/(t-A2/2+A3/2))*(t-A2/2)^i*x.derivative()/(2*y) - F12[i](t,y).derivative() - d12[i]*(1/t*x.derivative()/(2*y))).valuation() for i in range(rang)]
### [0, 0, 0, 0, 0]

################################################################################

# second pole reduction : x-coordinate is -A3/2
x = t - A3/2 # x-local-coordinate at the pole
y = g(x).sqrt() # y-local-coordinate at the pole
w = [l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^i*x.derivative()/(2*y) - F12[i](t-A3/2+A2/2,y).derivative() - d12[i]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) for i in range(rang)]
d13 = [w[i].residue()/(1/t*x.derivative()/(2*y)).residue() for i in range(rang)] # coefficients of 1/(x+A3/2)*dx/2y
w = [w[i] - d13[i]*(1/t*x.derivative()/(2*y)) for i in range(rang)]
ED = [0,0] + [(t*gder(x)-2*i*g(x))/(t^(i+1))*x.derivative()/(2*y) for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
RR.<u,v> = QQ[]
F13 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F13[i] += FF[i][j]*v/(u^(N-2-j)) # be careful with the first coordinate: u --> x + A3/2

################################################################################

# test
[(l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^i*x.derivative()/(2*y) - F12[i](t-A3/2+A2/2,y).derivative() - d12[i]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) - F13[i](t,y).derivative() - d13[i]*(1/t*x.derivative()/(2*y))).valuation() for i in range(rang)]
### [+Infinity, +Infinity, 0, 0, 0]

################################################################################

dinf120 = ((2*y)*(l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^2*x.derivative()/(2*y) - F12[2](t-A3/2+A2/2,y).derivative() - d12[2]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) - F13[2](t,y).derivative() - d13[2]*(1/t*x.derivative()/(2*y)))).list()[0]
dinf131 = ((2*y)*(l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^3*x.derivative()/(2*y) - F12[3](t-A3/2+A2/2,y).derivative() - d12[3]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) - F13[3](t,y).derivative() - d13[3]*(1/t*x.derivative()/(2*y)))).list()[1]
dinf130 = ((2*y)*(l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^3*x.derivative()/(2*y) - F12[3](t-A3/2+A2/2,y).derivative() - d12[3]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) - F13[3](t,y).derivative() - d13[3]*(1/t*x.derivative()/(2*y)))).list()[0] + dinf131*A3/2
dinf142 = ((2*y)*(l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^4*x.derivative()/(2*y) - F12[4](t-A3/2+A2/2,y).derivative() - d12[4]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) - F13[4](t,y).derivative() - d13[4]*(1/t*x.derivative()/(2*y)))).list()[2]
dinf141 = ((2*y)*(l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^4*x.derivative()/(2*y) - F12[4](t-A3/2+A2/2,y).derivative() - d12[4]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) - F13[4](t,y).derivative() - d13[4]*(1/t*x.derivative()/(2*y)))).list()[1] + dinf142*A3
dinf140 = ((2*y)*(l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^4*x.derivative()/(2*y) - F12[4](t-A3/2+A2/2,y).derivative() - d12[4]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) - F13[4](t,y).derivative() - d13[4]*(1/t*x.derivative()/(2*y)))).list()[0] + dinf141*A3/2 - dinf142*A3^2/4

################################################################################

[l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^i*x.derivative()/(2*y) == F12[i](t-A3/2+A2/2,y).derivative() + d12[i]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) + F13[i](t,y).derivative() + d13[i]*(1/t*x.derivative()/(2*y)) for i in range(2)]
l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^2*x.derivative()/(2*y) == F12[2](t-A3/2+A2/2,y).derivative() + d12[2]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) + F13[2](t,y).derivative() + d13[2]*(1/t*x.derivative()/(2*y)) + dinf120*x.derivative()/(2*y)
l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^3*x.derivative()/(2*y) == F12[3](t-A3/2+A2/2,y).derivative() + d12[3]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) + F13[3](t,y).derivative() + d13[3]*(1/t*x.derivative()/(2*y)) + dinf130*x.derivative()/(2*y) + dinf131*x*x.derivative()/(2*y)
l.truncate()(1/(t-A3/2+A2/2),1/t)*(t-A3/2)^4*x.derivative()/(2*y) == F12[4](t-A3/2+A2/2,y).derivative() + d12[4]*(1/(t-A3/2+A2/2)*x.derivative()/(2*y)) + F13[4](t,y).derivative() + d13[4]*(1/t*x.derivative()/(2*y)) + dinf140*x.derivative()/(2*y) + dinf141*x*x.derivative()/(2*y) + dinf142*x^2*x.derivative()/(2*y)
### [True, True]
### True
### True
### True

################################################################################

# RESULT:
# For i = 0,1, omega_i = dF12[i] + d12[i]*1/(x+A2/2)*dx/2y + dF13[i] + d13[i]*1/(x+A3/2)*dx/2y
# omega_2 = dF12[2] + d12[2]*1/(x+A2/2)*dx/2y + dF13[2] + d13[2]*1/(x+A3/2)*dx/2y + dinf120*dx/2y
# omega_3 = dF12[3] + d12[3]*1/(x+A2/2)*dx/2y + dF13[3] + d13[3]*1/(x+A3/2)*dx/2y + dinf130*dx/2y + dinf131*x*dx/2y
# omega_4 = dF12[4] + d12[4]*1/(x+A2/2)*dx/2y + dF13[4] + d13[4]*1/(x+A3/2)*dx/2y + dinf140*dx/2y + dinf141*x*dx/2y + dinf142*x^2*dx/2y

################################################################################

# 2) The component v_2
TT.<t1,t3> = PowerSeriesRing(K, 't1,t3', default_prec = 3*N)
l = TT(l1)*TT(l3)
g = fac2
gder = g.derivative()
TT.<t> = PowerSeriesRing(K, 't', default_prec = 2*N)

################################################################################

# first pole reduction : x-coordinate is -A1/2
x = t - A1/2 # x-local-coordinate at the pole
y = g(x).sqrt() # y-local-coordinate at the pole
w = [l.truncate()(1/t,1/(t-A1/2+A3/2))*(t-A1/2)^i*x.derivative()/(2*y) for i in range(rang)] # original forms around the first pole
d21 = [w[i].residue()/(1/t*x.derivative()/(2*y)).residue() for i in range(rang)] # coefficients of 1/(x+A1/2)*dx/2y
w = [w[i] - d21[i]*(1/t*x.derivative()/(2*y)) for i in range(rang)]
ED = [0,0] + [(t*gder(x)-2*i*g(x))/(t^(i+1))*x.derivative()/(2*y) for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
RR.<u,v> = QQ[]
F21 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F21[i] += FF[i][j]*v/(u^(N-2-j)) # be careful with the first coordinate: u --> x + A1/2

################################################################################

# test
[(l.truncate()(1/t,1/(t-A1/2+A3/2))*(t-A1/2)^i*x.derivative()/(2*y) - F21[i](t,y).derivative() - d21[i]*(1/t*x.derivative()/(2*y))).valuation() for i in range(rang)]
### [0, 0, 0, 0, 0]

################################################################################

# second pole reduction : x-coordinate is -A3/2
x = t - A3/2 # x-local-coordinate at the pole
y = g(x).sqrt() # y-local-coordinate at the pole
w = [l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^i*x.derivative()/(2*y) - F21[i](t-A3/2+A1/2,y).derivative() - d21[i]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) for i in range(rang)]
d23 = [w[i].residue()/(1/t*x.derivative()/(2*y)).residue() for i in range(rang)] # coefficients of 1/(x+A3/2)*dx/2y
w = [w[i] - d23[i]*(1/t*x.derivative()/(2*y)) for i in range(rang)]
ED = [0,0] + [(t*gder(x)-2*i*g(x))/(t^(i+1))*x.derivative()/(2*y) for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
RR.<u,v> = QQ[]
F23 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F23[i] += FF[i][j]*v/(u^(N-2-j)) # be careful with the first coordinate: u --> x + A3/2

################################################################################

# test
[(l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^i*x.derivative()/(2*y) - F21[i](t-A3/2+A1/2,y).derivative() - d21[i]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) - F23[i](t,y).derivative() - d23[i]*(1/t*x.derivative()/(2*y))).valuation() for i in range(rang)]
### [+Infinity, +Infinity, 0, 0, 0]

################################################################################

dinf220 = ((2*y)*(l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^2*x.derivative()/(2*y) - F21[2](t-A3/2+A1/2,y).derivative() - d21[2]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) - F23[2](t,y).derivative() - d23[2]*(1/t*x.derivative()/(2*y)))).list()[0]
dinf231 = ((2*y)*(l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^3*x.derivative()/(2*y) - F21[3](t-A3/2+A1/2,y).derivative() - d21[3]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) - F23[3](t,y).derivative() - d23[3]*(1/t*x.derivative()/(2*y)))).list()[1]
dinf230 = ((2*y)*(l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^3*x.derivative()/(2*y) - F21[3](t-A3/2+A1/2,y).derivative() - d21[3]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) - F23[3](t,y).derivative() - d23[3]*(1/t*x.derivative()/(2*y)))).list()[0] + dinf231*A3/2
dinf242 = ((2*y)*(l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^4*x.derivative()/(2*y) - F21[4](t-A3/2+A1/2,y).derivative() - d21[4]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) - F23[4](t,y).derivative() - d23[4]*(1/t*x.derivative()/(2*y)))).list()[2]
dinf241 = ((2*y)*(l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^4*x.derivative()/(2*y) - F21[4](t-A3/2+A1/2,y).derivative() - d21[4]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) - F23[4](t,y).derivative() - d23[4]*(1/t*x.derivative()/(2*y)))).list()[1] + dinf242*A3
dinf240 = ((2*y)*(l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^4*x.derivative()/(2*y) - F21[4](t-A3/2+A1/2,y).derivative() - d21[4]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) - F23[4](t,y).derivative() - d23[4]*(1/t*x.derivative()/(2*y)))).list()[0] + dinf241*A3/2 - dinf242*A3^2/4

################################################################################

# test
[l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^i*x.derivative()/(2*y) == F21[i](t-A3/2+A1/2,y).derivative() + d21[i]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) + F23[i](t,y).derivative() + d23[i]*(1/t*x.derivative()/(2*y)) for i in range(2)]
l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^2*x.derivative()/(2*y) == F21[2](t-A3/2+A1/2,y).derivative() + d21[2]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) + F23[2](t,y).derivative() + d23[2]*(1/t*x.derivative()/(2*y)) + dinf220*x.derivative()/(2*y)
l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^3*x.derivative()/(2*y) == F21[3](t-A3/2+A1/2,y).derivative() + d21[3]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) + F23[3](t,y).derivative() + d23[3]*(1/t*x.derivative()/(2*y)) + dinf230*x.derivative()/(2*y) + dinf231*x*x.derivative()/(2*y)
l.truncate()(1/(t-A3/2+A1/2),1/t)*(t-A3/2)^4*x.derivative()/(2*y) == F21[4](t-A3/2+A1/2,y).derivative() + d21[4]*(1/(t-A3/2+A1/2)*x.derivative()/(2*y)) + F23[4](t,y).derivative() + d23[4]*(1/t*x.derivative()/(2*y)) + dinf240*x.derivative()/(2*y) + dinf241*x*x.derivative()/(2*y) + dinf242*x^2*x.derivative()/(2*y)
### [True, True]
### True
### True
### True

################################################################################

# RESULT:
# For i = 0,1, omega_i = dF21[i] + d21[i]*1/(x+A1/2)*dx/2y + dF23[i] + d23[i]*1/(x+A3/2)*dx/2y
# omega_2 = dF21[2] + d21[2]*1/(x+A1/2)*dx/2y + dF23[2] + d23[2]*1/(x+A3/2)*dx/2y + dinf220*dx/2y
# omega_3 = dF21[3] + d21[3]*1/(x+A1/2)*dx/2y + dF23[3] + d23[3]*1/(x+A3/2)*dx/2y + dinf230*dx/2y + dinf231*x*dx/2y
# omega_4 = dF21[4] + d21[4]*1/(x+A1/2)*dx/2y + dF23[4] + d23[4]*1/(x+A3/2)*dx/2y + dinf240*dx/2y + dinf241*x*dx/2y + dinf242*x^2*dx/2y

################################################################################

# 3) The component v_3
TT.<t1,t2> = PowerSeriesRing(K, 't1,t2', default_prec = 3*N)
l = TT(l1)*TT(l2)
g = fac3
gder = g.derivative()
TT.<t> = PowerSeriesRing(K, 't', default_prec = 2*N)

################################################################################

# first pole reduction : x-coordinate is -A1/2
x = t - A1/2 # x-local-coordinate at the pole
y = g(x).sqrt() # y-local-coordinate at the pole
w = [l.truncate()(1/t,1/(t-A1/2+A2/2))*(t-A1/2)^i*x.derivative()/(2*y) for i in range(rang)] # original forms around the first pole
d31 = [w[i].residue()/(1/t*x.derivative()/(2*y)).residue() for i in range(rang)] # coefficients of 1/(x+A1/2)*dx/2y
w = [w[i] - d31[i]*(1/t*x.derivative()/(2*y)) for i in range(rang)]
ED = [0,0] + [(t*gder(x)-2*i*g(x))/(t^(i+1))*x.derivative()/(2*y) for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
RR.<u,v> = QQ[]
F31 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F31[i] += FF[i][j]*v/(u^(N-2-j)) # be careful with the first coordinate: u --> x + A1/2

################################################################################

# test
[(l.truncate()(1/t,1/(t-A1/2+A2/2))*(t-A1/2)^i*x.derivative()/(2*y) - F31[i](t,y).derivative() - d31[i]*(1/t*x.derivative()/(2*y))).valuation() for i in range(rang)]
### [0, 0, 0, 0, 0]

################################################################################

# second pole reduction : x-coordinate is -A2/2
x = t - A2/2 # x-local-coordinate at the pole
y = g(x).sqrt() # y-local-coordinate at the pole
w = [l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^i*x.derivative()/(2*y) - F31[i](t-A2/2+A1/2,y).derivative() - d31[i]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) for i in range(rang)]
d32 = [w[i].residue()/(1/t*x.derivative()/(2*y)).residue() for i in range(rang)] # coefficients of 1/(x+A2/2)*dx/2y
w = [w[i] - d32[i]*(1/t*x.derivative()/(2*y)) for i in range(rang)]
ED = [0,0] + [(t*gder(x)-2*i*g(x))/(t^(i+1))*x.derivative()/(2*y) for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
RR.<u,v> = QQ[]
F32 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F32[i] += FF[i][j]*v/(u^(N-2-j)) # be careful with the first coordinate: u --> x + A2/2

################################################################################

# test
[(l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^i*x.derivative()/(2*y) - F31[i](t-A2/2+A1/2,y).derivative() - d31[i]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) - F32[i](t,y).derivative() - d32[i]*(1/t*x.derivative()/(2*y))).valuation() for i in range(rang)]
### [+Infinity, +Infinity, 0, 0, 0]

################################################################################

dinf320 = ((2*y)*(l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^2*x.derivative()/(2*y) - F31[2](t-A2/2+A1/2,y).derivative() - d31[2]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) - F32[2](t,y).derivative() - d32[2]*(1/t*x.derivative()/(2*y)))).list()[0]
dinf331 = ((2*y)*(l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^3*x.derivative()/(2*y) - F31[3](t-A2/2+A1/2,y).derivative() - d31[3]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) - F32[3](t,y).derivative() - d32[3]*(1/t*x.derivative()/(2*y)))).list()[1]
dinf330 = ((2*y)*(l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^3*x.derivative()/(2*y) - F31[3](t-A2/2+A1/2,y).derivative() - d31[3]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) - F32[3](t,y).derivative() - d32[3]*(1/t*x.derivative()/(2*y)))).list()[0] + dinf331*A2/2
dinf342 = ((2*y)*(l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^4*x.derivative()/(2*y) - F31[4](t-A2/2+A1/2,y).derivative() - d31[4]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) - F32[4](t,y).derivative() - d32[4]*(1/t*x.derivative()/(2*y)))).list()[2]
dinf341 = ((2*y)*(l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^4*x.derivative()/(2*y) - F31[4](t-A2/2+A1/2,y).derivative() - d31[4]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) - F32[4](t,y).derivative() - d32[4]*(1/t*x.derivative()/(2*y)))).list()[1] + dinf342*A2
dinf340 = ((2*y)*(l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^4*x.derivative()/(2*y) - F31[4](t-A2/2+A1/2,y).derivative() - d31[4]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) - F32[4](t,y).derivative() - d32[4]*(1/t*x.derivative()/(2*y)))).list()[0] + dinf341*A2/2 - dinf342*A2^2/4

################################################################################

# test
[l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^i*x.derivative()/(2*y) == F31[i](t-A2/2+A1/2,y).derivative() + d31[i]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) + F32[i](t,y).derivative() + d32[i]*(1/t*x.derivative()/(2*y)) for i in range(2)]
l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^2*x.derivative()/(2*y) == F31[2](t-A2/2+A1/2,y).derivative() + d31[2]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) + F32[2](t,y).derivative() + d32[2]*(1/t*x.derivative()/(2*y)) + dinf320*x.derivative()/(2*y)
l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^3*x.derivative()/(2*y) == F31[3](t-A2/2+A1/2,y).derivative() + d31[3]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) + F32[3](t,y).derivative() + d32[3]*(1/t*x.derivative()/(2*y)) + dinf330*x.derivative()/(2*y) + dinf331*x*x.derivative()/(2*y)
l.truncate()(1/(t-A2/2+A1/2),1/t)*(t-A2/2)^4*x.derivative()/(2*y) == F31[4](t-A2/2+A1/2,y).derivative() + d31[4]*(1/(t-A2/2+A1/2)*x.derivative()/(2*y)) + F32[4](t,y).derivative() + d32[4]*(1/t*x.derivative()/(2*y)) + dinf340*x.derivative()/(2*y) + dinf341*x*x.derivative()/(2*y) + dinf342*x^2*x.derivative()/(2*y)
### [True, True]
### True
### True
### True

################################################################################

# RESULT:
# For i = 0,1, omega_i = dF31[i] + d31[i]*1/(x+A1/2)*dx/2y + dF32[i] + d32[i]*1/(x+A2/2)*dx/2y
# omega_2 = dF31[2] + d31[2]*1/(x+A1/2)*dx/2y + dF32[2] + d32[2]*1/(x+A2/2)*dx/2y + dinf320*dx/2y
# omega_3 = dF31[3] + d31[3]*1/(x+A1/2)*dx/2y + dF32[3] + d32[3]*1/(x+A2/2)*dx/2y + dinf330*dx/2y + dinf331*x*dx/2y
# omega_4 = dF31[4] + d31[4]*1/(x+A1/2)*dx/2y + dF32[4] + d32[4]*1/(x+A2/2)*dx/2y + dinf340*dx/2y + dinf341*x*dx/2y + dinf342*x^2*dx/2y

################################################################################

# INTEGRATION

def Int_vplus(i,z1,z2): # integral of omega_i on v_plus from z1 to z2
    x1 = L(z1[0])
    x2 = L(z2[0])
    exact_part1 = F1[i](x2 + A1/2) - F1[i](x1 + A1/2)
    exact_part2 = F2[i](x2 + A2/2) - F2[i](x1 + A2/2)
    exact_part3 = F3[i](x2 + A3/2) - F3[i](x1 + A3/2)
    exact_part = exact_part1 + exact_part2 + exact_part3
    third_kind_part1 = d1[i]*(Log(x2 + A1/2) - Log(x1 + A1/2))
    third_kind_part2 = d2[i]*(Log(x2 + A2/2) - Log(x1 + A2/2))
    third_kind_part3 = d3[i]*(Log(x2 + A3/2) - Log(x1 + A3/2))
    third_kind_part = third_kind_part1 + third_kind_part2 + third_kind_part3
    if i == 3:
        inf_part = dinf30*(x2 - x1)
    elif i == 4:
        inf_part = dinf40*(x2 - x1) + dinf41/2*(x2^2 - x1^2)
    else:
        inf_part = 0
    return exact_part + third_kind_part + inf_part

def Int_vminus(i,z1,z2): # integral of omega_i on v_minus from z1 to z2
    return -Int_vplus(i,z1,z2)

def NC(z,j): # points on v_j in the new coordinates
    x, y = L(z[0]), L(z[1])
    Z1, Z2, Z3 = x + A1/2, x + A2/2, x + A3/2
    l1 = Z1*(1 + C1/Z1^2).sqrt()
    l2 = Z2*(1 + C2/Z2^2).sqrt()
    l3 = Z3*(1 + C3/Z3^2).sqrt()
    if j == 1:
        l = l2*l3
    elif j == 2:
        l = l1*l3
    elif j == 3:
        l = l1*l2
    return (x,y/l)

def Int(i,j,z1,z2): # integral of omega_i on v_j from z1 to z2
    RR.<T> = K[]
    x1, y1 = NC(z1,j)
    x2, y2 = NC(z2,j)
    if j == 1:
        exact_part2 = F12[i](x2+A2/2,y2) - F12[i](x1+A2/2,y1)
        exact_part3 = F13[i](x2+A3/2,y2) - F13[i](x1+A3/2,y1)
        exact_part = exact_part2 + exact_part3
        pol12 = T^2 + (A2-A1)*T - C1
        pol13 = T^2 + (A3-A1)*T - C1
        r12, s12 = pol12.roots(multiplicities=False)
        r13, s13 = pol13.roots(multiplicities=False)
        T1 = x1 + y1 + A1/2
        T2 = x2 + y2 + A1/2
        third_kind_part_2 = d12[i]/(r12-s12)*(Log((T2-r12)/(T2-s12)) - Log((T1-r12)/(T1-s12)))
        third_kind_part_3 = d13[i]/(r13-s13)*(Log((T2-r13)/(T2-s13)) - Log((T1-r13)/(T1-s13)))
        third_kind_part = third_kind_part_2 + third_kind_part_3
        if i == 2:
            inf_part0 = dinf120/2*(Log(T2) - Log(T1))
            inf_part1 = 0
            inf_part2 = 0
        elif i == 3:
            inf_part0 = dinf130/2*(Log(T2) - Log(T1))
            inf_part1 = dinf131/4*((T2-A1*Log(T2)+C1/T2) - (T1-A1*Log(T1)+C1/T1))
            inf_part2 = 0
        elif i == 4:
            inf_part0 = dinf140/2*(Log(T2) - Log(T1))
            inf_part1 = dinf141/4*((T2-A1*Log(T2)+C1/T2) - (T1-A1*Log(T1)+C1/T1))
            inf_part2 = dinf142/8*((T2^2/2-2*A1*T2+(A1^2-2*C1)*Log(T2)-2*A1*C1/T2-C1^2/(2*T2^2)) - (T1^2/2-2*A1*T1+(A1^2-2*C1)*Log(T1)-2*A1*C1/T1-C1^2/(2*T1^2)))
        else:
            inf_part0 = 0
            inf_part1 = 0
            inf_part2 = 0
        inf_part = inf_part0 + inf_part1 + inf_part2
    elif j == 2:
        exact_part1 = F21[i](x2+A1/2,y2) - F21[i](x1+A1/2,y1)
        exact_part3 = F23[i](x2+A3/2,y2) - F23[i](x1+A3/2,y1)
        exact_part = exact_part1 + exact_part3
        pol21 = T^2 + (A1-A2)*T - C2
        pol23 = T^2 + (A3-A2)*T - C2
        r21, s21 = pol21.roots(multiplicities=False)
        r23, s23 = pol23.roots(multiplicities=False)
        T1 = x1 + y1 + A2/2
        T2 = x2 + y2 + A2/2
        third_kind_part_1 = d21[i]/(r21-s21)*(Log((T2-r21)/(T2-s21)) - Log((T1-r21)/(T1-s21)))
        third_kind_part_3 = d23[i]/(r23-s23)*(Log((T2-r23)/(T2-s23)) - Log((T1-r23)/(T1-s23)))
        third_kind_part = third_kind_part_1 + third_kind_part_3
        if i == 2:
            inf_part0 = dinf220/2*(Log(T2) - Log(T1))
            inf_part1 = 0
            inf_part2 = 0
        elif i == 3:
            inf_part0 = dinf230/2*(Log(T2) - Log(T1))
            inf_part1 = dinf231/4*((T2-A2*Log(T2)+C2/T2) - (T1-A2*Log(T1)+C2/T1))
            inf_part2 = 0
        elif i == 4:
            inf_part0 = dinf240/2*(Log(T2) - Log(T1))
            inf_part1 = dinf241/4*((T2-A2*Log(T2)+C2/T2) - (T1-A2*Log(T1)+C2/T1))
            inf_part2 = dinf242/8*((T2^2/2-2*A2*T2+(A2^2-2*C2)*Log(T2)-2*A2*C2/T2-C2^2/(2*T2^2)) - (T1^2/2-2*A2*T1+(A2^2-2*C2)*Log(T1)-2*A2*C2/T1-C2^2/(2*T1^2)))
        else:
            inf_part0 = 0
            inf_part1 = 0
            inf_part2 = 0
        inf_part = inf_part0 + inf_part1 + inf_part2
    elif j == 3:
        exact_part1 = F31[i](x2+A1/2,y2) - F31[i](x1+A1/2,y1)
        exact_part2 = F32[i](x2+A2/2,y2) - F32[i](x1+A2/2,y1)
        exact_part = exact_part1 + exact_part2
        pol31 = T^2 + (A1-A3)*T - C3
        pol32 = T^2 + (A2-A3)*T - C3
        r31, s31 = pol31.roots(multiplicities=False)
        r32, s32 = pol32.roots(multiplicities=False)
        T1 = x1 + y1 + A3/2
        T2 = x2 + y2 + A3/2
        third_kind_part_1 = d31[i]/(r31-s31)*(Log((T2-r31)/(T2-s31)) - Log((T1-r31)/(T1-s31)))
        third_kind_part_2 = d32[i]/(r32-s32)*(Log((T2-r32)/(T2-s32)) - Log((T1-r32)/(T1-s32)))
        third_kind_part = third_kind_part_1 + third_kind_part_2
        if i == 2:
            inf_part0 = dinf320/2*(Log(T2) - Log(T1))
            inf_part1 = 0
            inf_part2 = 0
        elif i == 3:
            inf_part0 = dinf330/2*(Log(T2) - Log(T1))
            inf_part1 = dinf331/4*((T2-A3*Log(T2)+C3/T2) - (T1-A3*Log(T1)+C3/T1))
            inf_part2 = 0
        elif i == 4:
            inf_part0 = dinf340/2*(Log(T2) - Log(T1))
            inf_part1 = dinf341/4*((T2-A3*Log(T2)+C3/T2) - (T1-A3*Log(T1)+C3/T1))
            inf_part2 = dinf342/8*((T2^2/2-2*A3*T2+(A3^2-2*C3)*Log(T2)-2*A3*C3/T2-C3^2/(2*T2^2)) - (T1^2/2-2*A3*T1+(A3^2-2*C3)*Log(T1)-2*A3*C3/T1-C3^2/(2*T1^2)))
        else:
            inf_part0 = 0
            inf_part1 = 0
            inf_part2 = 0
        inf_part = inf_part0 + inf_part1 + inf_part2
    return exact_part + third_kind_part + inf_part

################################################################################

# period integrals
per1 = [Int(i,1,P1,P2) + Int_vplus(i,P2,P3) + Int(i,2,P3,P4) + Int_vminus(i,P4,P1) for i in range(rang)]
per2 = [Int(i,2,P3,P4) + Int_vminus(i,P4,P5) + Int(i,3,P5,P6) + Int_vplus(i,P6,P3) for i in range(rang)]

################################################################################

path = [Int_vminus(i,S,P1) + Int(i,1,P1,P2) + Int_vplus(i,P2,R) for i in range(rang)] # image of path under tau = e1e2
for i in range(5):
    path[i] - 2/3*per1[i] + 1/3*per2[i]
### 2*a^34 + a^35 + a^38 + a^39 + 4*a^41 + 4*a^43 + O(a^44)
### 2*a^34 + a^35 + a^37 + a^38 + 3*a^39 + 3*a^40 + a^42 + a^43 + O(a^44)
### a^34 + 3*a^35 + 2*a^36 + 2*a^38 + a^40 + 3*a^41 + 4*a^42 + a^43 + O(a^44)
### 1 + 3*a^4 + a^8 + 3*a^12 + a^16 + 3*a^20 + a^24 + 3*a^28 + a^32 + 3*a^34 + 4*a^35 + 3*a^36 + a^37 + 4*a^38 + a^39 + 4*a^40 + 2*a^41 + a^42 + 4*a^43 + O(a^44)
### a^34 + 2*a^35 + 3*a^36 + 2*a^38 + 4*a^39 + a^40 + a^41 + 2*a^42 + a^43 + O(a^44)

################################################################################

path = [Int_vminus(i,S,P4) + Int(i,2,P4,P3) + Int_vplus(i,P3,R) for i in range(rang)] # image of path under tau = (-e4)(-e3)
for i in range(5):
    path[i] + 1/3*per1[i] + 1/3*per2[i]
### 2*a^34 + a^35 + a^38 + a^39 + 4*a^41 + 4*a^43 + O(a^44)
### 2*a^34 + a^35 + a^37 + a^38 + 3*a^39 + 3*a^40 + a^42 + a^43 + O(a^44)
### a^34 + 3*a^35 + 2*a^36 + 2*a^38 + a^40 + 3*a^41 + 4*a^42 + a^43 + O(a^44)
### 1 + 3*a^4 + a^8 + 3*a^12 + a^16 + 3*a^20 + a^24 + 3*a^28 + a^32 + 3*a^34 + 4*a^35 + 3*a^36 + a^37 + 4*a^38 + a^39 + 4*a^40 + 2*a^41 + a^42 + 4*a^43 + O(a^44)
### a^34 + 2*a^35 + 3*a^36 + 2*a^38 + 4*a^39 + a^40 + a^41 + 2*a^42 + a^43 + O(a^44)

################################################################################

path = [Int_vminus(i,S,P5) + Int(i,3,P5,P6) + Int_vplus(i,P6,R) for i in range(rang)] # image of path under tau = e5e6
for i in range(5):
    path[i] + 1/3*per1[i] - 2/3*per2[i]
### 2*a^34 + a^35 + a^38 + a^39 + 4*a^41 + 4*a^43 + O(a^44)
### 2*a^34 + a^35 + a^37 + a^38 + 3*a^39 + 3*a^40 + a^42 + a^43 + O(a^44)
### a^34 + 3*a^35 + 2*a^36 + 2*a^38 + a^40 + 3*a^41 + 4*a^42 + a^43 + O(a^44)
### 1 + 3*a^4 + a^8 + 3*a^12 + a^16 + 3*a^20 + a^24 + 3*a^28 + a^32 + 3*a^34 + 4*a^35 + 3*a^36 + a^37 + 4*a^38 + a^39 + 4*a^40 + 2*a^41 + a^42 + 4*a^43 + O(a^44)
### a^34 + 2*a^35 + 3*a^36 + 2*a^38 + 4*a^39 + a^40 + a^41 + 2*a^42 + a^43 + O(a^44)

################################################################################
