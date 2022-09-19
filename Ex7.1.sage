################################################################################

X = EllipticCurve([0,0,0, -1351755, 555015942]) # MW-group is Z x Z/6Z
P = X(330483/361,63148032/6859)
Q = X(2523,114912) # (Q)-(-Q) = P
R = X(219,16416) # (R)-(-R) = P
free = X(-1293,11880) # generator of the free part
tors = X(-501,33264) # generator of the torsion part
p = 43 # X has split multiplicative reduction at p
K = Qp(p,12)
XK = X.change_ring(K)
def Log(z): return K(z).log(p_branch=0, change_frac=True) # Iwasawa branch

################################################################################

RR.<x> = QQ[]
g = x - 507
gder = g.derivative()
N = 20 # truncation level
TT.<t> = PowerSeriesRing(K, 't', default_prec = 4*N)
l = (1-4635873/4*t)^(-1/2) # l starts with 1
UU = 0
for i in range(N):
    UU += l.list()[i]*((t^2)^(N-1-i))
U = UU/t^(2*N-1) # be careful with t: t --> x + 507/2

################################################################################

# 1) omega_i's
x = t - 507/2 # x-local-coordinate at the pole
y = g(x).sqrt() # y-local-coordinate at the pole
w =[U*x^i*x.derivative()/(2*y) for i in range(2)] # original forms around the pole
d = [w[i].residue()/(1/t*x.derivative()/(2*y)).residue() for i in range(2)] # coefficients of 1/(x+507/2)*dx/2y
w = [w[i] - d[i]/t*x.derivative()/(2*y) for i in range(2)] # original forms without residues
ED = [0,0] + [(t*gder(x)-2*i*g(x))/(t^(i+1))*x.derivative()/(2*y) for i in range(1,2*N-1)] # exact differentials
FF = [[] for i in range(2)]
for i in range(2):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
RR.<u,v> = QQ[]
F1 = [0 for i in range(2)] # meromorphic functions
for i in range(2):
    for j in range(2*N-2):
        F1[i] += FF[i][j]*v/(u^(2*N-2-j)) # be careful with the first coordinate: u --> x + 507/2
F1[1] = F1[1] + v

################################################################################

# TEST 1)
[U*x^i*x.derivative()/(2*y) == F1[i](t,y).derivative() + d[i]/t*x.derivative()/(2*y) for i in range(2)]
### [True, True]

# RESULT : omega_i = dF1[i] + d[i]/(x+507/2)*dx/2y

################################################################################

# 2) third kinds
# Let P = (a,b) and omega = b/(x-a)*dx/y. Then omega = dF + d0/(x-a)*dx/2y + d1/(x+507/2)*dx/2y.
def base(P):
    a, b = P[0], P[1]
    x = t - 507/2 # x-local-coordinate at the pole
    y = g(x).sqrt() # y-local-coordinate at the pole
    w = (b/(x-a))*U*x.derivative()/y # original form around the pole
    d1 = w.residue()/(1/t*x.derivative()/(2*y)).residue() # the coefficient of 1/(x+507/2)*dx/2y
    w = w - d1*(1/t*x.derivative()/(2*y)) # original form without residue
    ED = [0,0] + [(t*gder(x)-2*i*g(x))/(t^(i+1))*x.derivative()/(2*y) for i in range(1,2*N-1)] # exact differentials
    FF = []
    while w.valuation() < 0:
        val = w.valuation()
        FF.append(w.list()[0]/ED[-val].list()[0])
        w = w - (w.list()[0]/ED[-val].list()[0])*ED[-val]
    RR.<u,v> = QQ[]
    F = 0 # meromorphic function
    for i in range(2*N-2):
        F += FF[i]*v/(u^(2*N-2-i)) # be careful with the first coordinate: u--> x + 507/2
    d0 = (((b/(x-a))*U*x.derivative()/y-F(t,y).derivative()-d1/t*x.derivative()/(2*y))/(1/(x-a)*x.derivative()/(2*y))).list()[0] # the coefficient of 1/(x-a)*dx/2y
    return [F,d0,d1]

################################################################################

# TEST 2)
x = t - 507/2
y = g(x).sqrt()
for point in [P,Q,R,tors]:
    F, d0, d1 = base(point)
    a, b = point[0], point[1]
    (b/(x-a))*U*x.derivative()/y == F(t,y).derivative() + d0/(x-a)*x.derivative()/(2*y) + d1/t*x.derivative()/(2*y)
### True
### True
### True
### True

# RESULT : If P=(a,b), then base(P)=[F,d0,d1] with omega = b/(x-a)*dx/y = dF + d0/(x-a)*dx/2y + d1/(x+507/2)*dx/2y

################################################################################

def NC(z): # points in the new coordinates
    Z = z[0] + 507/2
    l = Z*(1-4635873/(4*Z^2)).sqrt()
    return (z[0],z[1]/l)
A = K(-3/2*507).sqrt()
def Int1(i,z1,z2): # integral of omega_i from z1 to z2
    z1 = XK(z1[0],z1[1])
    z2 = XK(z2[0],z2[1])
    exact_part = F1[i](NC(z2)[0]+507/2,NC(z2)[1]) - F1[i](NC(z1)[0]+507/2,NC(z1)[1])
    third_kind_part = d[i]/(2*A)*(Log((NC(z2)[1]-A)/(NC(z2)[1]+A)) - Log((NC(z1)[1]-A)/(NC(z1)[1]+A)))
    return exact_part + third_kind_part
def Int2(P,z1,z2): # integral of omega from z1 to z2, where omega = b/(x-a)*dx/y, P = (a,b)
    z1 = XK(z1[0],z1[1])
    z2 = XK(z2[0],z2[1])
    F, d0, d1 = base(P)
    exact_part = F(NC(z2)[0]+507/2,NC(z2)[1]) - F(NC(z1)[0]+507/2,NC(z1)[1])
    a = P[0]
    B = K(a-507).sqrt()
    third_kind_part0 = d0/(2*B)*(Log((NC(z2)[1]-B)/(NC(z2)[1]+B)) - Log((NC(z1)[1]-B)/(NC(z1)[1]+B)))
    third_kind_part1 = d1/(2*A)*(Log((NC(z2)[1]-A)/(NC(z2)[1]+A)) - Log((NC(z1)[1]-A)/(NC(z1)[1]+A)))
    return exact_part + third_kind_part0 + third_kind_part1
# D := (Q)-(-Q), w := y(Q)/(x-x(Q))*dx/y ===> w is a form of the third kind such that Res(w) = D
# psi(w) = c0*omega_0 + c1*omega_1, wD := w - c0*omega_0 ===> wD is the uniqe form of the third kind such that
# Res(wD) = D and psi(wD) is in W = <omega_1>. Easy computation by hand shows that c0 = integral of omega_1 from Q to -Q
def h_p(Q,R,S): # h_p((Q)-(-Q),(R)-(S)), Q,R,S are non-W'ss points
    c0 = Int1(1,Q,-Q)
    return Int2(Q,S,R) - c0*Int1(0,S,R)

################################################################################

# BEFORE MORE TESTS, we have to be sure that the points that we will consider lie in the component U_1
LLL = []

for i in range(6):
    for j in range(1,6):
        for k in range(1,6):
            LLL.append(K((i*P + j*tors)[0]-26).abs())
            LLL.append(K((i*P + k*tors)[0]-26).abs())
            LLL.append(K((i*R + j*tors)[0]-26).abs())
            LLL.append(K((i*R + k*tors)[0]-26).abs())
            LLL.append(K((i*Q + j*tors)[0]-26).abs())
            LLL.append(K((i*Q + k*tors)[0]-26).abs()) # for test 1

for i in range(1,10):
    for j in range(1,10):
        LLL.append(K((i*P)[0]-26).abs())
        LLL.append(K((i*R)[0]-26).abs())
        LLL.append(K((j*R)[0]-26).abs())
        LLL.append(K((j*Q)[0]-26).abs())
        LLL.append(K((i*P)[0]-26).abs())
        LLL.append(K((i*P)[0]-26).abs()) # for test 2

for i in range(1,6):
    for j in range(6,11):
        LLL.append(K((i*R)[0]-26).abs())
        LLL.append(K((j*Q)[0]-26).abs())
        LLL.append(K((i*-R)[0]-26).abs())
        LLL.append(K((j*-Q)[0]-26).abs()) # for test 3

for point in [2*tors,-tors,R+2*tors,R-2*tors,-2*tors,Q+2*tors,Q-5*tors]:
    LLL.append(K(point[0]-26).abs()) # for test 4

Set(LLL)
### {1}

################################################################################

# TEST 3) vanishing of omega_0 between points with torsion difference
LLL = []
for i in range(6):
    for j in range(1,6):
        for k in range(1,6):
            LLL.append(Int1(0,i*P + j*tors,i*P + k*tors))
            LLL.append(Int1(0,i*R + j*tors,i*R + k*tors))
            LLL.append(Int1(0,i*Q + j*tors,i*Q + k*tors))

Set(LLL)
### {O(43^12)}

################################################################################

# TEST 4) formal log implementation
def ForLog(z1,z2):
    formal_log = XK.formal_group().log()
    c = X.tamagawa_number(p)
    Z2 = c*(p-1)*z2
    Z1 = c*(p-1)*z1
    log2 = 1/(c*(p-1))*formal_log.truncate()(-Z2[0]/Z2[1])
    log1 = 1/(c*(p-1))*formal_log.truncate()(-Z1[0]/Z1[1])
    return log2 - log1
LLL = []
for i in range(1,10):
    for j in range(1,10):
        LLL.append(Int1(0,i*P,j*R) == ForLog(i*P,j*R))
        LLL.append(Int1(0,i*P,j*Q) == ForLog(i*P,j*Q))
        LLL.append(Int1(0,i*R,j*Q) == ForLog(i*R,j*Q))

Set(LLL)
### {True}

################################################################################

# TEST 5) symmetry of h_p (when S = -R)
LLL = []
for i in range(1,5):
    for j in range(6,10):
        LLL.append(h_p(j*Q,i*R,i*(-R)) == h_p(i*R,j*Q,j*(-Q)))

Set(LLL)
### {True}

################################################################################

# TEST 6) vanishing of Coleman-Gross height
h_p(Q,2*tors,-tors) + 2*Log(2) - Log(3) - Log(5)
h_p(Q,R+2*tors,R-2*tors) - Log(53) + Log(67)
h_p(R,tors,-2*tors) - 2*Log(2) + Log(3) + Log(5)
h_p(R,Q+2*tors,Q-5*tors) - Log(2) - Log(3) - Log(7) + Log(67)
### O(43^12)
### O(43^12)
### O(43^12)
### O(43^12)

################################################################################

# TEST 7) quadraticty of Coleman-Gross height
# Since 2R=2Q, we can't simply use h_p(2Q,2R,-2R) in order to compute p-part of Coleman-Gross height of 2*P;
# because in this case poles of the third kind form and endpoints coincide. We will represent 2P differently.
Rnew = 2*tors
Snew = 4*free
2*P == 2*Q-(-2*Q) == Rnew - Snew
### True

################################################################################

height_of_P = h_p(Q,R,-R) + 9*Log(2) # Coleman-Gross height of P
height_of_2P = h_p(2*Q,Rnew,Snew) - 2*Log(3) + Log(5) - 2*Log(7) + Log(17) + Log(89) + Log(4139) + Log(8612423) # Coleman-Gross height of 2*P
height_of_2P == 4*(height_of_P)
### True

