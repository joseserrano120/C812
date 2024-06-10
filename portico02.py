from fenics import *
mesh = Mesh("portico02.xml")
plot(mesh)

delta = 0.1
total_time = 60.0
k = 3.74

cara = MeshFunction(`size_t`, mesh, "portico02_facet_region.xml")

V = FunctionSpace(mesh, "P", 1)

ub = 360
bc1 = DirichletBC(V, Constant(ub), caras, 1715)
bc2 = DirichletBC(V, Constant(ub), caras, 1716)
bc = [bc1, bc2]

u = TrialFunction(V)
v = TestFunction(V)
un = Function(V)

Edificio = u*v*dx + deltat * k * dot( grad(u), grad(v))*dx - un * v * dx

a = lhs(Edificio)
l = rhs(Edificio)

vtkfile = File("solucion_fenics/edificio.pvd")

u = Function(V)
t = 0

num_iter = int(total_time/deltat)
for n in range(num_iter):
    t += deltat
    solve(a == l, u, bc)
    un.assing(u)
    if (n%5 == 0):
        print(n)
        vtkfile << (u, t)


