using ParametricPolyhedra
using Base.Test

# equilateral triangle
a = ParametricTriangle(1.0,1.0,1.0,pi/3,pi/3,pi/3)
@test typeof(a) == ParametricTriangle{Float64}
@test !isnull(a)
@test isvalid(a,atol=eps(1.0))

# right triangle
a = ParametricTriangle(sqrt(2),1.0,1.0,pi/2,pi/4,pi/4)
@test typeof(a) == ParametricTriangle{Float64}
@test !isnull(a)
@test isvalid(a,atol=eps(1.0))

# Not valid (angles == 5pi/4)
a = ParametricTriangle(sqrt(2),1.0,1.0,pi/2,pi/4,pi/2)
@test !isvalid(a,atol=eps(1.0))

e1 = edges(pi/3,pi/3,pi/3)
@test e1 == (sin(pi/3), sin(pi/3), sin(pi/3))
e1 = edges(pi/3,pi/3,pi/3,2)
@test e1 == (2*sin(pi/3), 2*sin(pi/3), 2*sin(pi/3))
e1 = edges(3*pi/5,pi/5,pi/5)
@test e1 == (0.9510565162951536,0.5877852522924731,0.5877852522924731)

# missing value, fill tests
a = ParametricTriangle(Nullable{Float64}(),1.0,1.0,pi/2,pi/4,pi/4)
@test isnull(a)
af = fill(a)
@test isvalid(af,atol=eps(10.0))
@test isapprox(get(af.a),sqrt(2))

a = ParametricTriangle(Nullable{Float64}(),Nullable{Float64}(),Nullable{Float64}(),
                       pi/2,pi/4,pi/4)
@test isnull(a)
af = fill(a)
@test isvalid(af,atol=eps(10.0))
