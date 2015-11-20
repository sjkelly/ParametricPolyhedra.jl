module ParametricPolyhedra

export ParametricTriangle, edges

type ParametricTriangle{T}
    # edge lengths
    a::Nullable{T}
    b::Nullable{T}
    c::Nullable{T}
    # angles (radians)
    alpha::Nullable{T}
    beta::Nullable{T}
    gamma::Nullable{T}
end

# functions to allow mixed constructors
makenullable(a::Nullable) = a
makenullable(a) = Nullable(a)

# mixed case
function ParametricTriangle(a, b, c,
                            alpha, beta, gamma)
    ParametricTriangle(makenullable(a),
                       makenullable(b),
                       makenullable(c),
                       makenullable(alpha),
                       makenullable(beta),
                       makenullable(gamma))
end

function Base.isnull(p::ParametricTriangle)
    isnull(p.a) || isnull(p.b) || isnull(p.c) ||
    isnull(p.alpha) || isnull(p.beta) || isnull(p.gamma)
end

"""
Test if a ParametricTriangle has a valid configuration.
"""
function Base.isvalid(p::ParametricTriangle)
    # underdetermined case
    isnull(p) && return false
    # otherwise check constraints since all values exist
    a = get(p.a)
    b = get(p.b)
    c = get(p.c)
    alpha = get(p.alpha)
    beta = get(p.beta)
    gamma = get(p.gamma)
    return a*cos(beta) + b*cos(alpha) - c == 0 &&
           b*sin(alpha) - a*sin(beta) == 0 &&
           alpha + beta + gamma - pi == 0
end

# version with isapprox for floats
function Base.isvalid{T<:AbstractFloat}(p::ParametricTriangle{T};
                                        rtol=sqrt(eps(T)),
                                        atol=zero(T))
    # underdetermined case
    isnull(p) && return false
    # otherwise check constraints since all values exist
    a = p.a.value
    b = p.b.value
    c = p.c.value
    alpha = p.alpha.value
    beta = p.beta.value
    gamma = p.gamma.value
    return isapprox(a*cos(beta) + b*cos(alpha) - c,0,
                    rtol=rtol,atol=atol) &&
           isapprox(b*sin(alpha) - a*sin(beta),0,
                    rtol=rtol,atol=atol) &&
           isapprox(alpha + beta + gamma - pi,0,
                    rtol=rtol,atol=atol)
end

"""
Given three angles, compute the edge lengths. The circumcircle radius
default is 1.0, but can be specified as a fourth argument.
"""
function edges{T}(alpha::T,beta::T,gamma::T,c=one(T))
    (c*sin(alpha),c*sin(beta),c*sin(gamma))
end

"""
Given angle `a` and circumcircle diameter, `c`,
compute the opposite edge length in the triangle.
"""
function _edge(a,c)
    c*sin(a)
end

"""
Given an underdetermined ParametricTriangle, compute the missing values
and return a new ParametricTriangle
"""
function Base.fill(p::ParametricTriangle)
    # all angles must be specified
    if isnull(p.alpha) || isnull(p.alpha) || isnull(p.gamma)
        error("Cannot fill in values for this triangle. All angles must be specified")
    end
    alpha = get(p.alpha)
    beta = get(p.beta)
    gamma = get(p.gamma)
    # no edges given, use circumcircle=1
    if isnull(p.a) && isnull(p.b) && isnull(p.c)
        e = edges(alpha,beta,gamma)
        return ParametricTriangle(e[1],e[2],e[3],p.alpha,p.beta,p.gamma)
    else
        # find the circumcircle
        D = !isnull(p.a) ? get(p.a)/sin(alpha) :
            !isnull(p.b) ? get(p.b)/sin(beta) :
            get(p.c)/sin(gamma) # one must be specified because of prior check
        # we only need to figure one side that is specified
        # so we can (re)compute the other two
        if !isnull(p.a)
            return ParametricTriangle(p.a, _edge(beta,D), _edge(gamma,D),
                                      p.alpha, p.beta, p.gamma)
        elseif !isnull(p.b)
            return ParametricTriangle(_edge(alpha,D), p.b, _edge(gamma,D),
                                      p.alpha, p.beta, p.gamma)
        elseif !isnull(p.c)
            return ParametricTriangle(_edge(alpha,D), _edge(beta,D), p.c,
                                      p.alpha, p.beta, p.gamma)
        end
    end
    # otherwise we broke it ... TODO remove once we start using this in loops
    error("filling triangle failed! What did you do???")
end

end # module
