using GeometryTypes
using PyPlot


res = 0.1

s = SignedDistanceField(HyperRectangle(Vec(0,0.),Vec(pi*1,pi*1)), res) do v
    implicit_triangle(3,3,3,v[1],v[2],pi/3)
end


"""
Find the value closest to zero and return the linear index of the element.
"""
function find_zeros{T}(mat::Array{T})
    x = typemax(T)
    ind = 0
    # find the value closest to zero
    @inbounds for i = 1:length(mat)
        val = abs(mat[i])
        if val < x
            ind = i
            x = val
        end
    end
    ind
end

i = find_zeros(s.data)

println("the value closest to zero occurs at index $(ind2sub(s.data, i)) with value $(s.data[i])")

# grab indices
x = Vec(ind2sub(s.data, i)...) # make a vector for now so we can find acutal vals
x = s.bounds.minimum + x.* res # scale for resolution

@show x

surf(s.data)

