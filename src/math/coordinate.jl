function vector3d(coordinate_a::Coordinate,coordinate_b::Coordinate)
    vector=zeros(3)
    vector[1]=coordinate_a.x-coordinate_b.x
    vector[2]=coordinate_a.y-coordinate_b.y
    vector[3]=coordinate_a.z-coordinate_b.z
    return vector
end

function dot(coordinate_a::Coordinate,coordinate_b::Coordinate)
    dx=coordinate_a.x*coordinate_b.x
    dy=coordinate_a.y*coordinate_b.y
    dz=coordinate_a.z*coordinate_b.z
    return dx+dy+dz
end

function sqr(coordinate_a::Coordinate,coordinate_b::Coordinate)
    return coordinate_a.x*coordinate_b.x+coordinate_a.y*coordinate_b.y+coordinate_a.z*coordinate_b.z
end

function sqr(vector::Vector{Float64})
    result=0
    for i in 1:length(a)
        result+=vector[i]^2
    end
    return result
end