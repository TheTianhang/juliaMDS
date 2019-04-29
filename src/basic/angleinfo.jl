struct angle
    tag_a::Int64
    tag_b::Int64
    tag_c::Int64
end

struct anglebasic
    angles::Vector{angle}
    tag_angle::Vector{Int64}
    idx_angle::Vector{Int64}
end

function get_Nangles(anglebasic::anglebasic)
    return length(bondinfo.bonds)
end

function get_NangleTypes(anglebasic::anglebasic)
    angletype=[]
    for angleii in 1:length(anglebasic.tag_angle)
        if anglebasic.tag_angled[angleii] in angletype
            continue
        else
            push!(angletype,anglebasic.tag_angle[angleii])
        end
    end
    return length(angletype) 
end

struct angleinfo
    angles::Vector{angle}
    tag_angle::Vector{Int64}
    idx_angle::Vector{Int64}
    Nangles::Int64
    NangleTypes::Int64
end
    
function get_angleinfo(anglebasic::anglebasic)
    Nangles=get_Nangles(anglebasic)
    NangleTypes=get_NangleTypes(anglebasic)
    return bondinfo(anglebasic.angles,anglebasic.tag_angle,anglebasic.idx_angle,Nangles,NangleTypes)
end
