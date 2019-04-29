const Reg_SingleInt64=r"\s+\d*\s+"
const Reg_SingleFloat=r"\s+\d+\.\d+\s*"
const Reg_DoubleFloat=r"\s+-?\d+\.\d+\s+-?\d+\.\d+"
const Reg_DoubleInt=r"\s+-?\d*\s+-?\d*"
const Reg_SinglePhrase=r"[a-z]+"
const Reg_DoublePhrases=r"[a-z]+\s[a-z]+"
const Reg_HeadCapitalPhrase=r"[A-Z][a-z]+"
const Reg_DoubleHeadCapitalPhrases=r"[A-Z][a-z]+\s[A-Z][a-z]+"

function fetch_NumberAtoms(lines::Array{String,1})::Int64
    for linecontent in lines
        Phrase=match(Reg_SinglePhrase,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="atoms"
                Data=match(Reg_SingleInt64,linecontent)
                return parse(Int,Data.match)
            end
        end
    end
end

function fetch_NumberBonds(lines::Array{String,1})::Int64
    for linecontent in lines
        Phrase=match(Reg_SinglePhrase,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="bonds"
                Data=match(Reg_SingleInt64,linecontent)
                return parse(Int,Data.match)
            end
        end
    end
    return 0
end

function fetch_NumberAngles(lines::Array{String,1})::Int64
    for linecontent in lines
        Phrase=match(Reg_SinglePhrase,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="angles"
                Data=match(Reg_SingleInt64,linecontent)
                return parse(Int,Data.match)
            end
        end
    end
    return 0
end

function fetch_NumberDihedrals(lines::Array{String,1})::Int64
    for linecontent in lines
        Phrase=match(Reg_SinglePhrase,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="dihedrals"
                Data=match(Reg_SingleInt64,linecontent)
                return parse(Int,Data.match)
            end
        end
    end
    return 0
end

function fetch_NumberImpropers(lines::Array{String,1})::Int64
    for linecontent in lines
        Phrase=match(Reg_SinglePhrase,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="impropers"
                Data=match(Reg_SingleInt64,linecontent)
                return parse(Int,Data.match)
            end
        end
    end
    return 0
end

function fetch_NumberAtomTypes(lines::Array{String,1})::Int64
    for linecontent in lines
        Phrase=match(Reg_DoublePhrases,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="atom types"
                Data=match(Reg_SingleInt64,linecontent)
                return parse(Int,Data.match)
            end
        end
    end
    return 0
end

function fetch_NumberBondTypes(lines::Array{String,1})::Int64
    for linecontent in lines
        Phrase=match(Reg_DoublePhrases,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="bond types"
                Data=match(Reg_SingleInt64,linecontent)
                return parse(Int,Data.match)
            end
        end
    end
    return 0
end

function fetch_NumberAngleTypes(lines::Array{String,1})::Int64
    for linecontent in lines
        Phrase=match(Reg_DoublePhrases,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="angle types"
                Data=match(Reg_SingleInt64,linecontent)
                return parse(Int,Data.match)
            end
        end
    end
    return 0
end

function fetch_NumberDihedralTypes(lines::Array{String,1})::Int64
    for linecontent in lines
        Phrase=match(Reg_DoublePhrases,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="dihedral types"
                Data=match(Reg_SingleInt64,linecontent)
                return parse(Int,Data.match)
            end
        end
    end
    return 0
end

function fetch_NumberImproperTypes(lines::Array{String,1})::Int64
    for linecontent in lines
        Phrase=match(Reg_DoublePhrases,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="improper types"
                Data=match(Reg_SingleInt64,linecontent)
                return parse(Int,Data.match)
            end
        end
    end
    return 0
end

function fetch_BoxSize_XAxis(lines::Array{String,1})
    for linecontent in lines
        Phrase=match(Reg_DoublePhrases,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="xlo xhi"
                data=match(Reg_DoubleFloat,linecontent)
                if data==nothing
                    data=match(Reg_DoubleInt,linecontent)
                end
                xlo,xhi=parse(Float64,split(data.match)[1]),parse(Float64,split(data.match)[2])
                
                return xlo,xhi
            end
        end
    end
    return 0
end

function fetch_BoxSize_YAxis(lines::Array{String,1})
    BoxSize=zeros(2)
    for linecontent in lines
        Phrase=match(Reg_DoublePhrases,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="ylo yhi"
                data=match(Reg_DoubleFloat,linecontent)
                if data==nothing
                    data=match(Reg_DoubleInt,linecontent)
                end
                ylo,yhi=parse(Float64,split(data.match)[1]),parse(Float64,split(data.match)[2])
                return ylo,yhi
            end
        end
    end
end

function fetch_BoxSize_ZAxis(lines::Array{String,1})
    BoxSize=zeros(2)
    for linecontent in lines
        Phrase=match(Reg_DoublePhrases,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="zlo zhi"
                data=match(Reg_DoubleFloat,linecontent)
                if data==nothing
                    data=match(Reg_DoubleInt,linecontent)
                end
                zlo,zhi=parse(Float64,split(data.match)[1]),parse(Float64,split(data.match)[2])
                return zlo,zhi
            end
        end
    end
end

function fetch_Masses(lines::Array{String,1},NumberAtomTypes::Int64)::Dict
    Masses=Dict()
    for (lineindex,linecontent) in enumerate(lines)
        Phrase=match(Reg_HeadCapitalPhrase,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="Masses"
                #println(Phrase.match)
                for atomtypeii in 1:NumberAtomTypes
                    #println(atomtypeii)
                    #println(lines[lineindex+1+atomtypeii])
                    data=parse(Float64,match(Reg_SingleFloat,lines[lineindex+1+atomtypeii]).match)
                    get!(Masses,atomtypeii,data)
                end
                return Masses
            else
                continue
            end
        end
    end
end

function fetch_BondConnectivity(lines::Array{String,1},NumberBonds::Int64)
    bonds=bond[]
    tag_bond=Int64[]
    for (lineindex,linecontent) in enumerate(lines)
        Phrase=match(Reg_HeadCapitalPhrase,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="Bonds"
                for bondii in 1:NumberBonds
                    bondindex,bondtype,atomindex1,atomindex2=split(lines[lineindex+1+bondii])
                    push!(bonds,bond(parse(Int,atomindex1),parse(Int,atomindex2)))
                    push!(tag_bond,parse(Int,bondtype))
                    #selectdim(BondConnectivity,1,bondii)[1]=parse(Int,bondindex)
                    #selectdim(BondConnectivity,1,bondii)[2]=parse(Int,bondtype)
                    #selectdim(BondConnectivity,1,bondii)[3]=parse(Int,atomindex1)
                    #selectdim(BondConnectivity,1,bondii)[4]=parse(Int,atomindex2)
                    #BondConnectivity[bondii,:,:][1][1]=parse(Int,bondindex)
                    #BondConnectivity[bondii,:,:][1][2]=parse(Int,bondtype)
                    #BondConnectivity[bondii,:,:][1][3]=parse(Int,atomindex1)
                    #BondConnectivity[bondii,:,:][1][4]=parse(Int,atomindex2)
                end
                return bonds,tag_bond
            else
                continue
            end
        end
        return bonds,tag_bond
    end
end

function fetch_AtomInfo(lines::Array{String,1},NumberAtoms::Int64,AtomStyle::String,masses::Dict)
    coordinates=coordinate[]
    tags_atom=tag_atom[]
    for (lineindex,linecontent) in enumerate(lines)
        Phrase=match(Reg_HeadCapitalPhrase,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="Atoms"
                if AtomStyle=="full"
                    for atomii in 1:NumberAtoms
                        atomindex,moleculeindex,atomtype,atomcharge,x,y,z=split(lines[lineindex+1+atomii])
                        #selectdim(AtomInfo,1,atomii)[1]=parse(Int,atomindex)
                        #selectdim(AtomInfo,1,atomii)[2]=parse(Int,moleculeindex)
                        #selectdim(AtomInfo,1,atomii)[4]=parse(Float64,atomcharge)
                        push!(tags_atom,tag_atom(parse(Int,atomtype)))
                        push!(coordinates,coordinate(parse(Float64,x),parse(Float64,y),parse(Float64,z),masses[parse(Int,atomtype)]))
                    end
                    return coordinates,tags_atom
                elseif AtomStyle=="bond"
                    for atomii in 1:NumberAtoms
                        atomindex,moleculeindex,atomtype,x,y,z=split(lines[lineindex+1+atomii])
                        #selectdim(AtomInfo,1,atomii)[1]=parse(Int,atomindex)
                        #selectdim(AtomInfo,1,atomii)[2]=parse(Int,moleculeindex)
                        #selectdim(AtomInfo,1,atomii)[3]=parse(Int,atomtype)
                        #selectdim(AtomInfo,1,atomii)[4]=parse(Float64,x)
                        #selectdim(AtomInfo,1,atomii)[5]=parse(Float64,y)
                        #selectdim(AtomInfo,1,atomii)[6]=parse(Float64,z)
                        push!(tags_atom,tag_atom(parse(Int,atomtype)))
                        push!(coordinates,coordinate(parse(Float64,x),parse(Float64,y),parse(Float64,z),masses[parse(Int,atomtype)]))
                    end
                    return coordinates,tags_atom
                else
                    error("unknow atom_style")
                end
            
            end
        end
    end
    return coordinates,tags_atom
end

function fetch_AngleConnectivity(lines::Array{String,1},NumberAngles::Int64)
    angles=angle[]
    tag_angle=Int64[]
    for (lineindex,linecontent) in enumerate(lines)
        Phrase=match(Reg_HeadCapitalPhrase,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="Angles"
                for angleii in 1:NumberAngles
                    angleindex,angletype,atomindex1,atomindex2,atomindex3=split(lines[lineindex+1+angleii])
                    push!(angles,angle(parse(Int,atomindex1),parse(Int,atomindex2),parse(Int,atomindex3)))
                    push!(tag_angle,parse(Int,angleindex))
                    #selectdim(AngleConnectivity,1,angleii)[1]=parse(Int,angleindex)
                    #selectdim(AngleConnectivity,1,angleii)[2]=parse(Int,angletype)
                    #selectdim(AngleConnectivity,1,angleii)[3]=parse(Int,atomindex1)
                    #selectdim(AngleConnectivity,1,angleii)[4]=parse(Int,atomindex2)
                    #selectdim(AngleConnectivity,1,angleii)[5]=parse(Int,atomindex3)
                end
                return angles,tag_angle
            else
                continue
            end
        end
        return angles,tag_angle
    end
end

function fetch_DihedralConnectivity(lines::Array{String,1},NumberDihedrals::Int64)::Array
    DihedralConnectivity=zeros(Int64,NumberDihedrals,6)
    for (lineindex,linecontent) in enumerate(lines)
        Phrase=match(Reg_HeadCapitalPhrase,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="Dihedrals"
                for dihedralii in 1:NumberDihedrals
                    dihedralindex,dihedraltype,atomindex1,atomindex2,atomindex3,atomindex4=split(lines[lineindex+1+dihedralii])
                    selectdim(DihedralConnectivity,1,dihedralii)[1]=parse(Int,dihedralindex)
                    selectdim(DihedralConnectivity,1,dihedralii)[2]=parse(Int,dihedraltype)
                    selectdim(DihedralConnectivity,1,dihedralii)[3]=parse(Int,atomindex1)
                    selectdim(DihedralConnectivity,1,dihedralii)[4]=parse(Int,atomindex2)
                    selectdim(DihedralConnectivity,1,dihedralii)[5]=parse(Int,atomindex3)
                    selectdim(DihedralConnectivity,1,dihedralii)[6]=parse(Int,atomindex4)
                end
                return DihedralConnectivity
            else
                continue
            end
        end
        return DihedralConnectivity
    end
end

function fetch_ImproperConnectivity(lines::Array{String,1},NumberImpropers::Int64)::Array
    ImproperConnectivity=zeros(Int64,NumberImpropers,6)
    for (lineindex,linecontent) in enumerate(lines)
        Phrase=match(Reg_HeadCapitalPhrase,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="Impropers"
                for improperii in 1:NumberImpropers
                    improperindex,impropertype,atomindex1,atomindex2,atomindex3,atomindex4=split(lines[lineindex+1+improperii])
                    selectdim(ImproperConnectivity,1,improperii)[1]=parse(Int,improperindex)
                    selectdim(ImproperConnectivity,1,improperii)[2]=parse(Int,impropertype)
                    selectdim(ImproperConnectivity,1,improperii)[3]=parse(Int,atomindex1)
                    selectdim(ImproperConnectivity,1,improperii)[4]=parse(Int,atomindex2)
                    selectdim(ImproperConnectivity,1,improperii)[5]=parse(Int,atomindex3)
                    selectdim(ImproperConnectivity,1,improperii)[6]=parse(Int,atomindex4)
                end
                return ImproperConnectivity
            else
                continue
            end
        end
        return ImproperConnectivity
    end
end

struct input
    coordinates::Vector{coordinate}
    #velocities::Vector{velocity}
    #images::Vector{image}
    tags_atom::Vector{tag_atom}
    bonds::Vector{bond}
    tag_bond::Vector{Int64}
    #ngles::Vector{angle}
    #ag_angle::Vector{Int64}
    boxsize::Vector{Float64}
    inittimestep::Int64
end


function Parse_Input(Filename::String;AtomStyle::String="bond")
    lines=readlines(Filename)
    
    NumberAtoms=fetch_NumberAtoms(lines)
    NumberBonds=fetch_NumberBonds(lines)
    NumberAngles=fetch_NumberAngles(lines)
    #println(NumberAngles)
    NumberDihedrals=fetch_NumberDihedrals(lines)
    NumberImpropers=fetch_NumberImpropers(lines)

    NumberAtomTypes=fetch_NumberAtomTypes(lines)
    NumberBondTypes=fetch_NumberBondTypes(lines)
    NumberAngleTypes=fetch_NumberAngleTypes(lines)
    NumberDihedralTypes=fetch_NumberDihedralTypes(lines)
    NumberImproperTypes=fetch_NumberImproperTypes(lines)

    xlo,xhi=fetch_BoxSize_XAxis(lines)
    ylo,yhi=fetch_BoxSize_YAxis(lines)
    zlo,zhi=fetch_BoxSize_ZAxis(lines)

    boxsize=[xlo,xhi,ylo,yhi,zlo,zhi]

    masses=fetch_Masses(lines,NumberAtomTypes)
    inittimestep=0
    coordinates,tags_atom=fetch_AtomInfo(lines,NumberAtoms,AtomStyle,masses)
    
    bonds,tag_bond=fetch_BondConnectivity(lines,NumberBonds)
    #ngles,tag_angle=fetch_AngleConnectivity(lines,NumberAngles)
    #ihedrals=fetch_DihedralConnectivity(lines,NumberDihedrals)
    #mproper=fetch_ImproperConnectivity(lines,NumberImpropers)

    return input(coordinates,tags_atom,bonds,tag_bond,boxsize,inittimestep)
    
end