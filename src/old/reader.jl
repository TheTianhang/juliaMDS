include("basic.jl")

module reader

function lammpsreader(data_file::AbstractString)
    lines=readlines(data_file)
    for l in eachline(data_file)
        sl = strip(l)
        if length(sl) == 0 || startswith(sl, ';')
            continue
        end
        ls=split(sl,r"\s+")
        if length(ls) > 1
            if ls[2]=="atoms"
                numberatoms=ls[1]
            elseif ls[2]=="bonds"
                numberbonds=ls[1]
            elseif ls[2]=="angles"
                numberangles=ls[1]
            elseif ls[2]=="dihedrals"
                numberdihedrals=ls[1]
            elseif ls[2]=="impropers"
                numberimpropers=ls[1]
            elseif ls[2]=="atom" & ls[3]=="types"
                numberatomtype=ls[1]
            elseif ls[2]=="bond" & ls[3]=="types"
                numberbondtype=ls[1]
            elseif ls[2]=="angle" & ls[3]=="types"
                numberangletype=ls[1]
            elseif ls[2]=="dihedral" & ls[3]=="types"
                numberdihedraltype=ls[1]
            elseif ls[2]=="improper" & ls[3]=="types"
                numberimpropertype=ls[1]
            elseif ls[3]=="xlo" & ls[4] == "xhi"
                xlo,xhi=ls[1:2]
            elseif ls[3]=="ylo" & ls[4] == "yhi"
                ylo,yhi=ls[1:2]
            elseif ls[3]=="zlo" & ls[4] == "zhi"
                zlo,zhi=ls[1:2]
            elseif ls[1]=="Masses"
                readMass(numberatomtype)
            elseif ls[1]=="Atoms"
                readAtoms(numberatoms)
            elseif ls[1]=="Bonds"
                readBonds(numberbonds)
            elseif ls[1]=="Angles"
                readAngles(numberangles)
            elseif ls[1]=="Dihedrals"
                readDihedrals(numberdihedrals)
            elseif ls[1]=="Impropers"
                readImpropers(numberimpropers)
            else
                println("ERROR for reading")
        end
    end
end

function readMass(numberatomtype::UInt64)

end

Reg_SingleInt64=r"\s+\d*\s+"
Reg_SingleFloat=r"\s+\d+\.\d+\s+"
Reg_DoubleFloat=r"\s+-?\d+\.\d+\s+-?\d+\.\d+"
Reg_SinglePhrase=r"[a-z]+"
Reg_DoublePhrases=r"[a-z]+\s[a-z]+"
Reg_HeadCapitalPhrase=r"[A-Z][a-z]+"
Reg_DoubleHeadCapitalPhrases=r"[A-Z][a-z]+\s[A-Z][a-z]+"

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
end

function fetch_BoxSize_XAxis(lines::Array{String,1})::Array{Float64}
    BoxSize=zeros(2)
    for linecontent in lines
        Phrase=match(Reg_DoublePhrases,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="xlo xhi"
                data=match(Reg_DoubleFloat,linecontent)
                BoxSize[1],BoxSize[2]=parse(Float64,split(data.match)[1]),parse(Float64,split(data.match)[2])
                return BoxSize
            end
        end
    end
end

function fetch_BoxSize_YAxis(lines::Array{String,1})::Array{Float64}
    BoxSize=zeros(2)
    for linecontent in lines
        Phrase=match(Reg_DoublePhrases,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="ylo yhi"
                data=match(Reg_DoubleFloat,linecontent)
                BoxSize[1],BoxSize[2]=parse(Float64,split(data.match)[1]),parse(Float64,split(data.match)[2])
                return BoxSize
            end
        end
    end
end

function fetch_BoxSize_ZAxis(lines::Array{String,1})::Array{Float64}
    BoxSize=zeros(2)
    for linecontent in lines
        Phrase=match(Reg_DoublePhrases,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="zlo zhi"
                data=match(Reg_DoubleFloat,linecontent)
                BoxSize[1],BoxSize[2]=parse(Float64,split(data.match)[1]),parse(Float64,split(data.match)[2])
                return BoxSize
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
                for atomtypeii in 1:NumberAtomTypes
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

function fetch_BondConnectivity(lines::Array{String,1},NumberBonds::Int64)::Array
    BondConnectivity=zeros(Int64,NumberBonds,4)
    for (lineindex,linecontent) in enumerate(lines)
        Phrase=match(Reg_HeadCapitalPhrase,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="Bonds"
                for bondii in 1:NumberBonds
                    bondindex,bondtype,atomindex1,atomindex2=split(lines[lineindex+1+bondii])
                    selectdim(BondConnectivity,1,bondii)[1]=parse(Int,bondindex)
                    selectdim(BondConnectivity,1,bondii)[2]=parse(Int,bondtype)
                    selectdim(BondConnectivity,1,bondii)[3]=parse(Int,atomindex1)
                    selectdim(BondConnectivity,1,bondii)[4]=parse(Int,atomindex2)
                    #BondConnectivity[bondii,:,:][1][1]=parse(Int,bondindex)
                    #BondConnectivity[bondii,:,:][1][2]=parse(Int,bondtype)
                    #BondConnectivity[bondii,:,:][1][3]=parse(Int,atomindex1)
                    #BondConnectivity[bondii,:,:][1][4]=parse(Int,atomindex2)
                end
                return BondConnectivity
            else
                continue
            end
        end
        return BondConnectivity
    end
end

function fetch_AtomInfo(lines::Array{String,1},NumberAtoms::Int64,AtomStyle::String)::Array
    for (lineindex,linecontent) in enumerate(lines)
        Phrase=match(Reg_HeadCapitalPhrase,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="Atoms"
                if AtomStyle=="full"
                    AtomInfo=zeros(NumberAtoms,7)
                    
                    for atomii in 1:NumberAtoms
                        atomindex,moleculeindex,atomtype,atomcharge,x,y,z=split(lines[lineindex+1+atomii])
                        selectdim(AtomInfo,1,atomii)[1]=parse(Int,atomindex)
                        selectdim(AtomInfo,1,atomii)[2]=parse(Int,moleculeindex)
                        selectdim(AtomInfo,1,atomii)[3]=parse(Int,atomtype)
                        selectdim(AtomInfo,1,atomii)[4]=parse(Float64,atomcharge)
                        selectdim(AtomInfo,1,atomii)[5]=parse(Float64,x)
                        selectdim(AtomInfo,1,atomii)[6]=parse(Float64,y)
                        selectdim(AtomInfo,1,atomii)[7]=parse(Float64,z)
                    end
                    return AtomInfo
                elseif AtomStyle=="bond"
                    AtomInfo=zeros(NumberAtoms,6)
                    for atomii in 1:NumberAtoms
                        atomindex,moleculeindex,atomtype,x,y,z=split(lines[lineindex+1+atomii])
                        selectdim(AtomInfo,1,atomii)[1]=parse(Int,atomindex)
                        selectdim(AtomInfo,1,atomii)[2]=parse(Int,moleculeindex)
                        selectdim(AtomInfo,1,atomii)[3]=parse(Int,atomtype)
                        selectdim(AtomInfo,1,atomii)[4]=parse(Float64,x)
                        selectdim(AtomInfo,1,atomii)[5]=parse(Float64,y)
                        selectdim(AtomInfo,1,atomii)[6]=parse(Float64,z)
                    end
                    return AtomInfo
                else
                    error("unknow atom_style")
                end
            
            end
        end
    end
    return zeros(NumberAtoms,1)
end

function fetch_AngleConnectivity(lines::Array{String,1},NumberAngles::Int64)::Array
    AngleConnectivity=zeros(Int64,NumberAngles,5)
    for (lineindex,linecontent) in enumerate(lines)
        Phrase=match(Reg_HeadCapitalPhrase,linecontent)
        if Phrase==nothing
            continue
        else
            if Phrase.match=="Angles"
                for angleii in 1:NumberAngles
                    angleindex,angletype,atomindex1,atomindex2,atomindex3=split(lines[lineindex+1+angleii])
                    selectdim(AngleConnectivity,1,angleii)[1]=parse(Int,angleindex)
                    selectdim(AngleConnectivity,1,angleii)[2]=parse(Int,angletype)
                    selectdim(AngleConnectivity,1,angleii)[3]=parse(Int,atomindex1)
                    selectdim(AngleConnectivity,1,angleii)[4]=parse(Int,atomindex2)
                    selectdim(AngleConnectivity,1,angleii)[5]=parse(Int,atomindex3)
                end
                return AngleConnectivity
            else
                continue
            end
        end
        return AngleConnectivity
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


end
