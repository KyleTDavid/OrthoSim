#define tree data structure and associated helper functions

module SimTree
    using Random

    #structure for working with tress
    mutable struct Clade
        node::Union{String,Nothing} #optional node name
        children::Union{Vector{Clade},Nothing} #optional children
        branchlength::Union{Vector{Real},Nothing} #child branch lengths, required if children present
        evoltype::Union{String,Nothing} #to track ortho/paralogs, speciation by default
        function Clade(a="",b=nothing,c=nothing,d="Speciation")
            if typeof(b) == Vector{Clade}
                if !(typeof(c) in [Vector{Float64},Vector{Int64},Vector{Real}])
                    error("children do not have branch lengths")
                elseif 0 in c
                    error("branch lengths cannot be 0")
                end
            elseif d != "Loss"
                d = "Leaf"
            end
            new(a,b,c,d)
        end
    end

    Base.copy(s::Clade) = Clade(s.node, s.children, s.branchlength, s.evoltype)


    #internal function to write Clade structures to newick format
    function _writenewick(x::Clade)
        newick="("
        if x.children[1].children == nothing
            newick = newick * x.children[1].node * ":" * string(x.branchlength[1]) * ","
        else newick = newick * _writenewick(x.children[1]) * string(x.children[1].node) * ":" * string(x.branchlength[1]) * ","
        end

        if x.children[2].children == nothing
            newick = newick * x.children[2].node * ":" * string(x.branchlength[2]) * ")"
        else newick = newick * _writenewick(x.children[2]) * string(x.children[2].node) * ":" * string(x.branchlength[2]) * ")"
        end
        return(newick)
    end

    #adds root info to _writenewick output
    function writenewick(clade)
        return(_writenewick(clade) * clade.node * ";")
    end

    #internal function called by readnewick
    function _readnewick(testarray, index)
        name = split(testarray[index][1],":")[1]
        if  index != length(testarray) &&
            testarray[index][2]+1 == testarray[index+1][2]
            child1 = split(testarray[index+1][1],":")
            branchlength1 = parse(Float64,child1[2])
            child2index = findnext(x -> x[2] == testarray[index][2]+1, testarray, index+2)
            child2 = split(testarray[child2index][1],":")
            branchlength2 = parse(Float64,child2[2])
            return Clade(name, [_readnewick(testarray, index+1), _readnewick(testarray, child2index)], [branchlength1, branchlength2])
        else
            return Clade(name)
        end
    end

    #reads newick format into Clade data structure
    function readnewick(newick::String)
        newicksplit = split(newick, r"(?=;)|(?=\()|(?<=\()|(?=\))|(?<=\))|,")
        nest=0
        newickarray = []
        for char in reverse(newicksplit)
            if char == "("
                nest-=1
            elseif char == ")"
                nest+=1
            else
                push!(newickarray, [char,nest])
            end
        end
        _readnewick(newickarray, 1)
    end

    function isleaf(x::Clade)
        if x.evoltype == "Leaf"
            return(true)
        else
            return(false)
        end
    end

    function isloss(x::Clade)
        if x.evoltype == "Loss"
            return(true)
        else
            return(false)
        end
    end

    function istip(x::Clade)
        if x.children == nothing
            return(true)
        else
            return(false)
        end
    end

    function isdegenerate(x::Clade)
        y = gettips(x)
        if all(t->isloss(t),y)
            return(true)
        else
            return(false)
        end
    end

    function getheight(x::Clade,h::Real=0)
        if isleaf(x)
            return 0
        elseif isleaf(x.children[1])
            h = h + x.branchlength[1]
            return h
        else
            h = h + x.branchlength[1]
            getheight(x.children[1],h)
        end
    end

    function tallestclade(x::Array{Clade})
        heightlist=[]
        for clade in x
            append!(heightlist,[getheight(clade)])
        end
        return x[findmax(heightlist)[2]]
    end

    #internal function called by nameleaves
    function _getnames(x::Array,E::Int64=1,l::Array=[])
        if length(x) <= 26^E
            AZlist = collect('A':'Z')
            AZlist = map(string,AZlist)
            for combo in Iterators.product(fill(AZlist,(1,E))...)
                push!(l,combo)
            end
            return(sort(l))
        else
            E += 1
            _getnames(x,E,l)
        end
    end

    function nameleaves(x::Clade,l::Array=[])
        index=0
        for node in getnodes(x)
            if isleaf(node)
                push!(l,node)
            end
        end
        names=_getnames(l)
        for leaf in l
            index+=1
            leaf.node=string("Species"*join(names[index]))
        end
        return(x)
    end

    function getnodes(clade::Clade,cladelist::Array=[])
        if istip(clade)
            push!(cladelist,clade)
            return(cladelist)
        else
            push!(cladelist,clade)
            getnodes(clade.children[1],cladelist)
            getnodes(clade.children[2],cladelist)
            return(cladelist)
        end
    end

    function gettips(clade::Clade,cladelist::Array=[])
        if istip(clade)
            push!(cladelist,clade)
            return(cladelist)
        else
            gettips(clade.children[1],cladelist)
            gettips(clade.children[2],cladelist)
            return(cladelist)
        end
    end

    function getleaves(clade::Clade,cladelist::Array=[])
        if isleaf(clade)
            push!(cladelist,clade)
            return(cladelist)
        else
            getleaves(clade.children[1],cladelist)
            getleaves(clade.children[2],cladelist)
            return(cladelist)
        end
    end

    #internal function called by treeprune
    function _treeprune(clade::Clade)
        for i in [1,2]
            if istip(clade)==false && istip(clade.children[i])==false
                if isdegenerate(clade.children[i].children[1])
                    _treeprune(clade.children[i])
                    clade.branchlength[i]+=clade.children[i].branchlength[2]
                    clade.children[i]=clade.children[i].children[2]
                elseif isdegenerate(clade.children[i].children[2])
                    _treeprune(clade.children[i])
                    clade.branchlength[i]+=clade.children[i].branchlength[1]
                    clade.children[i]=clade.children[i].children[1]
                else
                    _treeprune(clade.children[i])
                end
            end
        end
        return clade
    end

    #removes losses from Clade structure
    function treeprune(clade::Clade)
        if istip(clade.children[1]) && istip(clade.children[2])
            error("only <=2 lineages found")
        end
        if isdegenerate(clade.children[2])
            clade=clade.children[1]
            treeprune(clade)
        elseif isdegenerate(clade.children[1])
            clade=clade.children[2]
            treeprune(clade)
        else
            _treeprune(clade)
        end
    end


    export Clade, _writenewick, writenewick, _readnewick, readnewick, getnodes,
    isleaf, istip, _getnames, nameleaves, getheight, treeprune, _treeprune,
    gettips, getnodes

end


collect('α':'ω')
