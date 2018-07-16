#define tree data structure and associated helper functions
module SimTree

    mutable struct Clade
        node::Union{String,Void}
        children::Union{Vector{Clade},Void}
        branchlength::Union{Vector{Real},Void}
        evoltype::Union{String,Void}
        function Clade(a=nothing,b=nothing,c=nothing,d="Speciation")
            if typeof(a) == Void
                a = randstring(4)
            end
            if typeof(b) == Vector{Clade}
                if !(typeof(c) in [Vector{Float64},Vector{Int64},Vector{Real}])
                    error("children do not have branch lengths")
                elseif length(b) != 2
                    error("nodes must be bifurcating")
                elseif length(c) != 2
                    error("nodes must be bifurcating")
                elseif 0 in c
                    error("branch lengths cannot be 0")
                end
            elseif d != "Loss"
                d = "Leaf"
            end
            new(a,b,c,d)
        end
    end

    function isleaf(x::Clade)
        if x.evoltype == "Leaf"
            return(true)
        else
            return(false)
        end
    end

    function getheight(x::Clade,h::Real=0)
        if isleaf(x.children[1])
            h = h + x.branchlength[1]
        else
            h = h + x.branchlength[1]
            getheight(x.children[1],h)
        end
    end

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
        if isleaf(clade)
            push!(cladelist,clade)
            return(cladelist)
        else
            push!(cladelist,clade)
            getnodes(clade.children[1],cladelist)
            getnodes(clade.children[2],cladelist)
            return(cladelist)
        end
    end

    export Clade, getnodes, isleaf, _getnames, nameleaves, getheight

end
