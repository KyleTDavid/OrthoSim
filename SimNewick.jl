#read and write newick tree format
module SimNewick

    using DataStructures
    import Main.SimTree: Clade

    function _writenewickint(x::Clade)

        newick="("
        if x.children[1].children == nothing
            newick = newick * x.children[1].node * ":" * string(x.branchlength[1]) * ","
        else newick = newick * _writenewickint(x.children[1]) * string(x.children[1].node) * ":" * string(x.branchlength[1]) * ","
        end

        if x.children[2].children == nothing
            newick = newick * x.children[2].node * ":" * string(x.branchlength[2]) * ")"
        else newick = newick * _writenewickint(x.children[2]) * string(x.children[2].node) * ":" * string(x.branchlength[2]) * ")"
        end

        return(newick)
    end

    function writenewick(x)
        out = _writenewickint(x)
        return(out * x.node *";")
    end


    #function to read newick format into internal Clade data structure
    #currently trees must be rooted, bifurcating, with names at the tips
    #if present, names must be redundant
    function readnewick(x)
        #assign name to internal nodes when not provided
        newtestlist = split(x,r"(\):|\);)")
        newtest = ""
        index=0
        for n in newtestlist
            index += 1
            if index == endof(newtestlist)
                test = n
            else
                test = n * ")"randstring(4)":"
            end
            newtest = newtest*test
        end
        if newtest[end] == ':'
            newtest = replace(newtest,r"(?!$):$",";")
        end

        #get levels
        L = []
        nest = 2
        for C in newtest
            if string(C) == "("
                nest+=1
            elseif string(C) == ")"
                nest-=1
            elseif length(L) == 0
                push!(L,nest)
            elseif nest != L[end]
                push!(L,nest)
            end
        end

        #get substrings
        S = filter!(e -> e!="",split(newtest,r"[\(\)]"))
        rootname = split(S[end],";")[1]

        #put levels and substrings together into dictionary
        D = Dict()
        for i in (1:length(L))
            D[i] = [S[i],L[i]]
        end

        #get nodenames, branch lengths, and children for each node
        #i=level S=string_index sub=substring
        cladedict=OrderedDict()
        for i in unique(sort(L,rev=true)) #for each level
            LevelList = filter((k,v) -> v[2]==(i), D) #list of strings on current level
            UpLevelList = filter((k,v) -> v[2]==(i-1), D) #list of strings on upper level
            for S in keys(LevelList) #iterate through strings in level i
                sub = D[S][1]
                if contains(sub,",")
                    if any([k>S for k in keys(UpLevelList)])
                        NUM1 = minimum(keys(filter((k,v) -> k>S,UpLevelList)))
                        name = split(UpLevelList[NUM1][1],r"[:;]")[1]
                    else
                        name = ""
                    end
                    branch1 = float(split(split(sub,",")[1],":")[2])
                    childname1 = string(split(split(sub,",")[1],":")[1])
                    if split(sub,",")[2] == ""
                        NUM = minimum(keys(filter((k,v) -> k>S,LevelList)))
                        branch2 = float(split(split(LevelList[NUM][1],":")[2],",")[1])
                        childname2 = string(split(split(LevelList[NUM][1],":")[1],",")[1])
                    else
                        branch2 = float(split(split(sub,",")[2],":")[2])
                        childname2 = string(split(split(sub,",")[2],":")[1])
                    end
                    cladedict[name] = [Clade(name,nothing,[branch1,branch2]),[childname1,childname2]]
                end
            end
        end

        #iteratively build tree based on shared names between nodes and children
        for key in keys(cladedict)
            branch1 = cladedict[key][2][1]
            if branch1 in [k for k in keys(cladedict)]
                child1 = cladedict[branch1][1]
            else
                child1 = Clade(branch1,nothing,nothing)
            end
            branch2 = cladedict[key][2][2]
            if branch2 in [k for k in keys(cladedict)]
                child2 = cladedict[branch2][1]
            else
                child2 = Clade(branch2,nothing,nothing)
            end
            cladedict[key][1]=Clade(cladedict[key][1].node,[child1,child2],cladedict[key][1].branchlength)
        end
        return(cladedict[rootname][1])
    end

    export readnewick, _writenewickint, writenewick

end
