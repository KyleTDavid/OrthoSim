#generate species tree (add emote)
module SimSpecies

    using Distributions
    using Main.SimTree
    using Random


    #simple sampling approach, used to simulate trees with mrca age t birth rate b and death rate d
    function ssa(t::Real=10.0,λ::Real=0.1,μ::Real=0,c::Clade=Clade("",[Clade(),Clade()],[t,t]))
        birth = Exponential(1/λ)
        death = Exponential(1/μ)
        for branch in [1,2]
            split = rand(birth)
            kill = rand(death)
            if split < t && split < kill
                newbranch = t - split
                c.children[branch] = Clade("",[Clade(),Clade()],[newbranch,newbranch])
                c.branchlength[branch] = split
                ssa(newbranch,λ,μ,c.children[branch])
            elseif kill < t && kill < split
                newbranch = t - kill
                c.children[branch] = Clade("",nothing,nothing,"Loss")
                c.branchlength[branch] = kill
            end
        end
        nameleaves(c)
        return(c)
    end

    #make tree from list of specation times provided by bdsa()
    function _bdsatree(x,y=nothing)
        if y==nothing
            y=repeat([Clade()],length(x))
        end
        min=minimum(x)
        index=findmin(x)[2]
        if index==1
            if isleaf(y[index+1])
                y[index]=Clade("",[Clade(),Clade()],[min,min])
                x[index]=Inf
            else
                rightneighbor=y[index+1]
                y[index]=Clade("",[Clade(),rightneighbor],[min,min-getheight(rightneighbor)])
                x[index]=Inf
                deleteat!(y,index+1)
                deleteat!(x,index+1)
            end
        elseif index==length(x)
            if isleaf(y[index-1])
                y[index]=Clade("",[Clade(),Clade()],[min,min])
                x[index]=Inf
            else
                leftneighbor=y[index-1]
                y[index]=Clade("",[leftneighbor,Clade()],[min-getheight(leftneighbor),min])
                x[index]=Inf
                deleteat!(x,index-1)
                deleteat!(y,index-1)
            end
        else
            rightneighbor=y[index+1]
            leftneighbor=y[index-1]
            if isleaf(leftneighbor)==false && isleaf(rightneighbor)==false
                y[index]=Clade("",[leftneighbor,rightneighbor],[min-getheight(leftneighbor),min-getheight(rightneighbor)])
                x[index]=Inf
                deleteat!(x,[index-1,index+1])
                deleteat!(y,[index-1,index+1])
            elseif isleaf(leftneighbor)==true && isleaf(rightneighbor)==true
                y[index]=Clade("",[Clade(),Clade()],[min,min])
                x[index]=Inf
            elseif isleaf(leftneighbor)==false && isleaf(rightneighbor)
                y[index]=Clade("",[leftneighbor,Clade()],[min-getheight(leftneighbor),min])
                x[index]=Inf
                deleteat!(x,index-1)
                deleteat!(y,index-1)
            elseif isleaf(leftneighbor) && isleaf(rightneighbor)==false
                y[index]=Clade("",[Clade(),rightneighbor],[min,min-getheight(rightneighbor)])
                x[index]=Inf
                deleteat!(y,index+1)
                deleteat!(x,index+1)
            end
        end
        if length(y)==1
            return y[1]
        else
            _bdsatree(x,y)
        end
    end

    #bird-death sampling approach, used to simulate trees with n leaves and/or mrca age t birth rate b and death rate d as long as b>=d
    function bdsa(λ::Real=0.1,μ::Real=0,n::Integer=5,t::Union{Real,Nothing}=nothing)
        r0 = rand(Uniform(0,1))
        speciationtimes = zeros(n-1)
        if t == nothing
            t = (1/(λ-μ))*log((1-(μ/λ)*r0^(1/n))/(1-r0^(1/n)))
        end
        speciationtimes[1] = t
        for st in 2:n-1
            speciationtimes[st] = (1/(λ-μ))*log((λ-μ*ℯ^(-(λ-μ)*t)-μ*(1-ℯ^(-(λ-μ)*t))*rand(Uniform(0,1)))/
            (λ-μ*ℯ^(-(λ-μ)*t)-λ*(1-ℯ^(-(λ-μ)*t))*rand(Uniform(0,1))))
        end
        tree = _bdsatree(shuffle(speciationtimes))
        nameleaves(tree)
        return (tree)
    end

    function speciesgen(mu=1.0,lamda=0.0,time=nothing,n=nothing)
        if time != nothing && n == nothing
            ssa(time, mu, lamda)
        else
            bdsa(mu, lamda, n, time)
        end
    end

end

testarray = []
for i in 1:1000
    append!(testarray, rand(Exponential(1/0.2)))
end
mean(testarray)

using Main.SimTree.SimTree

println(Main.SimTree.writenewick(Main.SimSpecies.SimSpecies.bdsa(0.01, 0, 4, 100)))
println(Main.SimTree.writenewick(Main.SimSpecies.SimSpecies.bdsa(0.01, 0, 8, 100)))
sixteentaxa = Main.SimSpecies.SimSpecies.bdsa(0.01, 0, 16, 100)
println(Main.SimTree.writenewick(sixteentaxa))

println(Main.SimTree.writenewick(Main.SimSpecies.SimSpecies.ssa(100,0.01)))

writenewick(Clade([Clade(),Clade()],[2,3]))

writenewick(Clade("",[Clade("frog"),Clade("dog")],[2,3]))

# #include error from excessively large taxa
# function emojinamer(x::Clade,l::Array=[])
#     index=0
#     for node in getnodes(x)
#         if isleaf(node)
#             push!(l,node)
#         end
#     end
#     names=["\U1F600", "\U1F601", "\U1F602", "\U1F603", "\U1F604", "\U1F605", "\U1F606", "\U1F607", "\U1F608", "\U1F609", "\U1F60A",
#      "\U1F60B", "\U1F60C", "\U1F60D", "\U1F60E", "\U1F60F", "\U1F610", "\U1F611", "\U1F612", "\U1F613", "\U1F614", "\U1F615", "\U1F616",
#       "\U1F617", "\U1F618", "\U1F619", "\U1F61A", "\U1F61B", "\U1F61C", "\U1F61D", "\U1F61E", "\U1F61F", "\U1F620", "\U1F621", "\U1F622",
#        "\U1F623", "\U1F624", "\U1F625", "\U1F626", "\U1F627", "\U1F628", "\U1F629", "\U1F62A", "\U1F62B", "\U1F62C", "\U1F62D", "\U1F62E",
#         "\U1F62F", "\U1F630", "\U1F631", "\U1F632", "\U1F633", "\U1F634", "\U1F635", "\U1F636", "\U1F637"]
#     for leaf in l
#         index+=1
#         leaf.node=string(names[index])
#     end
#     return(x)
# end
