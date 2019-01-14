#generate species tree (add emote)
module SimSpecies

    using Distributions
    import Main.SimTree: Clade, nameleaves
    import Main.SimNewick: writenewick

    #simple sampling approach, used to simulate trees with mrca age t birth rate b and death rate d
    function ssa(t::Real=10.0,b::Real=0.1,d::Real=0,c::Clade=Clade(nothing,[Clade(),Clade()],[t,t]))
        birth = Exponential(1/b)
        death = Exponential(1/d)
        for branch in [1,2]
            split = rand(birth)
            kill = rand(death)
            if split < t && split < kill
                newbranch = t - split
                c.children[branch] = Clade(nothing,[Clade(),Clade()],[newbranch,newbranch])
                c.branchlength[branch] = split
                ssa(newbranch,b,d,c.children[branch])
            end
        end
        nameleaves(c)
        return(c)
    end

    #bird-death sampling approach, used to simulate trees with n leaves and/or mrca age t birth rate b and death rate d as long as b>=d
    function bdsa(b::Real=0.1,d::Real=0,n::Integer=5,t::Union{Real,Nothing}=nothing)
        r0 = rand(Uniform(0,1))
        speciationtimes = zeros(n-1)
        if t == nothing
            t = (1/(b-d))*log((1-(d/b)*r0^(1/n))/(1-r0^(1/n)))
        end
        speciationtimes[1] = t
        for st in 2:n-1
            speciationtimes[st] = (1/(b-d))*log((b-d*ℯ^(-(b-d)*t)-d*(1-ℯ^(-(b-d)*t))*rand(Uniform(0,1)))/
            (b-d*ℯ^(-(b-d)*t)-b*(1-ℯ^(-(b-d)*t))*rand(Uniform(0,1))))
        end
        return(speciationtimes)
    end


    function bdsa_tree(x,y=nothing)
        if y==nothing
            y=repeat([Clade()],length(x))
        end
        min=minimum(x)
        index=findmin(x)[2]
        if index==1
            if isleaf(y[index+1])
                y[index]=Clade(nothing,[Clade(),Clade()],[min,min])
                x[index]=Inf
            else
                rightneighbor=y[index+1]
                y[index]=Clade(nothing,[Clade(),rightneighbor],[min,min-getheight(rightneighbor)])
                x[index]=Inf
                deleteat!(y,index+1)
                deleteat!(x,index+1)
            end
        elseif index==length(x)
            if isleaf(y[index-1])
                y[index]=Clade(nothing,[Clade(),Clade()],[min,min])
                x[index]=Inf
            else
                leftneighbor=y[index-1]
                y[index]=Clade(nothing,[leftneighbor,Clade()],[min-getheight(leftneighbor),min])
                x[index]=Inf
                deleteat!(x,index-1)
                deleteat!(y,index-1)
            end
        else
            rightneighbor=y[index+1]
            leftneighbor=y[index-1]
            if isleaf(leftneighbor)==false && isleaf(rightneighbor)==false
                y[index]=Clade(nothing,[leftneighbor,rightneighbor],[min-getheight(leftneighbor),min-getheight(rightneighbor)])
                x[index]=Inf
                deleteat!(x,[index-1,index+1])
                deleteat!(y,[index-1,index+1])
            elseif isleaf(leftneighbor)==true && isleaf(rightneighbor)==true
                y[index]=Clade(nothing,[Clade(),Clade()],[min,min])
                x[index]=Inf
            elseif isleaf(leftneighbor)==false && isleaf(rightneighbor)
                y[index]=Clade(nothing,[leftneighbor,Clade()],[min-getheight(leftneighbor),min])
                x[index]=Inf
                deleteat!(x,index-1)
                deleteat!(y,index-1)
            elseif isleaf(leftneighbor) && isleaf(rightneighbor)==false
                y[index]=Clade(nothing,[Clade(),rightneighbor],[min,min-getheight(rightneighbor)])
                x[index]=Inf
                deleteat!(y,index+1)
                deleteat!(x,index+1)
            end
        end
        if length(y==1)
            return y[1]
        else
            bdsa_tree(x,y)
        end
    end

    test = bdsa()




    out=bdsa_tree(test)
    bdsa_tree(test,out[2])
    bdsa_tree(test,out[2])
    bdsa_tree(test,out[2])


    println(writenewick(out[2][1]))








out1


t1 = Clade(nothing,[Clade(),Clade()],[1,1])
t2 = Clade()
t3 = Clade()


hey=findfirst(x->isleaf(x),[t1,t2,t3])
hey

out1 = bdsa_tree(test2,out1[2])
println(writenewick(tallestclade(out1[2])))














out3 = bdsa_tree(test)
out4 = bdsa_tree(test)
out5 = bdsa_tree(test)
out6 = bdsa_tree(test)
out7 = bdsa_tree(test)
out8 = bdsa_tree(test)
out9 = bdsa_tree(test)






println(writenewick(out9[2][9]))






bdsa_tree(test)
bdsa_tree(test)
bdsa_tree(test)


tallestclade(out2[2])
getheight(out2[2][1])
out2[2][1]


all(x->x==Inf,out9[1])



testout = bdsa_tree(test)
println(writenewick(testout))


if all(x->x==Inf,x)
    return y[minindex]
else
    bdsa_tree(x,y)
end






        if findmin(y)[2]==1
            if typeof(x[2])==Clade
                splice!(x,minindex)
                insert!(x,minindex,Clade(nothing,[Clade(),x[2]],[min,min-getheight(x[2])]))
            else
                splice!(x,minindex)
                insert!(x,minindex,Clade(nothing,[Clade(),Clade()],[min,min]))
            end
        elseif findmin(y)[2]==length(y)
            if typeof(x[end-1])==Clade
                splice!(x,minindex)
                insert!(x,minindex,Clade(nothing,[Clade(),x[end-1]],[min,min-getheight(x[end-1])]))
            else
                splice!(x,minindex)
                insert!(x,minindex,Clade(nothing,[Clade(),Clade()],[min,min]))
            end
        elseif typeof(x[minindex-1])==Float64 && typeof(x[minindex+1])==Float64
            splice!(x,minindex)
            insert!(x,minindex,Clade(nothing,[Clade(),Clade()],[min,min]))
        elseif typeof(x[minindex-1])==Clade && typeof(x[minindex+1])==Clade
            splice!(x,minindex)
            insert!(x,minindex,Clade(nothing,[x[minindex-1],x[minindex+1]],[min-getheight(leftneighbor),min-getheight(rightneighbor)])
        elseif typeof(x[minindex-1])==Float64 && typeof(x[minindex+1])==Clade
            splice!(x,minindex)
            insert!(x,minindex,Clade(nothing,[Clade(),x[minindex+1]],[min,min-getheight(x[minindex+1])]))
        elseif typeof(x[minindex-1])==Clade && typeof(x[minindex+1])==Float64
            splice!(x,minindex)
            insert!(x,minindex,Clade(nothing,[Clade(),x[minindex+1]],[min-getheight(x[minindex-1]),min]))
        end
        splice!(y,minindex,Inf)
        return x
    end

fill(6,Clade())
repeat([Clade()],5)

using Distributions
test = bdsa()
test2 = copy(test)
out1 = bdsa_tree(test,test2)

out2 = bdsa_tree(out1,test2)

out3 = bdsa_tree(out2,test2)



if all(l->typeof(l)==Clade,x)
    return x[minindex]
else
    splice!(y,minindex,Inf)
    bdsa_tree(x,y)
end

list=[1,2,3]

append!([Inf],prepend!([Inf],list))



test

typeof(testarray[2])==Clade

    typeof(Clade())==Clade

    minimum(["frog",4])

    init_bdsa_tree(testarray)

    testarray=[1,2,4,42,5,23,4]

    append!(Clade(),testarray)
    splice!()

    splice!(testarray,2)
    splice!(testarray,Clade(nothing,[Clade(),Clade()],[minimum(testarray),minimum(testarray)]))

    Clade(nothing,[Clade(),Clade()],[minimum(testarray),minimum(testarray)])



    testclade = Clade(nothing,[Clade(),Clade()],[testarraysort[end]-testarraysort[end-1],testarraysort[end]-testarraysort[end-2]])
    filter!(i->i∉[testarraysort[end],testarraysort[end-1],testarraysort[end-2]],testarraysort)

    testclade

    foo(x=1, y=1, z) = x + y + z

        end
        for st in []
    end

    bdsa()

    [1:4]

    b = 0.1
    d = 0
    n = 5


    t = (1/(b-d))*log((1-(d/b)*r0^(1/n))/(1-r0^(1/n)))

    s = (1/(b-d))*log((b-d*e^(-(b-d)*t)-d*(1-e^(-(b-d)*t))*r1)/(b-d*e^(-(b-d)*t)-b*(1-e^(-(b-d)*t))*r1))





#include error from excessively large taxa
function emojinamer(x::Clade,l::Array=[])
    index=0
    for node in getnodes(x)
        if isleaf(node)
            push!(l,node)
        end
    end
    names=["\U1F600", "\U1F601", "\U1F602", "\U1F603", "\U1F604", "\U1F605", "\U1F606", "\U1F607", "\U1F608", "\U1F609", "\U1F60A",
     "\U1F60B", "\U1F60C", "\U1F60D", "\U1F60E", "\U1F60F", "\U1F610", "\U1F611", "\U1F612", "\U1F613", "\U1F614", "\U1F615", "\U1F616",
      "\U1F617", "\U1F618", "\U1F619", "\U1F61A", "\U1F61B", "\U1F61C", "\U1F61D", "\U1F61E", "\U1F61F", "\U1F620", "\U1F621", "\U1F622",
       "\U1F623", "\U1F624", "\U1F625", "\U1F626", "\U1F627", "\U1F628", "\U1F629", "\U1F62A", "\U1F62B", "\U1F62C", "\U1F62D", "\U1F62E",
        "\U1F62F", "\U1F630", "\U1F631", "\U1F632", "\U1F633", "\U1F634", "\U1F635", "\U1F636", "\U1F637"]
    for leaf in l
        index+=1
        leaf.node=string(names[index])
    end
    return(x)
end

testree2=genspeciesrate(10,0.2)
testree=genspeciestip2(5,50)
println(writenewick(emojinamer(testree2)))

importall SimNewick
importall SimTree
importall SimSpecies
using Distributions
println(writenewick(emojinamer(genspeciestip2(10,50))))
println(writenewick(emojinamer(genspeciesrate(10,0.1))))
