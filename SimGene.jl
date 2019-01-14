importall SimNewick
importall SimSpecies
importall SimTree

using Distributions


tree = readnewick("((A:4.2,B:3):2,C:3):4;")
bigtree = genspeciesrate(10,0.1)
tentree = genspeciestip(10,10)

println(writenewick(tentree))
println(writenewick(genetree(tentree,0.001,0.1)))

println(writenewick(genspeciestip(100,200)))

# +0 solution not super elegant, ask?
function genetree(S::Clade,L::Real,D::Real,G::Clade=Clade(nothing,[Clade(),Clade()],copy(S.branchlength)))
    for i in [1,2]
        lossrate = rand(Exponential(1/L))
        duplicationrate = rand(Exponential(1/D))
        if isleaf(S.children[i])
            G.children[i] = Clade()
        elseif lossrate > S.branchlength[i] && duplicationrate > S.branchlength[i]
                println("SPECIATION")
                G.children[i] = Clade(nothing,[Clade(),Clade()],copy(S.children[i].branchlength))
                genetree(S.children[i],L,D,G.children[i])
        elseif lossrate < duplicationrate
                println("LOSS")
                G.branchlength[i] = lossrate
                G.children[i] = Clade(nothing,nothing,nothing,"loss")
        elseif duplicationrate < lossrate
                println("DUPLICATION")
                G.branchlength[i] = duplicationrate
                println(duplicationrate, " ", S.branchlength[i])
                #needs to replicate descending gene tree
                G.children[i] = Clade(nothing,[S,S],[S.branchlength[i] - duplicationrate, S.branchlength[i] - duplicationrate],"duplication")
                genetree(S.children[i],L,D,G.children[i])
        end
    end
    return(G)
end

testS = genspeciestip(10,3)

testG = Clade(nothing,[Clade(),Clade()],copy(testS.branchlength))

testG.children[1] = Clade(nothing,[testS,testS],[testS.branchlength[1]-0.5,testS.branchlength[1]-0.5])

#function to copy tree topology and branch lengths
#extinction rate parameter?

println(writenewick(testG))

nameleaves(testG)

workspace()


"\ u000B7"

"🐖"

println("\U1F600")

namelist=["\U1F600", "\U1F601", "\U1F602", "\U1F603", "\U1F604", "\U1F605", "\U1F606", "\U1F607", "\U1F608", "\U1F609", "\U1F60A",
 "\U1F60B", "\U1F60C", "\U1F60D", "\U1F60E", "\U1F60F", "\U1F610", "\U1F611", "\U1F612", "\U1F613", "\U1F614", "\U1F615", "\U1F616",
  "\U1F617", "\U1F618", "\U1F619", "\U1F61A", "\U1F61B", "\U1F61C", "\U1F61D", "\U1F61E", "\U1F61F", "\U1F620", "\U1F621", "\U1F622",
   "\U1F623", "\U1F624", "\U1F625", "\U1F626", "\U1F627", "\U1F628", "\U1F629", "\U1F62A", "\U1F62B", "\U1F62C", "\U1F62D", "\U1F62E",
    "\U1F62F", "\U1F630", "\U1F631", "\U1F632", "\U1F633", "\U1F634", "\U1F635", "\U1F636", "\U1F637"]

println(namelist)

length(namelist)
