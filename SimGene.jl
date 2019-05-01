using Distributions
using Main.SimTree
using Main.SimSpecies

#fix deepcopy deepcopy
function genetree(S::Clade,SL::Real,DL::Real,D::Real)
    #name stuff
    for i in [1,2]
        if S.evoltype == "Speciation"
            lossrate = rand(Exponential(1/SL))
        else S.evoltype == "Duplication"
            lossrate = rand(Exponential(1/DL))
        end
        duplicationrate = rand(Exponential(1/D))
        if duplicationrate == min(duplicationrate, S.branchlength[i], lossrate)
                println("DUPLICATION")
                S.children[i] = Clade("",[deepcopy(S.children[i]),deepcopy(S.children[i])],[S.branchlength[i] - duplicationrate, S.branchlength[i] - duplicationrate],"Duplication")
                S.branchlength[i] = duplicationrate
                genetree(S.children[i],SL,DL,D)
        elseif lossrate == min(duplicationrate, S.branchlength[i], lossrate)
                println("LOSS")
                S.children[i] = Clade("",nothing,nothing,"Loss")
                S.branchlength[i] = lossrate
        elseif isleaf(S.children[i])
        elseif S.branchlength[i] == min(duplicationrate, S.branchlength[i], lossrate)
                println("SPECIATION")
                genetree(S.children[i],SL,DL,D)
        end
    end
    return(S)
end

fourtaxa = Main.SimSpecies.SimSpecies.bdsa(0.01, 0, 4, 10)
println(writenewick(fourtaxa))

geneout = genetree(deepcopy(fourtaxa), 0, 0, 0.1)
println(writenewick(geneout))

function getleafnames(clade::Clade,cladelist::Array=[])
    if istip(clade)
        if isleaf(clade)
            push!(cladelist,clade.node)
            return(cladelist)
        end
    else
        getleafnames(clade.children[1],cladelist)
        getleafnames(clade.children[2],cladelist)
        return(cladelist)
    end
end

function namegenes(clade::Clade, index::Int=1)
    if istip(clade) == true
        if isleaf(clade)
            clade.node = clade.node * "_" * string(index)
        end
    elseif clade.evoltype == "Duplication" #unique(getleafnames(clade.children[1])) == unique(getleafnames(clade.children[2]))
        namegenes(clade.children[1], index)
        namegenes(clade.children[2], index + 1)
    else
        namegenes(clade.children[1], index)
        namegenes(clade.children[2], index)
    end
    return(clade)
end

println(writenewick(namegenes(geneout)))

println(writenewick(geneout))

getleafnames(geneout.children[2].children[2])

for node in getnodes(geneout)
    println(getleafnames(node))
end


for node in getnodes(geneout)
    if istip(node) == false
    if unique(getleafnames(node.children[1])) == getleafnames(node.children[2])
        println(node.branchlength)
    end
    end
end





frog = 4
lizard = frog
frog += 1
lizard

#actually do need tree top for duplications
function genetreetest(G::Clade,SL::Real,DL::Real,D::Real,S::Clade=deepcopy(G))
    #name stuff
    for i in [1,2]
        if G.evoltype == "Speciation"
            lossrate = rand(Exponential(1/SL))
        elseif G.evoltype == "Duplication"
            lossrate = rand(Exponential(1/DL))
        end
        duplicationrate = rand(Exponential(1/D))
        #        if duplicationrate < S.branchlength[i]
        #            println("DUPLICATION")
        #            S.children[i] = Clade("",[S.children[i],S.children[i]],[S.branchlength[i] - duplicationrate, S.branchlength[i] - duplicationrate],"Duplication")
        #            S.branchlength[i] = duplicationrate
        #            genetreetest(S.children[i],SL,DL,D)
        #        end
        if duplicationrate == min(G.branchlength[i], duplicationrate, lossrate)
                println("DUPLICATION")
                G.children[i] = Clade("",[deepcopy(S.children[i]), deepcopy(S.children[i])],[G.branchlength[i] - duplicationrate, G.branchlength[i] - duplicationrate],"Duplication")
                G.branchlength[i] = duplicationrate
                genetreetest(G.children[i],SL,DL,D,S) #lalalookiehere
        elseif lossrate == min(G.branchlength[i], duplicationrate, lossrate)
                println("LOSS")
                G.children[i] = Clade("",nothing,nothing,"Loss")
                G.branchlength[i] = lossrate
                #can I really leave empty?
        #elseif isleaf(S.children[i]) || SimTree.isloss(S.children[i]) || isleaf(G.children[i]) || SimTree.isloss(G.children[i])
        else    #S.branchlength[i] == min(G.branchlength[i], duplicationrate, lossrate)
                println("SPECIATION")
                #maybe change branch lengths
                if isleaf(S.children[i])
                    G.children[i] = deepcopy(S.children[i])
                else
                    genetreetest(G.children[i],SL,DL,D,S.children[i])
                end

        end
    end
    return(G)
end

using Distributions
using Random
using Main.SimTree
using Main.SimSpecies

fourtaxa = Main.SimSpecies.SimSpecies.bdsa(0.01, 0, 4, 100)
println(writenewick(fourtaxa))

geneout = genetreetest(fourtaxa, 0, 0, .01)
println(writenewick(geneout))


function test(Clade1=Clade(),Clade2=Clade())
    Clade1.node = "frog"
    return(Clade1,Clade2)
end

testclade = Clade()
testcopy = copy(testclade)

test(testclade, testcopy)




fourtaxacopy = copy(fourtaxa)
treeprune(fourtaxa)
println(writenewick(fourtaxacopy))

a = 1
b = copy(a)
a = a+1
b

fourtaxa = Main.SimSpecies.SimSpecies.bdsa(0.01, 0, 4, 100)
fourtaxacopy = copy(fourtaxa)
fourtaxa = Clade()
fourtaxacopy
fourtaxa



copy(fourtaxa)

geneout[2]


println(writenewick(treeprune(geneout)))


1/3
1/0.01


species=speciesgen(1,0,nothing,4)
rand(Exponential(1/0))
rand(Exponential(1/1))
0.441 < species.branchlength[1]
gene = Clade(nothing,[Clade(),Clade()],copy(species.branchlength))
gene.branchlength[1] = 0.441
species.children[1] = Clade(nothing,[species.children[1],species.children[1]],[species.branchlength[1]-0.441,species.branchlength[1]-0.441],"duplication")
species.branchlength[1] = 0.441

rand(Exponential(1/0))
rand(Exponential(1/1))
0.623 < species.branchlength[2]
gene = Clade(nothing,[Clade(),Clade()],copy(species.branchlength))
gene.branchlength[2] = 0.623
species.children[2] = Clade(nothing,[species.children[2],species.children[2]],[species.branchlength[2]-0.623,species.branchlength[2]-0.623],"duplication")
species.branchlength[2] = 0.623

println(writenewick(species.children[1]))
println(writenewick(species))
println(writenewick(gene))


species=speciesgen(1,0,nothing,4)
println(writenewick(species))

GT = genetree(species,0.5,0.5)
println(writenewick(GT))

prin


println(writenewick(ssa(100.0,0.5,0.5)))

species.children[2] = Clade(nothing,[species.children[2],species.children[2]],[species.branchlength[2]-0.1,species.branchlength[2]-0.1],"duplication")
species.branchlength[2] = 0.1
println(writenewick(species))


species.branchlength


testG = Clade(nothing,[Clade(),Clade()],copy(testS.branchlength))

testG.children[1] = Clade(nothing,[testS,testS],[testS.branchlength[1]-0.5,testS.branchlength[1]-0.5])

#function to copy tree topology and branch lengths
#extinction rate parameter?

println(writenewick(testG))

nameleaves(testG)

workspace()


"\u000B7"

"🐖"

println("\U1F600")

namelist=["\U1F600", "\U1F601", "\U1F602", "\U1F603", "\U1F604", "\U1F605", "\U1F606", "\U1F607", "\U1F608", "\U1F609", "\U1F60A",
 "\U1F60B", "\U1F60C", "\U1F60D", "\U1F60E", "\U1F60F", "\U1F610", "\U1F611", "\U1F612", "\U1F613", "\U1F614", "\U1F615", "\U1F616",
  "\U1F617", "\U1F618", "\U1F619", "\U1F61A", "\U1F61B", "\U1F61C", "\U1F61D", "\U1F61E", "\U1F61F", "\U1F620", "\U1F621", "\U1F622",
   "\U1F623", "\U1F624", "\U1F625", "\U1F626", "\U1F627", "\U1F628", "\U1F629", "\U1F62A", "\U1F62B", "\U1F62C", "\U1F62D", "\U1F62E",
    "\U1F62F", "\U1F630", "\U1F631", "\U1F632", "\U1F633", "\U1F634", "\U1F635", "\U1F636", "\U1F637"]

println(namelist)

length(namelist)



test=Clade(nothing,[Clade(nothing,[Clade(),Clade(nothing,nothing,nothing,"Loss")],
[1,1]),Clade(nothing,[Clade(nothing,nothing,nothing,"Loss"),Clade(nothing,nothing,nothing,"Loss")],[1.0,1.0])],[1.0,1.0])


test=ssa(10,0.5,0.25)
println(writenewick(test))
test=treepruneint(test)

isdegenerate(test.children[1])
test = test.children[2]



println(writenewick(newtest))
treepruneint(test)
println(writenewick(test))


k=4
N=1000000
a=((k(k-1))/4*N) * exp(-(((k(k-1))/4*N)*t))
