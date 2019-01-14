using Distributions

#define substitution rate
Beta = 0.5

#define branch length
BranchLength = 1.0

#define substitution ratio (jukes cantor i think)
States = ["A","C","G","T"]

#define starting sequence
StartSequence = ["C","A","G"]

#draw exponential distribution of scale parameter (mean) 1/Beta and start clock at 0
d = Exponential(1/Beta)
StartTimePoint = 0

#decide when (if any) changes take place by drawing timepoint from exponential distribution
#decide what change takes place by drawing a A,C,G or T from a uniform distribution
#rinse and repeat until timepoint exceeds branch length
function evolve(x)
    TimePoint = StartTimePoint + rand(d)
    while TimePoint <= BranchLength
        x = sample(States)
        TimePoint += rand(d)
    end
    return(x)
end

#empty array to hold output sequence
EndSequence=[]

#run evolve function on each base in starting sequence and add to output
for i in range(1,length(StartSequence))
    append!(EndSequence,evolve(StartSequence[i]))
    println(EndSequence[i])
end
