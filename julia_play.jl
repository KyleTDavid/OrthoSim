#vector
A = [10,20,30] #row vector
A[1] = 5 #access first element (1-indexed)
A = [10;20;30] #column vector

#matrix
M = [1 2 3; 4 5 6; 7 8 9] #; for rows, \s for columns
M[1,2] = 3663 #access first row, second column
M' #transpose matrix
inv(M) #inverse matrix

#dictionary
D = Dict("name" => "Kyle", "Year" => [4,15,1994]) #create dictionary
D["Year"][3] #access values from key
D.count

#strings
text = "Some String" #create string
println(text[1]) #access first character
text[1] = "L" #strings are immutable

#l00ps
fact=1
for i in range(1,5) # foor loop
    fact=fact*i
end
print(fact)

fact=1 # while loop with conditionals
while fact <= 5
    if fact<=3
        print("cool\n")
    else
        print("hold up\n")
    end
fact+=1
end

function f(x,y) #basic function syntax
    x + y
end
f(1,1)

f(x,y) = x + y #alternatively
f(1,1)

function hypo(x,y)#function to get hypotenuse of right triangle
    x = abs(x)
    y = abs(y)
    if x > y
        r = y/x
        return x*sqrt(1+r*r)
    end
    if y == 0
        return zero(x)
    end
    r = x/y
    return y*sqrt(1+r*r)
end

hypo(12,5)

1 + 2 + 3 #operators are functions
+(1,2,3)

foo(a,b) =  a+b, a*b
foo(2,3)[2]

bar(a,b,x...) = (a,b,x) #variable number of arguments
bar(1,2,3,4,5,6)
x = (3,4)
bar(1,2,x...)

### Bayesian Linear Regression Model Tutorial (https://mambajl.readthedocs.io/en/latest/tutorial.html)###
using Mamba

## model specificiation
model = Model(
    y = Stochastic(1,(mu,s2) -> MvNormal(mu, sqrt(s2)),false),
    mu = Logical(1,(xmat,beta) -> xmat * beta,false),
    beta = Stochastic(1,()->MvNormal(2,sqrt(1000))),
    s2 = Stochastic(() -> InverseGamma(0.001, 0.001)),
    )

## Hybrid No-U-Turn and Slice Sampling Scheme
scheme1 = [NUTS(:beta),
           Slice(:s2, 3.0)]

## No-U-Turn Sampling Scheme
scheme2 = [NUTS([:beta, :s2])]

## User-Defined Samplers

Gibbs_beta = Sampler([:beta],
  (beta, s2, xmat, y) ->
    begin
      beta_mean = mean(beta.distr)
      beta_invcov = invcov(beta.distr)
      Sigma = inv(Symmetric(xmat' * xmat / s2 + beta_invcov))
      mu = Sigma * (xmat' * y / s2 + beta_invcov * beta_mean)
      rand(MvNormal(mu, Sigma))
    end
)

Gibbs_s2 = Sampler([:s2],
  (mu, s2, y) ->
    begin
      a = length(y) / 2.0 + shape(s2.distr)
      b = sum(abs2, y - mu) / 2.0 + scale(s2.distr)
      rand(InverseGamma(a, b))
    end
)

## User-Defined Sampling Scheme
scheme3 = [Gibbs_beta, Gibbs_s2]

#draw DAG
draw(model, filename="lineDAG.dot")

##data
line = Dict{Symbol, Any}(
    :x => [1,2,3,4,5],
    :y => [1,3,3,3,5])

line[:xmat] = [ones(5) line[:x]]

##initial values
inits = [
    Dict{Symbol, Any}(
    :y => line[:y],
    :beta => rand(Normal(0,1), 2),
    :s2 => rand(Gamma(1, 1))
    )
for i in 1:3
]

setsamplers!(model, scheme1)
sim1 = mcmc(model, line, inits, 10000, burnin=250,thin=2,chains=3)

setsamplers!(model, scheme2)
sim2 = mcmc(model, line, inits, 10000, burnin=250, thin=2, chains=3)

setsamplers!(model, scheme3)
sim3 = mcmc(model, line, inits, 10000, burnin=250, thin=2, chains=3)

gelmandiag(sim1,mpsrf=true,transform=true) |> showall
gewekediag(sim1) |> showall

using Gadfly
p = plot(sim1)
display(p)
gui()
draw(p, filename="summaryplot.svg")
