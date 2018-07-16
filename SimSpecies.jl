#generate species tree (add emote)
module SimSpecies

    using Distributions
    importall SimTree

    function genspeciesrate(t::Real=10.0,r::Real=0.1,c::Clade=Clade(nothing,[Clade(),Clade()],[t,t]))
        d=Exponential(1/r)
        #branch 1
        split = rand(d)
        if split < t
            newbranch = t-split
            c.children[1] = Clade(nothing,[Clade(),Clade()],[newbranch,newbranch])
            c.branchlength[1] = split
            genspeciesrate(newbranch,r,c.children[1])
        end

        #branch 2
        split = rand(d)
        if split < t
            newbranch = t-split
            c.children[2] = Clade(nothing,[Clade(),Clade()],[newbranch,newbranch])
            c.branchlength[2] = split
            genspeciesrate(newbranch,r,c.children[2])
        end
        nameleaves(c)
        return(c)
    end


    function sprout(dt::Real,c::Clade,extra::Real)
        Branch = rand([1,2])
        if dt < c.branchlength[Branch]
            if c.children[Branch].children == nothing
                c.children[Branch] = Clade(nothing,[Clade(),Clade()],[extra, extra])
                c.branchlength[Branch] = dt
            else
                oldnode = c.children[Branch]
                c.children[Branch] = Clade(nothing,[oldnode,Clade()],[(extra - getheight(oldnode)), extra])
                c.branchlength[Branch] = dt
            end
        else
            sprout((dt - c.branchlength[Branch]),c.children[Branch],extra)
        end
    end


    function genspeciestip(time::Real,n::Integer,outclade::Clade=Clade(nothing,[Clade(),Clade()],[time,time]))
        dis = Uniform(0,time)
        for i in 1:n-2
            dt = rand(dis)
            sprout(dt,outclade,time-dt)
        end
        nameleaves(outclade)
        return(outclade)
    end

    export genspeciesrate, genspeciestip, sprout

end
