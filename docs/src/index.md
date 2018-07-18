#EMST

```@docs
EMST.compute_emst
```

# Example usage
```
using EMST
x = rand(2,2000) # x = rand(d,n)
x_emst    = EMST.compute_emst(x)

# and plot..
using Plots
using GR

function plot_emst_2d(x,edges)
    Plots.gr()
    Plots.plot()
    Plots.scatter!(x[1,:],x[2,:],ms=2.0,marker=(stroke(0,:gray)))
    for zr in 1:size(edges,1)
        Plots.plot!( x[1,edges[zr,:]] , x[2,edges[zr,:]] )
    end
    Plots.scatter!() # for some reason with this the plot shows im juno..
end

plot_emst_2d(x,x_emst)
```
