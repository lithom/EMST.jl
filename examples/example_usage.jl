using EMST
using Plots
using GR
function plot_emst_2d(x,edges)
    Plots.gr()
    Plots.plot()
    Plots.scatter!(x[1,:],x[2,:],ms=3.0,marker=(stroke(0,:gray)))
    for zr in 1:size(edges,1)
      Plots.plot!( x[1,edges[zr,:]] , x[2,edges[zr,:]] , linecolor = :gray)
    end
    Plots.scatter!(legend=false) # for some reason with this the plot shows im juno.. :)
end
function plot_emst_3d(x,edges)
    Plots.gr()
    Plots.plot()
    Plots.scatter!(x[1,:],x[2,:],x[3,:],ms=3.0,marker=(stroke(0,:gray)))
    for zr in 1:size(edges,1)
      Plots.plot!( x[1,edges[zr,:]] , x[2,edges[zr,:]] , x[3,edges[zr,:]] , linecolor = :gray)
    end
    Plots.scatter!(legend=false) # for some reason with this the plot shows im juno.. :)
end


x_test  = rand(2,400)
tic(); e_test  = EMST.compute_emst(x_test;nmin=64); toc();
EMST.verify_emst(x_test,e_test,200)
plot_emst_2d(x_test,e_test)

x_test  = rand(3,400)
tic(); e_test  = EMST.compute_emst(x_test;nmin=64); toc();
EMST.verify_emst(x_test,e_test,200)
plot_emst_3d(x_test,e_test)


Plots.savefig("C:\\Temp\\emst_2d.png")
