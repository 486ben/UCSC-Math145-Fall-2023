using GLMakie
GLMakie.activate!()
Makie.inline!(true)                                                             # If you are not using VSCode, probably set this to false.
using Roots

function itg(f::Function, x::Real, n::Integer)                  # function that returns an orbit of x under f of size n
    out = zeros(Real, n)
    out[1] = x
    for i in 2:n
        out[i] = f(out[i-1])
    end
    out
end
function cobwebPlot(f::Function, x::Real, n::Integer, r, colorshift = false, spectrum::Symbol = :gist_heat)
    fig = Figure()                                                                                                      # define figure
    ax = Axis(fig[1,1], limits = (r[1], r[end],r[1], r[end])) 
    lines!(ax, r, f, color = :blue)                                                                                     # draw function
    lines!(ax, r, r, color = :black)                                                                                    # draw identity
    l = itg(f, x, n+1)
    extraShift = [0, l[2]]
    diags = Array{Array}(undef, n+1)
    shifts = Array{Array}(undef, n)
    for i in 1:(n+1)
        diags[i] = [l[i],l[i]]
    end
    c(t) = colorshift ? 1 - sqrt(t/(2n)) : 0.6
    for i in 1:n
        shifts[i] = [l[i],l[i+1]]
    end
    lines!(ax, diags[1], extraShift, color = [c(0), c(0)], colorrange = (0,1), colormap = spectrum)                     # draw initial vertical line
        for i in 2:n
            lines!(ax, shifts[i-1], diags[i], color = [c(2i),c(2i+1)], colorrange = (0,1), colormap =  spectrum)        # draw horizontal
            lines!(ax, diags[i], shifts[i], color = [c(2i+1),c(2i+2)], colorrange = (0,1), colormap =  spectrum)        # draw vertical
        end
        lines!(ax, shifts[n], diags[n+1],color = [c(2n),c(2n)], colorrange = (0,1), colormap =  spectrum)               # draw final horizontal 
    fig
end


# In class we saw that x=0 is a repelling fixed point. 
# Let us check (through a list of values) that any point as closed as I want to 0 is repelled. In fact,

g(x,) = 2*x*(1-x)                                  

smallrand(n) = rand()/n                                                         # picks a pseudo-random number between 0 and 1/n


println("Orbits for small values of x")
println("Less than 0.1:")
println(itg(g, smallrand(10), 20)) 
println("Less than 0.01:")
println(itg(g, smallrand(100), 20))   
println("Less than 0.001:") 
println(itg(g, smallrand(1000), 20))       

# Now, let us check that g "behaves" like the function "2x" for points really close to zero by plotting the ratio of the image of a point with the point. 

ratio(a) = g(a)/a
a = smallrand(1000)
println("Ratio g(a)/a for a random 'a' less than 0.001:")
println(ratio(a))

# We saw in class that x=1/2 was an attracting fixed point of the above map g. 
# Moreover, we proved that basin(1/2)=(0,1). Check this by 
# 1) taking random points between (0,1) an see the orbit tends to 1/2 through cobwebplot and/or through a list . 
# 2) taking random points smaller than 0 and bigger than 1, again through cobweb plot and/or through a list. 

z = rand()                                                                      # Generates a pseudorandom number between 0 and 1
cobwebPlot(g, z, 10, 0:0.001:1)


f(x) = (3x - x^3)/2
#We saw in class that 1 is an attracting fixed point of f. 
#Through a list check that any initial value in (0, √3) is in the basin of 1. 

w = rand()*(√3)
println(itg(f, w, 20)) 

range = -2.5:0.01:2.5
fig = Figure(); ax = Axis(fig[1,1])                                             # We write this on one line using the ; operator so that we do not 
ylims!(ax, (-2.5, 2.5))                                                         # accidentally generate a blank image
lines!(ax,range,f)
fig

# Notice that f(-2) is smaller than √3
println(f(-2))

# I explained in class why from these, we can conclude the interval I=(-2,-√3) 
# is in the basin of 1 (because f(I) is contained in (0, √3)). Check it through a list.

b = rand(-2:0.0001:√3)

println(itg(f, b, 20))

# Exercise: Do the same analysis with the fixed attracting point -1 with initial interval (-√3 ,0).

k = rand(-√3:0.0001:1)
println(itg(f, k, 20))
# Notice that f(2) is bigger than -√3
println(f(2))

# Therefore we can conclude that I=(Sqrt[3],2) is in the basin of -1  (because f(I) is contained in (-Sqrt[3],0)). Check it through a list.

c = rand(-√3:0.0001:2) 
itg(f, c, 20)

h(x) = 3.3x*(1-x)

range2 = 0:0.001:1
fig2 = Figure(); ax2 = Axis(fig2[1,1])
ylims!(ax2, (0, 1))
lines!(ax2,range2,h)
lines!(ax2,range2,range2)
fig2

# Check through a list that almost any initial condition 
# between 0 and 1  will be converging to two points. 
d = rand(0:0.001:1) 
println(itg(h, d, 40))

find_zeros((x -> h(x) - x), 0, 1)                                                       # finds zeros of h(x) - x; that is, fixed points of h
find_zeros((x -> h(h(x)) - x), 0, 1)                                                    # finds zeros of h^2(x) - x; that is, fixed points of h^2

j(x) = 3.5x*(1-x)


# Do the same analysis  for J[x]=3.5 x(1-x). It should show a periodic 4-orbit (0.500884,0.874997,0.38282,0.826941)
println(itg(j, d, 40))

find_zeros((x -> j(x) - x), 0, 1) 
find_zeros((x -> j(j(x)) - x), 0, 1) 
find_zeros((x -> j(j(j(x))) - x), 0, 1) 
find_zeros((x -> j(j(j(j(x)))) - x), 0, 1) 