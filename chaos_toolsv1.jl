#MAKE SURE THIS IS IN THE SAME FOLDER AS THE LAB YOU ARE RUNNING
using GLMakie
using Distributions
GLMakie.activate!() # figures will appear as second windows. Set to true to have figures in the VSCode window, but you will not be able to edit sliders
"""
    itg(f::Function, x::Real, n::Integer)

Compute the first n entries in the orbit of x under a real-valued function f.
The first entry is x itself.


# Examples
```julia-repl
julia> itg((x -> 2x), 1, 4)
4-element Vector{Real}:
   1
   2
   4
   8
```
"""
function itg(f::Function, x::Real, n::Integer)                  # function that returns an orbit of x under f of size n
    out = zeros(Real, n)
    out[1] = x
    for i in 2:n
        out[i] = f(out[i-1])
    end
    out
end

"""
fPower(f::Function, n::Integer)

Given a function f and a natural number n, returns the composed function f^n.


# Examples
```julia-repl
julia> f(x) = 2x
f (generic function with 1 method)

julia> g = fPower(f, 10)
f ∘ f ∘ f ∘ f ∘ f ∘ f ∘ f ∘ f ∘ f ∘ f

julia> g(1)
1024
```
"""
function fPower(f::Function, n::Integer)
   return (n <= 1) ? f : f ∘ fPower(f, n-1)
end

# if colorshift is set to true, the cobweb lines will get darker in further iterations
# to highlight the eventual destinations rather than beginning noise
#
# r is the range: give it in the form a:b:c. Here a is the beginning of the range,
# c is the end of the range, and b is the step size or resolution

"""
cobwebPlot(f::Function, x::Real, n::Integer, r::StepRangeLen, colorshift = false)

Return the cobweb plot of a real-valued 1D function f with starting point x and
n iterations over range r. If colorshift is set to true, color the cobweb
line itself in different shades in later iterations according to the parameter spectrum to emphasize convergence. Otherwise
all such lines are a monochrome.
"""
function cobwebPlot(f::Function, x::Real, n::Integer, r::StepRangeLen, colorshift = false, spectrum::Symbol = :gist_heat)
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

"""
cobwebPlots(f::Function, x::Real, n::Integer, r::StepRange, spectrum::Symbol = :gist_rainbow)

Return the cobweb plot of a real-valued 1D function f with all starting points in x and
n iterations over range r. Each line will have a different color according to the 
color range given by parameter spectrum.
"""
function cobwebPlots(f::Function, x::Vector{Real}, n::Integer, r::StepRange, spectrum::Symbol = :gist_rainbow)
    m = length(x)
    c(j) = (j-0.8)/m 
    
    fig = Figure()                                                                                                      # define figure
    ax = Axis(fig[1,1], limits = (r[1], r[end],r[1], r[end])) 
    lines!(ax, r, f, color = :blue)                                                                                     # draw function
    lines!(ax, r, r, color = :black)    
    for j in 1:m                                                                                                        # draw identity
        l = itg(f, x[j], n+1)
        extraShift = [0, l[2]]
        diags = Array{Array}(undef, n+1)
        shifts = Array{Array}(undef, n)
        for i in 1:(n+1)
            diags[i] = [l[i],l[i]]
        end
        for i in 1:n
            shifts[i] = [l[i],l[i+1]]
        end
        lines!(ax, diags[1], extraShift, color = [c(j), c(j)], colorrange = (0,1), colormap = spectrum)                 # draw initial vertical line
            for i in 2:n
                lines!(ax, shifts[i-1], diags[i], color = [c(j),c(j)], colorrange = (0,1), colormap = spectrum)         # draw horizontal
                lines!(ax, diags[i], shifts[i], color = [c(j),c(j)], colorrange = (0,1), colormap = spectrum)           # draw vertical
            end
            lines!(ax, shifts[n], diags[n+1],color = [c(j),c(j)], colorrange = (0,1), colormap = spectrum)              # draw final horizontal 
    end
    fig
end

"""
cobwebScattershot(f::Function, mean::Real, deviation::Real, n::Integer, shots::Integer, r::StepRange, spectrum::Symbol = :gist_rainbow)

Return the cobweb plot of a real-valued 1D function f with starting points in a normal distrubution, and
n iterations over range r. Each line will have a different color according to the 
color range given by parameter spectrum.
"""
function cobwebScattershot(f::Function, mean::Real, deviation::Real, n::Integer, shots::Integer, r::StepRange, spectrum::Symbol = :gist_rainbow)
    d = Normal(mean, deviation)
    x = rand(d, shots)
    cobwebPlots(f, x, n, r, spectrum)
end

"""
cobwebPlotStartBar(f::Function, x::Real, n::Integer, r::StepRange, colorshift = false, spectrum::Symbol = :gist_heat)

Return the cobweb plot of a real-valued 1D function f with starting point given by a slider and
n iterations over range r. If colorshift is set to true, color the cobweb
line itself in different shades in later iterations according to the parameter spectrum to emphasize convergence. Otherwise
all such lines are a monochrome.
"""
# same as previous function, but now the starting point is controlled by a slider bar
function cobwebPlotStartBar(f::Function, n::Integer, r, colorshift = false, spectrum::Symbol = :gist_heat)
    fig = Figure()                                              
    ax = Axis(fig[1,1], limits = (r[1], r[end],r[1], r[end]))  
    lines!(ax, r, f, color = :blue)                                            
    lines!(ax, r, r, color = :black) 
    
    lsgrid = labelslidergrid!(fig,
    ["Starting Point"],
        Ref(LinRange(r)); # same range for every slider via broadcast
        formats = [x -> "$(round(x, digits = 4))" ],
        tellheight = true)
        
    s1 = lsgrid.sliders[1]  
    fig[3,1] = lsgrid.layout                                             # Here a slider is defined
    x = s1.value                                        
    it = lift(x) do t                                                                                   # Don't think too hard about the syntax here 
        itg(f, t, n+1)                                                                                  # It's all in terms of "Observables", which are *evil*
    end
    extraShift = lift(it) do l
        [0,l[2]]
    end
    b = Button(fig, label = "Colorshift", tellheight = true, tellwidth=false)
   
    fig[2, 1] = b
    running = Node(false)

    function toggle_running()
        running[] = !running[] # or more complex logic
        value = running[]
    end
    value = running[]
    on(b.clicks) do x
        toggle_running()
    end
    c(t) = colorshift ? 1 - sqrt(t/(2n)) : 0.6
    diags = Array{Observable}(undef, n+1)
    shifts = Array{Observable}(undef, n)
    for i in 1:(n+1)
        diags[i] = lift(it) do l
            [l[i],l[i]]
        end
    end
    for i in 1:n
        shifts[i] = lift(it) do l
            [l[i],l[i+1]]
        end
    end
    cshift = Array{Observable}(undef, 2n +2)
    for i in 1:(2n+2)
        cshift[i] = lift(running) do o
            [o ? 1 - sqrt(i/(2n)) : 0.6, o ? 1 - sqrt((i+1)/(2n)) : 0.6]
        end
    end
    begincolor = lift(running) do o
        o ? [0,0] : [0.6,0.6]
    end
    endcolor = lift(running) do o
        o ? [1,1] : [0.6,0.6]
    end

    lines!(ax, diags[1], extraShift, color = begincolor, colorrange = (0,1), colormap = spectrum)                       # draw initial vertical line
        for i in 2:n
            lines!(ax, shifts[i-1], diags[i], color = cshift[2i], colorrange = (0,1), colormap = spectrum)           # draw horizontal
            lines!(ax, diags[i], shifts[i], color = cshift[2i+1], colorrange = (0,1), colormap = spectrum)           # draw vertical
        end
        lines!(ax, shifts[n], diags[n+1], color = endcolor, colorrange = (0,1), colormap = spectrum)                 # draw final horizontal 
    fig
end

"""
cobwebPlotDoubleSlider(f::Function, x::Real, n::Integer, r::StepRange, colorshift = false, spectrum::Symbol = :gist_heat)

Return the cobweb plot with 
n iterations over range r of a real-valued 1D function in a real-parametrized family of functions F with starting point and parameter given by sliders
ranging over r and fr respectively. 
If colorshift is set to true, color the cobweb
line itself in different shades in later iterations according to the parameter spectrum to emphasize convergence. Otherwise
all such lines are a monochrome.
"""
function cobwebPlotDoubleSlider(F::Function, n::Integer, r, fr, colorshift = false, spectrum::Symbol = :gist_heat)
    fig = Figure()                                              # define figure
    ax = Axis(fig[1,1], limits = (r[1], r[end],r[1], r[end]))  

    lsgrid = labelslidergrid!(fig, ["Starting Point", "Parameter"], [r, fr]; tellheight = true)
    s1 = lsgrid.sliders[1]  
    s2 = lsgrid.sliders[2]  
    fig[3,1] = lsgrid.layout 


    b = Button(fig, label = "Colorshift", tellheight = true, tellwidth=false)
   
    fig[2, 1] = b
    running = Node(false)

    function toggle_running()
        running[] = !running[] # or more complex logic
        value = running[]
    end
    value = running[]
    on(b.clicks) do x
        toggle_running()
    end
    f = lift(s2.value) do g
        x -> F(g, x)
    end 

    c = lift(running) do o
        (t -> o ? 1 - sqrt(t/(2n)) : 0.6)
    end
    
    lines!(ax, r, f, color = :blue)                                            
    lines!(ax, r, r, color = :black)                 
    x = s1.value                                        
    it = @lift(itg(t -> F($(s2.value), t), $x, n+1))                                                                        #this in particular is painful syntax.
    extraShift = lift(it) do l
        [0,l[2]]
    end
    diags = Array{Observable}(undef, n+1)
    shifts = Array{Observable}(undef, n)
    for i in 1:(n+1)
        diags[i] = lift(it) do l
            [l[i],l[i]]
        end
    end
    for i in 1:n
        shifts[i] = lift(it) do l
            [l[i],l[i+1]]
        end
    end

    cshift = Array{Observable}(undef, 2n +2)
    for i in 1:(2n+2)
        cshift[i] = lift(running) do o
            [o ? 1 - sqrt(i/(2n)) : 0.6, o ? 1 - sqrt((i+1)/(2n)) : 0.6]
        end
    end
    begincolor = lift(running) do o
        o ? [0,0] : [0.6,0.6]
    end
    endcolor = lift(running) do o
        o ? [1,1] : [0.6,0.6]
    end

    lines!(ax, diags[1], extraShift, color = begincolor, colorrange = (0,1), colormap = spectrum)                       # draw initial vertical line
        for i in 2:n
            lines!(ax, shifts[i-1], diags[i], color = cshift[2i], colorrange = (0,1), colormap = spectrum)           # draw horizontal
            lines!(ax, diags[i], shifts[i], color = cshift[2i+1], colorrange = (0,1), colormap = spectrum)           # draw vertical
        end
        lines!(ax, shifts[n], diags[n+1], color = endcolor, colorrange = (0,1), colormap = spectrum)                 # draw final horizontal 
    fig
end

"""
bifurcation(F::Function, insideIterations::Integer, range::StepRangeLen, startRange::StepRangeLen)

Return the bifurcation of a real-parametrized family of functions F. The parameter varies over the parameter range.
Searches for orbits over the range startRange using insideIterations number of iterations.
"""
function bifurcation(F::Function, insideIterations::Integer, range::StepRangeLen, startRange::StepRangeLen)
    function periodFinder(f, n, r)
        bigarray = zeros(Float32, (length(r), n))
        for i in 1:length(r)
            bigarray[i, :] = itg(f, r[i], n)
        end
        union(bigarray[:, end-3:end])
    end 
    fig = Figure()
    ax = Axis(fig[1,1])
    xlims!(ax, (range[1], range[end]))   
    orbits = [0.0]
    rs = [0.0]
    for i in 1:length(range)
        of = periodFinder((x -> F(range[i], x)), insideIterations, startRange)
        for o in of
            append!(orbits, [o])
            append!(rs, [range[i]])
        end
    end
    scatter!(ax, rs, orbits, color = :black, markersize = 1)
    fig
end