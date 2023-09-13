using GLMakie                                                       # old friend!
GLMakie.activate!()



#= 
        This is the Mandelbrot set function. Let c be the complex number x0 + i*y0, and let 
    z be the complex number x + i*y. This function returns z^2 + c. All the fractals in this
    program use this function in some way.
=#
Mandel(x0, y0, x, y) = [x^2 - y^2, 2x*y] + [x0, y0]                  



#=
        Fix a point c on the complex plane. Let f_c be the function 
    z -> z^2 + c. The number c is in the Mandelbrot set if (f_c)^n(0) 
    does not go to infinity. The below function checks if the orbit 
    goes to infinity by looking up to n steps. If it diverges, it returns
    the number of steps it takes for the point 0 to have an absolute value above 2. 
    The reason to return this value is to color the graph later.
=#
function mandelLeaveTime(x, y, n, bound)
    F(v) = Mandel(x, y, v[1], v[2])
    point = [x,y]
    b2 = bound^2
    if x^2 + y^2 >=  b2
        return 0
    end
    for i in 1:n
        point = F(point)
        if point[1]^2 + point[2]^2 >=  b2
            return i
        end
    end
    return n+1
end

# some functions which make color gradients work
makeColorTracker(n) = t -> ((n+ 1 -t)/(n+1))^6
c = makeColorTracker(300)

cSquare(x,y) = [x^2 - y^2, 2x*y]


#=
        The Julia set at a point c is defined similarly. Define f_c as before.
    The Julia set consists of the points z such that (f_c)^n(z) does not go
    to infinity. The value returned is decided in the same way as before.
=#
function juliaLeaveTime(x0, y0, x, y, n, bound)
    F(v) = Mandel(x0, y0, v[1],v[2])
    point = [x,y]
    b2 = bound^2
    if x^2 + y^2 >=  b2
        return 0
    end
    for i in 1:n
        point = F(point)
        if point[1]^2 + point[2]^2 >=  b2
            return i
        end
    end
    return n+1
end



#=
    The Burning Ship fractal is defined similarly to the Mandelbrot set, except that we 
take Mandel(x0, y0, |x|, -|y|) instead. See this link for more info.

https://en.wikipedia.org/wiki/Burning_Ship_fractal

It's my personal favorite fractal as of right now.
=#
bShip(x0, y0, x, y) = Mandel(x0, y0, abs(x), -abs(y)) 
function bShipLeaveTime(x, y, n, bound)
    F(v) = bShip(x, y, v[1], v[2])
    point = [x,y]
    b2 = bound^2
    if x^2 + y^2 >=  b2
        return 0
    end
    for i in 1:n
        point = F(point)
        if point[1]^2 + point[2]^2 >=  b2
            return i
        end
    end
    return n+1
end


# The Burning Ship has Julia sets in just the same way as the Mandelbrot set.

function bJuliaLeaveTime(x0, y0, x, y, n, bound)
    F(v) = bShip(x0, y0, v[1], v[2])
    point = [x,y]
    b2 = bound^2
    if x^2 + y^2 >=  b2
        return 0
    end
    for i in 1:n
        point = F(point)
        if point[1]^2 + point[2]^2 >=  b2
            return i
        end
    end
    return n+1
end


# Generate a graph of the mandelbrot set with x values in the range xr,
# y values in the range yr, and a maximum number of n iterations.

function mandelbrot(xr, yr, n, gradient = :thermal) 
    c = makeColorTracker(n)
    GLMakie.heatmap(xr, yr, (x,y) -> c(mandelLeaveTime(x, y, n, 2)), colormap = gradient)
end

# Generate a graph of the Burning Ship with x values in the range xr,
# y values in the range yr, and a maximum number of n iterations.

function burningShip(xr, yr, n, gradient = :thermal) 
    c = makeColorTracker(n)
    GLMakie.heatmap(xr, yr, (x,y) -> c(bShipLeaveTime(x, y, n, 2)), colormap = gradient)
end


# Generate a graph of the Julia set with c = (real + i*imaginary),
# x values in the range xr, y values in the range yr, and a maximum number of n iterations.

function juliaSet(real, imaginary, xr, yr, n, gradient = :thermal)
    c = makeColorTracker(n)
    GLMakie.heatmap(xr, yr, (x,y) -> c(juliaLeaveTime(real, imaginary, x, y, n, 2)), colormap = gradient)
end

# Generate a graph of the Bunrning Ship's Julia set with c = (real + i*imaginary),
# x values in the range xr, y values in the range yr, and a maximum number of n iterations.

function bJuliaSet(real, imaginary, xr, yr, n, gradient = :thermal)
    c = makeColorTracker(n)
    GLMakie.heatmap(xr, yr, (x,y) -> c(bJuliaLeaveTime(real, imaginary, x, y, n, 2)), colormap = gradient)
end

function geneMandelLeaveTime(M, x, y, n, bound)
    F(v) = M(x, y, v[1], v[2])
    point = [x,y]
    b2 = bound^2
    if x^2 + y^2 >=  b2
        return 0
    end
    for i in 1:n
        point = F(point)
        if point[1]^2 + point[2]^2 >=  b2
            return i
        end
    end
    return n+1
end

function geneMandel(M, xr, yr, n, gradient = :thermal) 
    c = makeColorTracker(n)
    GLMakie.heatmap(xr, yr, (x,y) -> c(geneMandelLeaveTime(M, x, y, n, 2)), colormap = gradient)
end

# WARNING: THESE COMPUTATIONS TAKE A LOT OF CPU TIME AND RAM
# START WITH SMALL STEP SIZES IN THE RANGE OR YOU WILL REGRET IT


# The classic image of the Mandelbrot set
mandelbrot(-2:0.003:0.5, -1.2:0.003:1.2, 100)

# The Mandelbrot set zoomed way in 
# Here an alternate gradient is used. You can find these gradient names at 
#
# https://docs.juliahub.com/AbstractPlotting/6fydZ/0.12.10/generated/colors.html
#
mandelbrot(-0.467:0.00001:-0.4625, 0.552:0.00001:0.556, 700, :seaborn_rocket_gradient)


# The most iconic Julia set
juliaSet(-0.79, 0.16, -1.5:0.003:1.5, -1:0.003:1, 200)

# A fun Julia set exhibiting 13-fold spiral galaxy shapes. Try picking values of 
# c close to the border of the Mandelbrot set!
juliaSet(-0.408, 0.595, -1.5:0.003:1.5, -1:0.003:1, 200,:PuBu_3)

# The same set zoomed in near a spiral center
juliaSet(-0.408, 0.595, -0.55:0.001:0.25, 0.2:0.001:0.5, 400, :acton)


# The Burning Ship
burningShip(-2:0.01:1, -0.5:0.01:2, 200)

# A smaller, much more exciting, burning ship within the burning ship
burningShip(-1.8:0.0002:-1.7, -0.01:0.0002:0.1, 200, :acton)


3
# A field of swords within the burning ship. Liked it so much I 
# made (a higher-res version of) it my desktop background
burningShip(-1.8:0.0001:-1.67, -0.01:0.00004:0.03, 200, :PuBu_3)

# A chaotic and messy burning ship julia set
bJuliaSet(-1, 0.8, -2:0.004:2, -1.1:0.004:1.1,  200, :bone_1)