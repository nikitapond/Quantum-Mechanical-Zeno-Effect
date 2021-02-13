### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ f279a816-6de3-11eb-0444-f7c49853dbbe
begin
	using Plots; pyplot();
	using FFTW; 
	using PlutoUI;
	using NumericalIntegration
	using Statistics
end

# ╔═╡ e2240790-6de3-11eb-36e0-9f1d39807148
md"""
Our aim in this notebook, is to use the FFT method to solve the Schrodinger equation for an electron on a 2d lattice.

The Schoringer equation in 2 dimensions can be written as:

$$(1) i\hbar \frac{\partial\psi}{\partial t} = -\frac{\hbar}{2m} \left[\frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2}\right]\psi + V \psi $$ 

To solve this, we use the Split-step fourier method

$$(2) \tilde{\psi}(u,v,t) = \frac{1}{\sqrt{2\pi}} \int\int_{-\infty}^{\infty}\psi(x,y,t)e^{-i(ux + vy)} dx dy$$

We can write our normal wave function in terms of its inverse fourier transform:

$$(3) \psi(x,y,t) = \frac{1}{\sqrt{2\pi}} \int\int_{-\infty}^{\infty}\psi(x,y,t)e^{-i(ux + vy)}$$

We can then substitute (1) into (3)


$$(4) i\hbar \frac{\partial \tilde{\psi}}{\partial t} = \frac{\hbar^2}{2m}(u^2 + v^2) \tilde{\psi} + V\left(i\frac{\partial}{\partial u}, i\frac{\partial}{\partial v}\right)\tilde{\psi}$$

We see that in fourier space, the first part relating to the Hamiltonian is trivial to solve:

$$(5)    i\hbar \frac{\partial \tilde{\psi}}{\partial t} = \frac{\hbar^2}{2m}(u^2 + v^2) \tilde{\psi}$$

$$(6) \tilde{\psi}(u,v,t+\Delta t) = \tilde{\psi}(u,v,t)e^{-i\hbar(u^2 + v^2)\Delta t/2m}$$

Similarly, our equation relating to the potential can be easily solved in the normal space

$$(7)i\hbar \frac{\partial\psi}{\partial t} = V(x,y) \psi$$
$$(8)\psi(x,y,t+\Delta t) = \psi(x,y,t) e^{-iV(x,y) \Delta t/\hbar}$$

From this it is clear, a numerical solution for the schrodinger equation can be found by solving this equation as follows:


"""

# ╔═╡ be5db710-6dee-11eb-1577-d37790b7d137


# ╔═╡ 51937f0c-6de4-11eb-30ce-e7397fef82c9
function MakePotential(total_N, grid_N, dx, a, V0)

	
			
	grid_count = Int(total_N/ grid_N)
	
	charges = []
	
	for i in 1:grid_count
		for j in 1:grid_count
			push!(charges, ((i-0.5)  * grid_N, (j-0.5) * grid_N))
		end
	end
    pot = fill(V0,(total_N, total_N)) * 1.0
	
    #return pot
	for g in 1:(grid_count^2)
		
		for x in 1:total_N
			for y in 1:total_N

				x_d = (x - charges[g][1])
				x_d -= total_N * round(x_d/total_N)
				
				y_d = (y - charges[g][2])
				y_d -= total_N * round(y_d/total_N)
				r_sqr =  a * (x_d^2 + y_d^2)

				pot[x,y] -= 1/(r_sqr+1)
			end
		end
	end
	

    
    return pot
end

# ╔═╡ 53958f70-6de4-11eb-2264-0b67a45869ba

struct Lattice
    grid_N # Size of each grid point 
    total_N #Total size of lattice
    dx;
    x; y;
    
    Pot
    Psi
    
    function Lattice(grid_size, dx, a, V0)
        grid_N = Int(1/dx)
        total_N = grid_N * grid_size
        x = Array(1:total_N) .* dx
        y = Array(1:total_N) .* dx
        Pot = MakePotential(total_N, grid_N, dx, a, V0)
        return new(grid_N,total_N,dx,x,y,Pot, nothing)
    end
    
end

# ╔═╡ 58db724c-6de4-11eb-0037-d31ec876b8c6
function PlotLatticePotential(lat::Lattice)
	plot(lat.x, lat.y, lat.Pot, st=:surface)
end

# ╔═╡ 0ba6b1c8-6de5-11eb-239c-cf607119a69b
begin
	lat = Lattice(3, 0.01, .001, 2)
	PlotLatticePotential(lat)
end

# ╔═╡ 76c14ae8-6df3-11eb-06f0-5b4094088ea2
function normalise(psi_xy, dx)
	
	N = size(psi_xy)[1]
	sum = 0
	for x in 1:N
		for y in 1:N
			sum += psi_xy[x,y] * dx^2
		end
	end
	psi_xy ./= sum
	return sum
end

# ╔═╡ 5e7184ba-6df0-11eb-197d-fde81eeabf5e
function gauss_xy(N, psi_xy, a, x0, y0, u0, v0, dx)
	
	d = 1/N
	
	x0 *=dx
	y0 *= dx
	
	for x_ in 1:N
		for y_ in 1:N
			x = x_ * dx
			y = y_ * dx
			psi_xy[x_,y_] = sqrt(a * sqrt(pi)) .* exp(-((x - x0)^2 + (y - y0)^2)/(2*a^2) + 1im*(x*u0 + y*v0)) 
			
		end
	end
	normalise(psi_xy, dx)
end

# ╔═╡ 0ee479fc-6df0-11eb-0b47-91b8aeea4c81
mutable struct System
	lat::Lattice
	psi_xy; psi_uv;
	
	xy_evolve_half; xy_evolove
	uv_evolve;
	function System(lattice, dt)
		N = lattice.total_N
		N_half = Int(N/2)
		psi_xy = zeros(N,N)
		
		gauss_xy(lattice.total_N, psi_xy, .1, N_half, N_half, 0,0, lattice.dx)
		return new(lattice, psi_xy)
	end
	
end

# ╔═╡ 812d7c3c-6dfa-11eb-25c3-6b2f9052ddaf
function PlotSystem(sys::System)
	#probability function
	psi_xy = real(sys.psi_xy .* conj.(sys.psi_xy))
	PlotLatticePotential(sys.lat)
	psi_xy_max = max(psi_xy ...)
	
	pot_mean = mean(sys.lat.Pot)
	
	scaled_psi_xy = (psi_xy.-0.1) ./ psi_xy_max
	
	plot!(sys.lat.x, sys.lat.y, scaled_psi_xy, st=:surface)
	plot!(zlim=(0, psi_xy_max))
end

# ╔═╡ 32265eca-6dfb-11eb-24ed-43fba7ff5395
begin
	system = System(lat, 0.01)
	PlotSystem(system)
end

# ╔═╡ 6f3106ec-6df3-11eb-3350-d7093af3b322
function fourier(sys::System)
	return 0
end

# ╔═╡ 24f40bba-6df3-11eb-1b58-ab9a60b14028
function Evolve(sys::System)
	sys.psi_xy .*= sys.xy_evolve_half
	sys.psi_uv = 0
end

# ╔═╡ 2d914e24-6df1-11eb-1d38-45bce1c2459c
begin
	N = 100
	mid = Int(N/2)
	psi_xy = zeros(N, N)
	gauss_xy(N, psi_xy, .1, mid, mid, 0,0, 0.001)
	normalise(psi_xy, 1/N)
	psi_xy = real.(psi_xy .* conj(psi_xy))
	plot(psi_xy, st=:surface)

end

# ╔═╡ ac206678-6df2-11eb-10d6-33d640fa0e85
i = normalise(psi_xy, 0.01)

# ╔═╡ Cell order:
# ╠═e2240790-6de3-11eb-36e0-9f1d39807148
# ╠═be5db710-6dee-11eb-1577-d37790b7d137
# ╠═f279a816-6de3-11eb-0444-f7c49853dbbe
# ╠═51937f0c-6de4-11eb-30ce-e7397fef82c9
# ╠═53958f70-6de4-11eb-2264-0b67a45869ba
# ╠═58db724c-6de4-11eb-0037-d31ec876b8c6
# ╠═0ba6b1c8-6de5-11eb-239c-cf607119a69b
# ╠═5e7184ba-6df0-11eb-197d-fde81eeabf5e
# ╠═0ee479fc-6df0-11eb-0b47-91b8aeea4c81
# ╠═812d7c3c-6dfa-11eb-25c3-6b2f9052ddaf
# ╠═32265eca-6dfb-11eb-24ed-43fba7ff5395
# ╠═6f3106ec-6df3-11eb-3350-d7093af3b322
# ╠═76c14ae8-6df3-11eb-06f0-5b4094088ea2
# ╠═24f40bba-6df3-11eb-1b58-ab9a60b14028
# ╠═2d914e24-6df1-11eb-1d38-45bce1c2459c
# ╠═ac206678-6df2-11eb-10d6-33d640fa0e85
