### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ fd20825c-6ba0-11eb-29ae-7d2924873c44
begin
	using Plots
	using FFTW
	using PlutoUI
	using NumericalIntegration
end

# ╔═╡ 17c4c1c2-6ba1-11eb-0256-61adfd9db8b8
begin
	function make_psi_x(psi_x, k, x, dx)
		return (psi_x .* exp.(-1im * k[1] .*x) * sqrt(2 * pi) / dx)
	end
	function get_psi_x(psi_mod_x, k, x, dx)
		return psi_mod_x .*exp.(1im * k[1] .*x)* sqrt(2 * pi) / dx
	end
	function make_psi_k(psi_k, x, dk, N)
		return psi_k .* exp.(1im * x[1] .*dk .* Array(0:N-1) )
	end
	function get_psi_k(psi_mod_k, x, dk, N)
		return psi_mod_k .* exp.(-1im * x[1] .* dk * Array(0:N-1))
	end
	function compute_k_from_x(psi_mod_x)
		return fft(psi_mod_x)
	end
	function compute_x_from_k(psi_mod_k)
		return ifft(psi_mod_k)
	end
end

# ╔═╡ 5023f006-6ba1-11eb-0654-d503a0de2782
mutable struct Simulation
	x # Float array
	v_x #Float array
	psi_x #Complex array

	hbar #Float
	m #Float
	t0 #Float
	k #Float array

	dx
	dk

	dt

	x_evolve_half
	x_evolve
	k_evolve

	psi_k 
	psi_mod_x
	psi_mod_k

	function Simulation(x, psi_x0, V_x, k0 = -28, hbar = 1., m = 1., t0 = 0., dt=0.01)

		N = size(x)[1]

		dx = x[2] - x[1];
		dk = 2 * pi / (N * dx);

		if(k0 == nothing)
			k0 = - 0.5 * N * dk
		end

		k =range(1, N, length=N).*dk .+ k0
		x_evolve_half = exp.(-0.5im .* V_x ./ hbar .* dt)
		x_evolve = x_evolve_half .* x_evolve_half
		k_evolve =  exp.(-0.5im .* hbar ./ m .* (k.*k) .* dt)

		psi_mod_x = make_psi_x(psi_x0, k, x,dx)

		psi_k = compute_k_from_x(psi_mod_x)

		psi_mod_k = make_psi_k(psi_k, x, dk,N)

		new(x,V_x, psi_x0, hbar, m,t0, k, dx, dk, dt, x_evolve_half, x_evolve, k_evolve, psi_k, psi_mod_x, psi_mod_k)

	#return sim
end

end

# ╔═╡ 5598e85c-6ba1-11eb-0bc6-a938d9f2e716
function sim_iteration(sim::Simulation)
	sim.psi_mod_x .*= sim.x_evolve_half

	sim.psi_mod_k = compute_k_from_x(sim.psi_mod_x)
	sim.psi_mod_k .*= sim.k_evolve;
	sim.psi_mod_x = compute_x_from_k(sim.psi_mod_k)
	sim.psi_mod_x .*= sim.x_evolve

	sim.psi_mod_x .*= sim.x_evolve_half
	#return get_psi_x(sim.psi_mod_x, sim.k, sim.x, sim.dx)
	return sim.psi_mod_x
end

# ╔═╡ a5ba68cc-6ba3-11eb-12a6-c9c9b223aeee
function to_real_normal(x, psi)
		real = abs.(psi .* conj(psi))
		normal = integrate(x, real)
		return real ./ normal
	end
	

# ╔═╡ caaf137e-6ba2-11eb-2221-0fafb4c0db80
"""
Contains global characteristics that will be shared by all experiments 
"""
struct GlobalArgs
	N
	dx
	dt
	
	hbar
	mass
	
	vel_0
	
	t_steps
	
end

# ╔═╡ 5b990e74-6ba1-11eb-0245-d1ed590e4631
function gauss_x(x, a, x0, k0)
	
	exp_part = (-0.5 * (((x .- x0) /a) .* ((x .- x0) /a)) + 1im .*x .*k0)
	const_part = (a * sqrt(pi)) ^ (-0.5)
	
	return const_part * exp.(exp_part)
	#return ((a * sqrt(pi)) ^ (-0.5) * arr_exp(-0.5 * ((x .- x0) /a) ^ 2 + 1im .*x .*k0))
end

	

# ╔═╡ 6cb6b014-6ba1-11eb-0a39-97628294a03e
function createDoubleWell(n, wellWidth, v0, vWell, vBoundry=9999999)
	
	V = fill(v0, (n))
	
	half = floor(Int64, n/2)
	quart = floor(Int64,half/2)
	
	
	halfWidth = floor(Int64,wellWidth/2)
	
	for i in 1:wellWidth
		
		V[quart-halfWidth + i] = vWell
		V[half + quart - halfWidth + i] = vWell
	end
	
	
	for i in 1:quart
		V[i] = vBoundry
		V[n-i + 1] = vBoundry
	end
	
	return V, quart
end

# ╔═╡ 58d14f20-6ba2-11eb-3ccb-9368d3dda110
@bind run CheckBox()

# ╔═╡ 6e8af61a-6ba2-11eb-149e-03aad5e4afbf
if(run)
	gArgs = GlobalArgs(1000, 0.1, 0.1, 1, 1, 10, 10000)
	x = (Array(1:gArgs.N) .- 0.5 * gArgs.N) .* gArgs.dx
	Vwell, x0 = createDoubleWell(gArgs.N, 10, 50, -10, 1000000000)
	
	x0 = x[x0]
	
	
	L = gArgs.hbar / sqrt(2*gArgs.mass*gArgs.vel_0)
	a = 3 * L
	
	p0= sqrt(2 * gArgs.mass * 0.2 * gArgs.vel_0)
	dp2 = p0 * p0 /80.0
	d = gArgs.hbar /sqrt(2*dp2)
	k0 = p0/gArgs.hbar
	v0 = p0/gArgs.mass
	
	psi_x0 = gauss_x(x,d, x0, k0)

	sim = Simulation(x, psi_x0, Vwell)
	psi_x_results = []

	
	for i in 1:gArgs.t_steps
		
		res = to_real_normal(x, sim_iteration(sim))
		append!(psi_x_results,res)
		
	end
	
end

# ╔═╡ 42059cba-6ba6-11eb-11df-07631e950957
@bind plot CheckBox()

# ╔═╡ b1de7dd4-6ba1-11eb-27a6-5109c282876b
begin
	V,  = createDoubleWell(1000, 100, 0, -200)
	plot(V, ylim=(-200, 10))
end

# ╔═╡ afbec052-6ba3-11eb-297c-1fc3a09f99db
if(plot)
	
	anim = Animation()
	N = gArgs.N
		#p = plot(x, v_x, size=(600,300))
		#plot!([0], [results[1:N])
	for t in 1:10:gArgs.t_steps-1

		
		plot(x, psi_x_results[(t*N +1):(t+1)*N])
		plot!(x, Vwell)
		plot(ylim = (-10, 50))
		frame(anim)
	end
	gif(anim, "test_1.gif", fps=200)
	
end

# ╔═╡ Cell order:
# ╠═fd20825c-6ba0-11eb-29ae-7d2924873c44
# ╠═17c4c1c2-6ba1-11eb-0256-61adfd9db8b8
# ╠═5023f006-6ba1-11eb-0654-d503a0de2782
# ╠═5598e85c-6ba1-11eb-0bc6-a938d9f2e716
# ╠═a5ba68cc-6ba3-11eb-12a6-c9c9b223aeee
# ╠═caaf137e-6ba2-11eb-2221-0fafb4c0db80
# ╠═5b990e74-6ba1-11eb-0245-d1ed590e4631
# ╠═6cb6b014-6ba1-11eb-0a39-97628294a03e
# ╠═b1de7dd4-6ba1-11eb-27a6-5109c282876b
# ╠═58d14f20-6ba2-11eb-3ccb-9368d3dda110
# ╠═6e8af61a-6ba2-11eb-149e-03aad5e4afbf
# ╠═42059cba-6ba6-11eb-11df-07631e950957
# ╠═afbec052-6ba3-11eb-297c-1fc3a09f99db
