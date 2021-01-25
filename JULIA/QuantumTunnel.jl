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

# ╔═╡ 1b746900-5c13-11eb-131a-35dbeb769750
using Plots

# ╔═╡ 8c40b6f0-5ca6-11eb-166d-7bcc06b57824
using FFTW


# ╔═╡ 793e1e50-5ccc-11eb-1a48-970cca724a39
using InteractiveUtils


# ╔═╡ c745fe40-5cd3-11eb-0834-cdc44d353dd3
using NumericalIntegration

# ╔═╡ 31eac102-5cb0-11eb-29f1-2b5dbde52ad7
function arr_mult(arr1, arr2)
	
	result = nothing
	if(eltype(arr1) <: Complex || eltype(arr2) <: Complex)
		result =  zeros(size(arr1)) * 1im
	else
		result = zeros(size(arr1))	
	end
	
	for en in enumerate(arr1)
		
		result[en[1]] = en[2] * arr2[en[1]]
	
	end
	return result
end


# ╔═╡ 41d16512-5cb5-11eb-2a78-dfc1caae2aaf
function arr_exp(arr)
	result = nothing
	if(eltype(arr) <: Complex)
		result = zeros(size(arr)) * 1im
	else
		result = zeros(size(arr))
	end
	for en in enumerate(arr)
		result[en[1]] = exp(en[2])
	end
	return result
end

# ╔═╡ 824c3340-5cb5-11eb-2065-0f921693205d
function arange_n(N)
		return Array(0:N-1)
	end

# ╔═╡ 8e0ef730-5cba-11eb-10ad-19b4ccb48665
# function set_psi_x(sim::Simulation)
# 		sim.psi_mod_x = (sim.psi_x .* arr_exp(-1im * sim.k[1] .*sim.x) * sqrt(2 * pi) / sim.dx
# 		return
# 	end
function make_psi_x(psi_x, k, x, dx)
	return (psi_x .* arr_exp(-1im * k[1] .*x) * sqrt(2 * pi) / dx)
	end

# ╔═╡ bc888c70-5cba-11eb-3d3e-19bcf1a4aea1
# function get_psi_x(sim::Simulation)
# 		return sim.psi_mod_x .*arr_exp(1im * sim.k[1] .*sim.x)* sqrt(2 * pi) / sim.dx
# 	end
function get_psi_x(psi_mod_x, k, x, dx)
		return psi_mod_x .*arr_exp(1im * k[1] .*x)* sqrt(2 * pi) / dx
	end

# ╔═╡ e5727c90-5cba-11eb-13aa-2d0efa1b45de
# function set_psi_k(sim::Simulation)
# 		sim.psi_mod_k = (sim.psi_k .* arr_exp(1im * sim.x[1] .*sim.dk .* arange_n(sim.N) )
# 		return
# 	end
function make_psi_k(psi_k, x, dk, N)
		return psi_k .* arr_exp(1im * x[1] .*dk .* arange_n(N) )
		
	end

# ╔═╡ 196807e2-5cbb-11eb-389b-5becb60e5b3f
function get_psi_k(psi_mod_k, x, dk, N)
		return psi_mod_k .* arr_exp(-1im * x[1] .* dk * arange_n(N))
		
	end

# ╔═╡ 80f9b3e2-5cbb-11eb-025b-6dc128074fb8
begin
# 	function compute_k_from_x(sim::Simulation)
# 		sim.psi_mod_k = fft(sim.psi_mod_x)
# 		return
# 	end
	
# 	function compute_x_from_k(sim::Simulation)
# 		sim.psi_mod_x = ifft(sim.psi_mod_k)
# 	end
	function compute_k_from_x(psi_mod_x)
		return fft(psi_mod_x)
		
	end
	
	function compute_x_from_k(psi_mod_k)
		return ifft(psi_mod_k)
	end
end

# ╔═╡ a49b2bf0-5cb9-11eb-130c-e163102987dc
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
		x_evolve_half = arr_exp(-0.5im .* V_x ./ hbar .* dt)
		x_evolve = x_evolve_half .* x_evolve_half
		k_evolve = arr_exp(-0.5im .* hbar ./ m .* (k.*k) .* dt)

		psi_mod_x = make_psi_x(psi_x0, k, x,dx)

		psi_k = compute_k_from_x(psi_mod_x)

		psi_mod_k = make_psi_k(psi_k, x, dk,N)

		new(x,V_x, psi_x0, hbar, m,t0, k, dx, dk, dt, x_evolve_half, x_evolve, k_evolve, psi_k, psi_mod_x, psi_mod_k)

	#return sim
end

end

# ╔═╡ 7473bfa0-5cb9-11eb-3814-8762c9ee9bf8
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

# ╔═╡ 38ff7850-5cd8-11eb-305b-45fafce7dad0
	

# ╔═╡ 91ae28e2-5cb8-11eb-21f7-0f633f56f889
function gauss_x(x, a, x0, k0)
	
	exp_part = (-0.5 * (((x .- x0) /a) .* ((x .- x0) /a)) + 1im .*x .*k0)
	const_part = (a * sqrt(pi)) ^ (-0.5)
	
	return const_part * arr_exp(exp_part)
	#return ((a * sqrt(pi)) ^ (-0.5) * arr_exp(-0.5 * ((x .- x0) /a) ^ 2 + 1im .*x .*k0))
end

	

# ╔═╡ 2aa0f9e0-5cde-11eb-1d73-e33f91f5193d
function getMomentum(psi_x, x, a)
	const_part = (a * sqrt(pi)) ^ (-0.5)

	exp_part = psi_x ./ const_part
	
	image_part = convert(Array{Float64}, imag.(psi_x))
	
	sqrt_mod = sqrt.(image_part.^2) 
	
	k0_arr = log.(sqrt_mod)
	
	
	return k0_arr
end
	
	
	

# ╔═╡ ee256ca0-5cb8-11eb-3625-575f35c0b028
function integrate_till_p(x, fx, p)
	#println("start")
	dx = x[2]-x[1]
	sum = 0
	for i in 1:size(fx)[1]
		#println("$p, $sum")
		sum += fx[i]*dx
		if(sum >= p)
			x_val = x[i]
			#println("Random p was: $p, found at i=$i, with x=$x_val")
			return x[i]
		end
	end
	return Inf
end
	

# ╔═╡ 9affb2a0-5ca5-11eb-061e-655f4b276155
begin
	
	#function integrate(x, fx)
	
	function to_real_normal(x, psi)
		real = abs.(psi .* conj(psi))
		normal = integrate(x, real)
		return real ./ normal
	end
	
	function measure(sim::Simulation, d)
		#Convert our wave function to a probability dsitribution
		prob_x = to_real_normal(sim.x, sim.psi_mod_x)
		#do same for k-space
		psi_k = get_psi_k(sim.psi_mod_k, x, sim.dk, size(sim.x)[1])
		prob_k = to_real_normal(sim.k, psi_k)
		
		#generate a random value ∈ [0,1] , integrate wave function
		#up until point x0 such that integrand = p		
		p = rand()
		x0 = integrate_till_p(sim.x, prob_x, p)
		k0 =  integrate_till_p(sim.k, prob_k, p)
		
		#println(k0)
		#Find the index of this x coordinate
		
		
		psi_x0 = gauss_x(sim.x, d, x0, k0)
		
		sim.psi_mod_x = make_psi_x(psi_x0, sim.k, sim.x,sim.dx)
		sim.psi_k = compute_k_from_x(sim.psi_mod_x)
		sim.psi_mod_k = make_psi_k(sim.psi_k, sim.x, sim.dk,size(x)[1])
		
		return x0
		
	end
	
	function getXIndex(x, x_find)
		
		for i in 1:size(x)[1]-1
			if(x[i] < x_find && x[i+1] > x_find)
				return i
			end
		end
		return -1
	end
		
	
	N = 1000
	dx = 0.1
	x = (Array(1:N) .- 0.5 * N) .* dx
	
	hbar = 1
	mass = 1
	
	V0 = 40
	L = hbar / sqrt(2*mass*V0)
	a = 3 * L
	x0 = -60 * L
	v_x = zeros(N)
	
	start = convert(Int64, N/2 + N*0.05)
	width = convert(Int64, N * 0.01)
	#width = 1
	for i in start:start+width
		v_x[i] = 4
	end
	
	boundary = 15
	bound_v = 999999999999999999
	for i in 1:boundary
		v_x[i] =bound_v
		v_x[N+1-i] = bound_v
	end
	#v_x[1] = 9999999999
	#v_x[N] = 9999999999
	#v_x[1] = 99999
	#v_x[1200] = 20
	
	
	p0= sqrt(2 * mass * 0.2 * V0)
	dp2 = p0 * p0 /80.0
	d = hbar /sqrt(2*dp2)
	k0 = p0/hbar
	v0 = p0/mass

	
	psi_x0 = gauss_x(x,d, x0, k0)
	
	x0_index = getXIndex(x, x0)
	
	iter = 1500
	psi_x_results = []
	psi_k_results = []
	
	sim = Simulation(x, psi_x0, v_x)
	
	measure_every = 100
	for i in 1:iter
		#println(i)
		res = to_real_normal(x, sim_iteration(sim))
		if( i % measure_every == 0)
			measure(sim, d)
			
		end
		append!(psi_x_results,res)
		res_k_raw = get_psi_k(sim.psi_mod_k, x, sim.dk, N)
		res_k = to_real_normal(x, res_k_raw)
		append!(psi_k_results,res_k)

	end
	
	
	
# 	for i in 1:1000
# 		res =to_real(sim_iteration(sim))
# 		plot(x,res)
		
# 		println("start")
# 		println(res)
# 		println("end")
# 		break
# 	end
end
	


# ╔═╡ 63353e60-5d51-11eb-0ac9-773a6fe9a16e


# ╔═╡ 36d44c30-5cde-11eb-2f61-63f32207ecb2
plot(sim.k, abs.(sim.psi_k))

# ╔═╡ 415bd862-5ce0-11eb-3d63-0dc3de28cbb2
abs.(sim.psi_k)[x0_index]

# ╔═╡ 21484402-5ccc-11eb-0df3-fdafe0fdc18e
begin
	@bind should_run html"<input type='checkbox'>"
end

# ╔═╡ f6a74980-5cd5-11eb-1331-b5e38071fcb7
let
	if(should_run)
		
		
		potential = (v_x .* 0.1)[boundary+1:N-boundary-2]
		
		pot_x = x[boundary+1:N-boundary-2]
		
		anim = Animation()
		#p = plot(x, v_x, size=(600,300))
		#plot!([0], [results[1:N])
		for t in 0:10:iter-1
			
			plot([x,x, pot_x], [psi_x_results[(t*N +1):(t+1)*N],psi_k_results[(t*N +1):(t+1)*N], potential], layout=(2,1), size=(600,300))
			frame(anim)
		end
		gif(anim, "lol.gif", fps=200)
		#scatter!([0], [sin,cos])
		# for i in 0:0.1:π
		# 	p[3] = [x], [sin(i)]
		# 	p[4] = [i], [cos(i)]
		# 	frame(anim)
		# end
		#gif(anim)
	end
end

# ╔═╡ be5f68b0-5d51-11eb-1d00-b5752101ee89
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
	
	v_height
	v_start
	v_width
	
	outer_boundary_height
	out_boundary_width
	
	t_steps
	
end

# ╔═╡ 389f9d40-5d55-11eb-1c9d-01537d1548c0
struct AnimatableResult
	x
	v_x
	v_boundary
	N
	t_steps
	psi_x_res
	psi_k_res
end

# ╔═╡ 8bffe9e0-5d55-11eb-2865-77e638194262
function plotAnimationResult(anim::AnimatableResult)
	potential = (anim.v_x .* 0.1)[anim.v_boundary+1:anim.N-anim.v_boundary-2]
		
	pot_x = anim.x[anim.v_boundary+1:anim.N-anim.v_boundary-2]

	anim_ = Animation()
	#p = plot(x, v_x, size=(600,300))
	#plot!([0], [results[1:N])
	for t in 0:10:anim.t_steps-1

		plot([anim.x,anim.x, pot_x], [anim.psi_x_res[(t*anim.N +1):(t+1)*anim.N],anim.psi_k_res[(t*anim.N +1):(t+1)*anim.N], potential], layout=(2,1), size=(600,300))
		frame(anim_)
	end
	gif(anim_, "lol.gif", fps=200)
end

# ╔═╡ 600a20b0-5d52-11eb-1b39-ad226a76aea4
function RunSimulation(gArgs::GlobalArgs, measure_every, iterations, return_result=false, debug=false)
	
	
	x = (Array(1:gArgs.N) .- 0.5 * gArgs.N) .* gArgs.dx

	
	L = gArgs.hbar / sqrt(2*gArgs.mass*gArgs.vel_0)
	a = 3 * L
	x0 = -60 * L
	v_x = zeros(N)
	#width = 1
	for i in gArgs.v_start:gArgs.v_start+gArgs.v_width
		v_x[i] = gArgs.v_height
	end
	
	for i in 1:gArgs.out_boundary_width
		v_x[i] =gArgs.outer_boundary_height
		v_x[gArgs.N+1-i] = gArgs.outer_boundary_height
	end
	
	
	v_x_start = x[ gArgs.v_start]
	v_x_end = x[gArgs.v_start + gArgs.v_width]
	
	
	p0= sqrt(2 * gArgs.mass * 0.2 * gArgs.vel_0)
	dp2 = p0 * p0 /80.0
	d = gArgs.hbar /sqrt(2*dp2)
	k0 = p0/gArgs.hbar
	v0 = p0/gArgs.mass
	
	#How many measurements we will take for each iteration
	measure_ratio = gArgs.t_steps / measure_every
	total_measurements = convert(Int64,ceil(measure_ratio))
	
	#We store our results as an array
	#Every time we measure our wave, we will note if our system is in the original state (undecayed = 0), or if it has passed the potential barrier (decayed = 1)
	#This v
	results = zeros(total_measurements)
	
	#only used if return_result = true
	psi_x_results = []
	psi_k_results = []
	
	iteration_10_perc = ceil(iterations * 0.1)
	
	#Run number of iterations
	for it in 1:iterations-1
		if(it % iteration_10_perc == 0)
			perc = round(it/iterations * 100, digits=2)
			println("$perc complete")
		end
		measure_num  = 0
		psi_x0 = gauss_x(x,d, x0, k0)
		#TODO - allow for variable dt
		sim = Simulation(x, psi_x0, v_x)

		for t in 1:gArgs.t_steps
			#run a single iteration
			res = sim_iteration(sim)
			
			#If we wish to return a result, we only do so for our first experiment
			if(return_result && it==1)
				append!(psi_x_results, to_real_normal(x, res))
				res_k_raw = get_psi_k(sim.psi_mod_k, x, sim.dk, N)
				res_k = to_real_normal(x, res_k_raw)
				append!(psi_k_results, res_k)
			end
			
			#If we wish to make a measurement this iteration
			if(t % measure_every == 0)
				measure_num += 1
				#Measure our quantum system, collect its new position
				wave_pos = measure(sim, d)
				if(wave_pos > v_x_end)
					results[measure_num] += 1
				end
			end
			
		end
		#Perform a final measurement
		wave_pos = measure(sim, d)
		
		if(wave_pos > v_x_end)
			if(debug)
				println("Quantum tunnel yes")
			end
			
			results[total_measurements] += 1
		elseif(debug)
			println("Quantum tunnel no")
		end
		
	end
	
	results ./= iterations
	
	println(size(psi_x_results))
	animResult = AnimatableResult(x, v_x, gArgs.out_boundary_width,gArgs.N, gArgs.t_steps, psi_x_results, psi_k_results)
	
	return animResult, results
end


# ╔═╡ 848d52f0-5d60-11eb-2fa3-93a45df6bff8
begin
	gArgs = GlobalArgs(1000, 0.1, 0.1, 1, 1, 40, 4.25, 550, 50, 99999999999, 15, 1250)
	iterations = 100
end

# ╔═╡ e35b2eb0-5d51-11eb-1e88-092528496ad9
begin
	#gArgs = GlobalArgs(1000, 0.1, 0.1, 1, 1, 40, 3.5, 550, 50, 99999999999, 15, 1500)
	
	 
	anim, results_no_measure = RunSimulation(gArgs,	25000, iterations, true, true)
	plotAnimationResult(anim)
	
end


# ╔═╡ f6bae632-5d51-11eb-004f-7bcad0682abe
results_no_measure # measured probability of tunneling occuring after all iterations

# ╔═╡ 31fc0100-5d5e-11eb-206c-cf1c1969b692
begin
	#gArgs = GlobalArgs(1000, 0.1, 0.1, 1, 1, 40, 3.5, 550, 50, 99999999999, 15, 1500)
	
	 
	anim_100, results_ever_100 = RunSimulation(gArgs,	100, iterations, true, true)
	plotAnimationResult(anim_100)
	
end


# ╔═╡ c45fbe10-5d5e-11eb-3683-e5901bba715f
begin
	result_100_x = Array(1:100:gArgs.t_steps)
	plot(result_100_x, results_ever_100)
end

# ╔═╡ Cell order:
# ╠═1b746900-5c13-11eb-131a-35dbeb769750
# ╠═8c40b6f0-5ca6-11eb-166d-7bcc06b57824
# ╠═793e1e50-5ccc-11eb-1a48-970cca724a39
# ╠═c745fe40-5cd3-11eb-0834-cdc44d353dd3
# ╠═31eac102-5cb0-11eb-29f1-2b5dbde52ad7
# ╠═41d16512-5cb5-11eb-2a78-dfc1caae2aaf
# ╠═824c3340-5cb5-11eb-2065-0f921693205d
# ╠═8e0ef730-5cba-11eb-10ad-19b4ccb48665
# ╠═bc888c70-5cba-11eb-3d3e-19bcf1a4aea1
# ╠═e5727c90-5cba-11eb-13aa-2d0efa1b45de
# ╠═196807e2-5cbb-11eb-389b-5becb60e5b3f
# ╠═80f9b3e2-5cbb-11eb-025b-6dc128074fb8
# ╠═a49b2bf0-5cb9-11eb-130c-e163102987dc
# ╠═7473bfa0-5cb9-11eb-3814-8762c9ee9bf8
# ╠═38ff7850-5cd8-11eb-305b-45fafce7dad0
# ╠═91ae28e2-5cb8-11eb-21f7-0f633f56f889
# ╠═2aa0f9e0-5cde-11eb-1d73-e33f91f5193d
# ╠═ee256ca0-5cb8-11eb-3625-575f35c0b028
# ╠═9affb2a0-5ca5-11eb-061e-655f4b276155
# ╠═63353e60-5d51-11eb-0ac9-773a6fe9a16e
# ╠═36d44c30-5cde-11eb-2f61-63f32207ecb2
# ╠═415bd862-5ce0-11eb-3d63-0dc3de28cbb2
# ╠═21484402-5ccc-11eb-0df3-fdafe0fdc18e
# ╠═f6a74980-5cd5-11eb-1331-b5e38071fcb7
# ╠═be5f68b0-5d51-11eb-1d00-b5752101ee89
# ╠═389f9d40-5d55-11eb-1c9d-01537d1548c0
# ╠═8bffe9e0-5d55-11eb-2865-77e638194262
# ╠═600a20b0-5d52-11eb-1b39-ad226a76aea4
# ╠═848d52f0-5d60-11eb-2fa3-93a45df6bff8
# ╠═e35b2eb0-5d51-11eb-1e88-092528496ad9
# ╠═f6bae632-5d51-11eb-004f-7bcad0682abe
# ╠═31fc0100-5d5e-11eb-206c-cf1c1969b692
# ╠═c45fbe10-5d5e-11eb-3683-e5901bba715f
