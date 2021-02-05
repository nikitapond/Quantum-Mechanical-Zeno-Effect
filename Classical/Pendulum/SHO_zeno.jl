### A Pluto.jl notebook ###
# v0.12.20

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

# ╔═╡ dd3575ba-66d0-11eb-36e8-add1165c2771
begin
	using Plots
	using PlutoUI
end

# ╔═╡ eb46ca7a-66d0-11eb-176c-47e881b4c84d
md"""
We have a SHO whos equation of motion takes the form:
$$m\ddot{x} = - kx$$
"""

# ╔═╡ d9abdb18-66d3-11eb-14d1-73230eb99ca8
function randDisturbance(distMax)
	return (rand() - 0.5) * 2 * distMax
end

# ╔═╡ 1e76464e-66d3-11eb-1063-b94c153373d5
function calculateAmplitute(results)
	turning_points = []
	start_ = Int(size(results)[1] * 0.25)
	for t in start_:size(results)[1]
		
		#if the product of velocities is less than 0, this means we have a sign change.
		if(results[t-1, 2] * results[t, 2] < 0)
			append!(turning_points, 0.5 * (results[t-1, 1] + results[t, 1]))
		end
		
	end
	
	sum = 0
	for t in turning_points
		sum += abs(t)
	end
	sum /= size(turning_points)[1]
	return sum
end

# ╔═╡ a6e77172-66d5-11eb-028a-f35ea97ff8db
mutable struct DataStore
	results;
	minX; maxX; avgAmp;
	function DataStore(results)
		minX = min(results...)
		maxX = max(results...)
		avgAmp = calculateAmplitute(results)
		return new(results, minX, maxX, avgAmp)
	end
end

# ╔═╡ 13e31524-66d1-11eb-35d6-59f4a54822c9
function runSHO(k, m, x0, tSteps, dt, measureFreq=-1, measureMax=0)
	
	results = zeros(tSteps,2 )
	if(measureFreq == -1)
		measureFreq = tSteps
	end
	v = 0
	x = x0
	for t in 1:tSteps
		acc = (- k * x) / m # ma = -kx
		v += acc * dt
		
		
		if(t % measureFreq  == 0)
			v += randDisturbance(measureMax)
		end
		
		x += v * dt
		
		results[t,1] = x
		results[t,2] = v
	end
	return DataStore(results)
end

# ╔═╡ b6db2302-66d1-11eb-36ef-4f41bac392bf
let
	
	resExp1 = runSHO(1, 1, 5, 10000, 0.001)
	
	
	animExp1 = Animation()
	plot(title="test")
	plot!(xlim = (resExp1.minX,resExp1.maxX), ylim=(-1,1))
	let
		t_jump = 100
		for t in 1:t_jump:size(resExp1.results)[1]
			plot((resExp1.results[t], 0), seriestype = :scatter)
			plot!(xlim = (resExp1.minX,resExp1.maxX), ylim=(-1,1))

			trail = t>t_jump ? resExp1.results[t-t_jump:t, 1] : resExp1.results[1:t, 1]
			trail_y = zeros(size(trail)[1])
			plot!(trail, trail_y)

			frame(animExp1)
		end
	end
	gif(animExp1, "sho.gif", fps=200)
end

# ╔═╡ a4719eb8-66d3-11eb-0c65-b54caffb263c
begin
	
	resExp2 = runSHO(1, 1, 5, 100000, 0.001, 1000, 2)
	
	
	animExp2 = Animation()
	plot(title="test")
	plot!(xlim = (resExp2.minX,resExp2.maxX), ylim=(-1,1))

	t_jump = 100
	for t in 1:t_jump:size(resExp2.results)[1]
		plot((resExp2.results[t], 0), seriestype = :scatter)
		plot!(xlim = (resExp2.minX,resExp2.maxX), ylim=(-1,1))
		
		trail = t>t_jump ? resExp2.results[t-t_jump:t, 1] : resExp2.results[1:t, 1]
		trail_y = zeros(size(trail)[1])
		plot!(trail, trail_y)
		
		frame(animExp2)
	end
	gif(animExp2, "sho_measure_test.gif", fps=200)
end

# ╔═╡ 4b2d0d22-66d4-11eb-2ddc-df63075dd642
resExp2.avgAmp

# ╔═╡ 5672d938-66d6-11eb-2541-f7124c3065bd
md"""
We see that by adding a regular disturbance, we increase the amplitute of the system
"""

# ╔═╡ a518f330-66d9-11eb-045c-057641b36c26
@bind testDiffFreq CheckBox()

# ╔═╡ 4d6b1fc6-66d6-11eb-3c5d-bd81f4c8b301
if(testDiffFreq)
	
	let
		total_steps = 640000 #640 thousand steps each
		measure_freq = [-1,6400,3200,1600,800,400, 200, 100]
		distMax = 0.5

		iterations = 50

		averageAmps = zeros(size(measure_freq)[1])
		#we iterate all our measurement frequecies
		for (i, mFreq) in enumerate(measure_freq)
			Threads.@threads for it in 1:iterations
				resExp_constDistMax = runSHO(1,1, 5, total_steps, 0.001, mFreq, distMax)
				averageAmps[i] += resExp_constDistMax.avgAmp
			end
		end

		averageAmps ./=iterations

		
		x = Array(1:size(measure_freq)[1])
		p = plot(x, averageAmps)
		plot!(xticks = (x, measure_freq))

		p
		
	end
end

# ╔═╡ d61fa410-66de-11eb-0180-d5189907779b
if(testDiffFreq)
	
	let
		total_steps = 128000 #128 thousand steps each
		measure_freq = [-1,6400,3200,1600,800,400, 200, 100]
		distMax = 0.5

		iterations = 2000

		averageAmps = zeros(size(measure_freq)[1])
		Threads.@threads for it in 1:iterations
		#we iterate all our measurement frequecies
			for (i, mFreq) in enumerate(measure_freq)
			
				resExp_constDistMax = runSHO(1,1, 5, total_steps, 0.001, mFreq, distMax)
				averageAmps[i] += resExp_constDistMax.avgAmp
			end
		end

		averageAmps ./=iterations

		
		x = Array(1:size(measure_freq)[1])
		p = plot(x, averageAmps)
		plot!(xticks = (x, measure_freq))

		p
		
	end
end

# ╔═╡ c31cf4ca-66da-11eb-147d-3ba9706317e8
@bind testDiffMax CheckBox()

# ╔═╡ c8bd49c0-66da-11eb-2c2d-a78dd78053f1
if(testDiffMax)
	
	let
		total_steps = 128000 #128 thousand steps each
		measure_freq = 800
		dist_maxes = Array(1:25) * 0.02 

		iterations = 2000

		averageAmps = zeros(size(dist_maxes)[1])
		#we iterate all our measurement frequecies
		for (i, distMax) in enumerate(dist_maxes)
			Threads.@threads for it in 1:iterations
				resExp_constDistMax = runSHO(1,1, 5, total_steps, 0.001, measure_freq, distMax)
				averageAmps[i] += resExp_constDistMax.avgAmp
			end
		end

		averageAmps ./=iterations

		let
			p = plot(dist_maxes, averageAmps)
			#plot!(xticks = (x, dist_maxes))

			p
		end
	end
end

# ╔═╡ 807a894c-66db-11eb-2495-479d8b75b6b9
begin
	test = zeros(100)
	for t in 1:100
		test[t] = randDisturbance(1)
	end
	test
end

# ╔═╡ Cell order:
# ╠═dd3575ba-66d0-11eb-36e8-add1165c2771
# ╠═eb46ca7a-66d0-11eb-176c-47e881b4c84d
# ╠═d9abdb18-66d3-11eb-14d1-73230eb99ca8
# ╠═1e76464e-66d3-11eb-1063-b94c153373d5
# ╠═a6e77172-66d5-11eb-028a-f35ea97ff8db
# ╠═13e31524-66d1-11eb-35d6-59f4a54822c9
# ╠═b6db2302-66d1-11eb-36ef-4f41bac392bf
# ╠═a4719eb8-66d3-11eb-0c65-b54caffb263c
# ╠═4b2d0d22-66d4-11eb-2ddc-df63075dd642
# ╠═5672d938-66d6-11eb-2541-f7124c3065bd
# ╠═a518f330-66d9-11eb-045c-057641b36c26
# ╠═4d6b1fc6-66d6-11eb-3c5d-bd81f4c8b301
# ╠═d61fa410-66de-11eb-0180-d5189907779b
# ╠═c31cf4ca-66da-11eb-147d-3ba9706317e8
# ╠═c8bd49c0-66da-11eb-2c2d-a78dd78053f1
# ╠═807a894c-66db-11eb-2495-479d8b75b6b9
