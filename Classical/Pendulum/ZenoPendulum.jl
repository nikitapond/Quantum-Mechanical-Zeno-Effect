### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ 699e3322-63cd-11eb-33b2-2d5cbc8c23e8
begin
	using Plots
	using PlutoUI
end

# ╔═╡ a8adbf42-63cd-11eb-39e5-0f70827c32c0
mutable struct Pendulum
	l; theta; mass; g;
	
	omega;
	
	
	function Pendulum(l, theta, mass)
		return new(l,theta,mass,-9.81, 0)
	end
	
end

# ╔═╡ e57f67d0-63cd-11eb-33cd-0fb1fb14894f
function Update(pend::Pendulum, dt, omega_dist=0)
	
	omega_dot = pend.mass * pend.g * sin(pend.theta)
	
	pend.omega += omega_dot * dt + omega_dist
	
	pend.theta += pend.omega  * dt

end

# ╔═╡ 8a40d79e-63cf-11eb-03ef-fb2bde23e405
function ToCartesian(pend::Pendulum)
	x = pend.l * sin(pend.theta)
	y = -pend.l * cos(pend.theta)
	
	return [x,y]
end

# ╔═╡ 57cc3cdc-63d1-11eb-2a34-f9534d79461a
let
	pend = Pendulum(2, pi/2, 1)
	
	t_steps = 20000
	dt = 0.0001
	
	positions = zeros(t_steps, 2)
	
	
	
	
	for t in 1:t_steps
		xy = ToCartesian(pend)
		positions[t,:] = xy	
		Update(pend, dt)
	end
	result = positions
	
	anim = Animation()
	plot(xlim = (-2,2), ylim=(-2,0))
	
	t_jump = 200
	
	for t in 1:t_jump:t_steps
		plot((result[t,1], result[t,2]), seriestype = :scatter)
		
		trail = t>t_jump ? result[t-t_jump:t,:] : result[1:t,:]
		
		plot!(trail[:,1], trail[:,2])
		
		plot!(xlim = (-2,2), ylim=(-2,0))

		frame(anim)
	end
	gif(anim, "pendulum.gif", fps=200)
end

# ╔═╡ 6b3aede0-63d3-11eb-2f88-29c00081012b
function CalculateAverageAmplitute(thetas, omegas)
	
	turning_points = []
	
	for t in 2:size(thetas)[1]
		
		#if the product of angular velocities is less than 0, this means we have a sign change.
		if(omegas[t-1] * omegas[t] < 0)
			append!(turning_points, 0.5 * (thetas[t-1] + thetas[t]))
		end
		
	end
	
	sum = 0
	for t in turning_points
		sum += abs(t)
	end
	sum /= size(turning_points)[1]
	return sum

end

# ╔═╡ 195ab17e-63d4-11eb-000d-55e78e07fe33
let
	pend = Pendulum(2, pi/2, 1)
	
	t_steps = 20000
	dt = 0.0001
	
	positions = zeros(t_steps, 2)
	thetas = zeros(t_steps)
	omegas = zeros(t_steps)
	
	
	
	for t in 1:t_steps
		xy = ToCartesian(pend)
		positions[t,:] = xy	
		thetas[t] = pend.theta
		omegas[t] = pend.omega
		Update(pend, dt)
	end
	amd = CalculateAverageAmplitute(thetas, omegas) #pi/2, as we expect
end

# ╔═╡ c2f0c23c-63d4-11eb-138a-01b47769ab9b
function randOmega(max_)
	return (rand()-0.5) * 2 * max_
end

# ╔═╡ a2dca45c-63d4-11eb-0ca8-93520431a664
let
	
	
	pend = Pendulum(2, pi/2, 1)
	
	t_steps = 20000
	dt = 0.0001
	
	results = zeros(t_steps, 2)
	thetas = zeros(t_steps)
	omegas = zeros(t_steps)
	
	dist_freq = 100
	max_dist = .1
	for t in 1:t_steps
		
		
		
		xy = ToCartesian(pend)
		results[t,:] = xy	
		thetas[t] = pend.theta
		omegas[t] = pend.omega
		
		if(t % dist_freq == 0)
			Update(pend, dt, randOmega(max_dist))
		else
			Update(pend, dt)
		end
	end
	
		
		
	
	anim = Animation()
	plot(xlim = (-2,2), ylim=(-2,2))
	
	t_jump = 500
	
	for t in 1:t_jump:t_steps
		plot((results[t,1], results[t,2]), seriestype = :scatter)
		
		trail = t>t_jump ? results[t-t_jump:t,:] : results[1:t,:]
		
		plot!(trail[:,1], trail[:,2])
		
		plot!(xlim = (-2,2), ylim=(-2,2))

		frame(anim)
	end
	#gif(anim, "pendulum.gif", fps=200)
	
	amd = CalculateAverageAmplitute(thetas, omegas) #pi/2, as we expect
	gif(anim, "pendulum.gif", fps=200)

	
end

# ╔═╡ 7feb2d5c-63d5-11eb-053b-f10737c0aa8b
begin
	
	t_steps = 200000
	dt = 0.0001
	
	dist_freq = [100, 500, 1000, 2000, t_steps+1]
	max_dist = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
	
	average_amps = zeros(size(dist_freq)[1],size(max_dist)[1])
	
	iterations = 100
	
	
	Threads.@threads for i in 1:size(dist_freq)[1]
		freq = dist_freq[i]
		for (j, dist) in enumerate(max_dist)
			for it in 1:iterations
				pend = Pendulum(2, pi/2, 1)

				#results = zeros(t_steps, 2)
				thetas = zeros(t_steps)
				omegas = zeros(t_steps)
				for t in 1:t_steps

					thetas[t] = pend.theta % (2*pi)
					omegas[t] = pend.omega

					if(t % freq == 0)
						Update(pend, dt, randOmega(dist))
					else
						Update(pend, dt)
					end
				end
				average_amps[i,j] += CalculateAverageAmplitute(thetas, omegas)
			end
		end
	end
	average_amps./= iterations
	
end

# ╔═╡ 470b3c7c-63d6-11eb-04a5-438e39f33789
begin
	hm = heatmap(average_amps)
	plot(hm)
	yticks = ["100", "500", "1000", "2000", "no measure"]
	
	plot!(yticks =(1:5, yticks))
	plot!(xticks = (1:6, max_dist))
end

# ╔═╡ f4fbe2a8-63dd-11eb-2dea-4d26f56e3a66
average_amps[2, 5]

# ╔═╡ 7d60bc56-63dc-11eb-06c3-cf9a95fa71cd
begin
	plot(average_amps[3,:], xticks = max_dist)
	
	
end

# ╔═╡ 871b15e6-63de-11eb-0043-b38c69b0b747
begin

	
	max_dist_1000 = Array(0.05:0.05:0.4)
	
	average_amps_1000 = zeros(size(max_dist_1000)[1])
	iterations_1000 = 100
	freq = 1000
	Threads.@threads for i in 1:size(max_dist_1000)[1]
			
		dist = max_dist_1000[i]
		for it in 1:iterations_1000
			pend = Pendulum(2, pi/2, 1)

			#results = zeros(t_steps, 2)
			thetas = zeros(t_steps)
			omegas = zeros(t_steps)
			for t in 1:t_steps

				thetas[t] = pend.theta % (2*pi)
				omegas[t] = pend.omega

				if(t % freq == 0)
					Update(pend, dt, randOmega(dist))
				else
					Update(pend, dt)
				end
			end
			average_amps_1000[i] += CalculateAverageAmplitute(thetas, omegas)
		end
	end
	average_amps_1000./= iterations
	
end

# ╔═╡ ce16fbf8-63df-11eb-3e69-99ffa418f133
plot(average_amps_1000, xticks=(1:size(max_dist_1000)[1], max_dist_1000))

# ╔═╡ Cell order:
# ╠═699e3322-63cd-11eb-33b2-2d5cbc8c23e8
# ╠═a8adbf42-63cd-11eb-39e5-0f70827c32c0
# ╠═e57f67d0-63cd-11eb-33cd-0fb1fb14894f
# ╠═8a40d79e-63cf-11eb-03ef-fb2bde23e405
# ╠═57cc3cdc-63d1-11eb-2a34-f9534d79461a
# ╠═6b3aede0-63d3-11eb-2f88-29c00081012b
# ╠═195ab17e-63d4-11eb-000d-55e78e07fe33
# ╠═c2f0c23c-63d4-11eb-138a-01b47769ab9b
# ╠═a2dca45c-63d4-11eb-0ca8-93520431a664
# ╠═7feb2d5c-63d5-11eb-053b-f10737c0aa8b
# ╠═470b3c7c-63d6-11eb-04a5-438e39f33789
# ╠═f4fbe2a8-63dd-11eb-2dea-4d26f56e3a66
# ╠═7d60bc56-63dc-11eb-06c3-cf9a95fa71cd
# ╠═871b15e6-63de-11eb-0043-b38c69b0b747
# ╠═ce16fbf8-63df-11eb-3e69-99ffa418f133
