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

# ╔═╡ 5b4c4458-5f38-11eb-0998-db9dfa182660
begin
	using Plots;	pyplot()

	using LaTeXStrings
#	using Plots.PlotMeasures
	using PlutoUI
end

# ╔═╡ 6ef5c61e-5fb5-11eb-3c1a-f16d8288c036
md"""
# Classical Zeno Effect Seen in the Coherence of Light Sources
In this notebook, we wish to repeat the results of [1].
In this paper, the authors perform an example of the Zeno effect in the classical limit, showing that use of slits to cause a distrubance on a beam of light can cause the measured intensity over a set distance to increase.

To do this, we have 4 equations we need to find:\

1.    $J(x_1, x_2) = \langle E(x_1)E^*(x_2)\rangle = \frac{1}{2a} exp\left[\frac{(x_1-x_2)^2}{d^2}\right]$

2.    $P = \int_{-a}^a J(x,x) dx$

3.    $\mu_g = \frac{1}{P}\cdot\left[\int\int_{-a}^{a} |J(x_1,x_2)|^2\right]^{1/2}$

4.    $J(x_1',x_2') = \frac{k}{2\pi z} \cdot \int\int_{-a}^a J(x_1, x_2) \times K(x_1' - x_1)K^*(x_2'-x_2) dx_1 dx_2$

Where:\

$K(x) = exp\left[\frac{ikx^2}{2z}\right]$

We solve this set of equation numerically. 

[1] M.A. Porras, A. Luis, and I. Gonzalo, Classical Zeno dynamics in the light emitted by an extended, partially
coherent source, Phys. Rev. A 88, 052101
"""

# ╔═╡ cbd62ed8-5f3c-11eb-15bb-131bdf424007
"""
Struct containg all variables and data arrays needed for this calculation.
# Arguments 
- n::Int : Dimensionality of the system
- zmax::Float : Distance from source slit to detector
- nr::Int : Number of slits the beam is to pass through, excluding souce and detector slit.
- dd::Float - 'transverse coherence length is proportional to dd'
"""
mutable struct Data
	n; zmax; nr; dd;
	
	a; b; c; d; e;
	
	x;
	
	z; step; ss;
	
	function Data(n, zmax, nr, dd)
		
		a = complex(zeros(n,n))
		b = complex(zeros(n,n))
		c = complex(zeros(n,n))
		d = complex(zeros(n,n))
		e = complex(zeros(n,n))
		
		x = zeros(n) * 0.0
		
		z = zmax/nr
		#step seems to be our 'dx'
		step = 2/(n-1)
		
		#This is the only 
		ss = step * (1-1im)/(2*sqrt(pi * z))
		
		return new(n,zmax,nr,dd,a,b,c,d,e,x,z,step,ss)
	end
end

# ╔═╡ 31566570-5f3d-11eb-2bed-91e88dad9bfd
"""
Impliments equation (2)
"""
function calculatePower(dat::Data)
	power = 0
	
	for i in 1:dat.n-1
		power += dat.a[i,i] * dat.step
	end
	#This currently can return a complex number. Is this correct?
	return Float64(abs(power))
end

# ╔═╡ 2c65ab14-5f3d-11eb-3f8d-67f64d0bc471
"""
Impliments equation (3)
"""
function calculateCoherence(dat::Data, power)
	coh = 0;
	
	for i in 1:dat.n
		for j in 1:dat.n
			#
			coh += abs(dat.a[i,j])^2 * 2 * dat.step^2
		end
	end
	
	coh = Float64(sqrt(coh)/abs(power))
	return coh;
end

# ╔═╡ ac2673b0-5f3f-11eb-21c0-9d6bce434a99
function getCentreIntensity(dat::Data)
	return abs(dat.a[convert(Int64, (dat.n+1)/2), convert(Int64, (dat.n+1)/2)])
end

# ╔═╡ 5bfe1386-5f3d-11eb-03e8-75f1ba6e9db9
function initExperimentValues(dat::Data)
	for i in 1:dat.n
		dat.x[i] = -1 + 2 * (i-1)/(dat.n-1);
	end
	
	for i in 1:dat.n
		for j in 1:dat.n
			
			#dat.a = J(x1,x2) -> Equation 1 (a=1)
			dat.a[i,j] = exp(- (dat.x[i] - dat.x[j])^2 / dat.dd^2)/2
			
			#Equation for K(x) and K*(x)
			#We have also included the numberical 
			dat.b[i,j] = exp(-1im * (dat.x[i]-dat.x[j])^2 / (2*dat.z)) * dat.ss
			dat.d[i,j] = conj(dat.b[i,j])
		end
	end
end

# ╔═╡ c9463304-5f3a-11eb-30e6-e1882e01e062
"""
Runs the calculation, iterating over all slits within this experiment to
calculate the final result at the detector
"""
function iterateOverSlits(dat::Data, print_all=false)
	#iterate all slits
	for m in 1:dat.nr
		
		#Fill the first octant of e
		for i in 1:convert(Int64, (dat.n+1)/2)
			for j in 1:i
				dat.e[i,j] = 0
				for k in 1:dat.n
					dat.c[k,j] = 0
					for l in 1:dat.n
						dat.c[k,j ] += dat.a[k,l] * dat.b[l,j]
					end
					dat.e[i,j] += dat.d[i,k] * dat.c[k,j]
				end
			end
		end
		 
		#Fill the octant below
		for i in convert(Int64, (dat.n+3)/2):dat.n
			for j in 1:(dat.n-i+1)
				dat.e[i,j] = 0
				for k in 1:dat.n
					dat.c[k,j] = 0
					for l in 1:dat.n
						dat.c[k,j] += dat.a[k,l] * dat.b[l,j]
					end
					dat.e[i,j] += dat.d[i,k]* dat.c[k,j]
				end
			end
		end
		
		#Fill the final quadrants
		for i in convert(Int64, (dat.n+3)/2):dat.n
			for j in (dat.n-i+2):i
				dat.e[i,j] = conj(dat.e[dat.n+1-j, dat.n+1-i])
			end
		end
		
		#finalising data	
		for i in 1:dat.n
			for j in 1:dat.n
				if( j<= i)
					dat.a[i,j] = dat.e[i,j]
				else
					dat.a[i,j] = conj(dat.e[j,i])
				end
			end
		end
	end
	
	return
end


# ╔═╡ 28168b1a-5f3e-11eb-0f4a-493b9b0399e8
"""
Runs the experiment on the supplied data. We first initialise this data.
We then calculate the relevent start values, before running the calculation by calling 'iterateOverSlits'.
Once this is complete, we re-calculate relevent end values before returning all.
"""
function runExperiment(dat::Data)
	#Set up 
	initExperimentValues(dat)
	
	#calculate start values
 	startPower = calculatePower(dat)
 	startCoh = calculateCoherence(dat, startPower)
 	startIntensity = getCentreIntensity(dat)
	startA = copy(dat.a)
	
	#run calculation
	iterateOverSlits(dat)
	
	#calculate end values
 	endPower = calculatePower(dat)
	endCoh = calculateCoherence(dat, endPower)
 	endIntensity = getCentreIntensity(dat)
	endA = copy(dat.a)
	
	return startPower, startCoh, startIntensity, endPower, endCoh, endIntensity, startA, endA

end

# ╔═╡ 5a86cac8-5f3f-11eb-32c0-41ea6c56e911
"""
Runs a set of different experiments, with all varaibles constant except
the total number of slits 'nr'. 
# Arguments
- 'n::Integer' : The dimensionality to solve over
- 'zmax::Float' : The distance between soruce slit and detector
- 'dd::Float' : Not really sure
"""
function runMultipleExperiments(n, zmax, dd, min_nr, max_nr)
	
	results = zeros(max_nr-min_nr + 1, 7) * 1im
	#iterate for different slit counts
	for nr in min_nr:max_nr
		#create data with this nr
		dat = Data(n,zmax, nr, dd)
		
		#run calculation
		sPow, sCoh, sInt, ePow, eCoh, eInt = runExperiment(dat)
		
		#store all values
		results[nr - min_nr + 1, 1] = nr
 		results[nr - min_nr + 1, 2] = sPow;
 		results[nr - min_nr + 1, 3] = sCoh;
 		results[nr - min_nr + 1, 4] = sInt;
		
 		results[nr - min_nr + 1, 5] = ePow;
 		results[nr - min_nr + 1, 6] = eCoh;
 		results[nr - min_nr + 1, 7] = eInt;
		
	end
	return results
	
end

# ╔═╡ 33803354-607c-11eb-15ec-d71ffaa07000
begin
	struct GlobalArgs
		n; zmax; dd;
	end
end

# ╔═╡ 49a99602-607c-11eb-24da-6900e24db0ae
exp_1_args = GlobalArgs(51, 0.5, 0.1)

# ╔═╡ 5a461e62-607d-11eb-310a-c3fac838e1c9
min_max_nr_exp_1 = (1, 10)

# ╔═╡ bf547a4e-5f3b-11eb-17ed-f1685ba4801d
@bind run_exp_1 CheckBox()

# ╔═╡ 29cb5eb4-6159-11eb-1fba-3164718f1193
md"""
Press the toggle button above to start the calculation
"""

# ╔═╡ 313fa49a-5f40-11eb-3e48-a90c4320c965
if(run_exp_1)
	
	#n = 51, zmax = 0.5, dd=0.1, nr_min=1, nr_max=5
	results = runMultipleExperiments(exp_1_args.n, exp_1_args.zmax, exp_1_args.dd, min_max_nr_exp_1[1], min_max_nr_exp_1[2]);
	"Experiment 1 run succesfully for slit counts between $(min_max_nr_exp_1[1]) and $min_max_nr_exp_1[2]"
	
end

# ╔═╡ 8d995928-5fb8-11eb-28ec-e90a1a5bf8de
md"""
We have generated and stored our results as a 2d array. Each row represents a different experiment, with the collomns representing (from left to right):\

1.    NR - The number of slits for this experiment
2.    Start Power (Complex?)
3.    Start Coherence
4.    Start Intensity
5.    End Power (Complex)
6.    End Coherence
7.    End Intensity

We can now try and plot these, such that we can compare them to the original paper.
"""

# ╔═╡ 3b6f627e-6081-11eb-2b90-4feedc83d256
function saveFigToDir(dir, fname, fig)
	splitstr = split(dir, "/")

	for s in splitstr


		if(!isdir(s))
			mkdir(s)
		end
		cd(s)
	end
	savefig(fig, fname)
	for i in splitstr
		cd("..")
	end
end


# ╔═╡ 87a90dbc-5fb8-11eb-1ff7-15608ff9fb10
if(run_exp_1)	
	all_N = convert(Array{Int}, results[:,[1]]);
	
	#Take the absolute value of our power at the end of the experiment
	abs_end_power = abs.(results[:, [5]]);
	p_N_over_p_1 =  abs_end_power ./ abs_end_power[1];
	
	coh_n = abs.(results[:,[6]]);
	
	const_N = 4
	
	p_const_over_p_1 = abs_end_power ./ abs_end_power[const_N+1];
	
	
	"Plot values succesfully extracted"
end

# ╔═╡ 292ab05a-5fbe-11eb-2419-51a6c705dec7
let
	if(run_exp_1)
		plt = plot(all_N, p_N_over_p_1, title="Power gain for different slit counts")
		plot!(xlabel="N - Number of slits", ylabel =L"P^{(N)}/P^{(1)}")
		saveFigToDir("figs/exp1", "pn_over_p1.pdf", plt)
		plt
	end
end

# ╔═╡ 80e20acc-5fbf-11eb-0960-21686dd875d8
let
	if(run_exp_1)
		plt = plot(all_N, coh_n, title="Final global coherence for different slit counts")
		plot!(xlabel= "N - Number of slits", ylabel=L"\mu_g^{(N)}")
		saveFigToDir("figs/exp1", "coh_n_against_N.pdf", plt)
		plt
	end
end

# ╔═╡ bd94c1d0-5ff1-11eb-2863-ef94e26d9cc8
md"""
# Start and end coherence heatmaps for different slit counts
We now have some preliminary results that match our original paper, we can try and explore further.
For example, in the FORTRAN code we requested, it seems as though they choose to store the entire matrix A before and after running the experiment. We can do one better, choosing to plot it.
In the cell below, we print a set of heat maps that represent the coherence of the light after passing through different slit counts. We did not plot these in this notebook, instead choosing to save them. They can be found in 'figs/heatmaps/'
"""

# ╔═╡ 6cf6803c-607e-11eb-3bfc-69879c07310e
exp_2_args = GlobalArgs(51, 0.5, 0.2)

# ╔═╡ f3c717b4-6233-11eb-3f44-b95ab2db7951
let
	
	data_test = Data(exp_2_args.n, exp_2_args.zmax, 1, exp_2_args.dd)
	initExperimentValues(data_test)
	
	#calculate start values
 	startPower = calculatePower(data_test)
 	startCoh = calculateCoherence(data_test, startPower)
end

# ╔═╡ 0a5b3306-614f-11eb-08d6-07ac92e07c51
@bind run_heatmaps CheckBox()

# ╔═╡ 8a0550d6-5ff2-11eb-33b5-79c6ac0b4c63
if run_heatmaps
		
	file_path = "figs/heatmaps"
	file_path_no_cbar = "figs/heatmaps/nocbar"
	cbar_min = 0
	cbar_max = 0.045
	
	all_slits = append!([1], Array(5:5:25))
	#try between 1 and
	for (i, slit) in enumerate(all_slits)
		dat = Data(exp_2_args.n, exp_2_args.zmax, slit, exp_2_args.dd)
		
		sPow, sCoh, sInt, ePow, eCoh, eInt, startA, endA = runExperiment(dat)
		
		#If we are on our first slit count, also output
		if(i == 1)
			startA = abs.(startA) .^ 2 #convert from intensity to coherence like
			hm = heatmap(startA, aspect_ratio=:equal,  clim=(cbar_min, cbar_max))
			plt = plot(hm)
			plot!(xlabel= L"x_1", ylabel=L"x_2")

			saveFigToDir(file_path, "coherence_heatmap_source.pdf", plt)
			
				#save without cbar
			hm_no_cbar = heatmap(startA, aspect_ration=:equal, legend=:none, clim=(cbar_min, cbar_max))
			plot_no_cbar = plot(hm_no_cbar)
			plot!(xlabel= L"x_1", ylabel=L"x_2")
			saveFigToDir(file_path_no_cbar, "coherence_heatmap_$slit _slits_nocbar.pdf", plot_no_cbar)
			#savefig("$file_path /coherence_heatmap_source.pdf")
		end
		
		endA = abs.(endA) .^ 2#convert from intensity to coherence like

		#save file with cbar
		hm = heatmap(endA, aspect_ration=:equal, clim=(cbar_min, cbar_max))
		plt = plot(hm)
		plot!(xlabel= L"x_1", ylabel=L"x_2")
		saveFigToDir(file_path, "coherence_heatmap_$slit _slits.pdf", plt)
		
		#save without cbar
		hm_no_cbar = heatmap(endA, aspect_ration=:equal, legend=:none, clim=(cbar_min, cbar_max))
		plot_no_cbar = plot(hm_no_cbar)
		plot!(xlabel= L"x_1", ylabel=L"x_2")
		saveFigToDir(file_path_no_cbar, "coherence_heatmap_$slit _slits_nocbar.pdf", plot_no_cbar)

	end
		
end

# ╔═╡ ffa9ae64-6156-11eb-393b-2bf51a61773e
md"""
# Power gain for a single slit count with respect to different starting coherence
Below, we calculate and plot the way that the power increase ratio for a single slit count $P^{(4)}/P^{(1)}$ for different starting coherences of light.
"""

# ╔═╡ 9757caca-6083-11eb-1123-49788ed0e63d
exp_3_args = GlobalArgs(101, 0.5, 0)

# ╔═╡ 820aa19c-6f6a-11eb-0591-6d5de3b3587b
@bind exp_3_title CheckBox()

# ╔═╡ 267db182-6152-11eb-3284-f9f1ca65646a
@bind run_exp_3 CheckBox()

# ╔═╡ 405c5032-5ff8-11eb-1423-3392b17d1945
let
	if(run_exp_3)
		pyplot()


		#lets try changing dd while keeping everything else constant

		num_slits = 4
		dd_vals = [0.025,0.05,0.075, 0.1, 0.125, 0.15, 0.2, 0.3, 0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1]
		tot_len = size(dd_vals)
		#dd_vals = Array(range(0.001,stop=1,length=tot_len))
		coherence_results = zeros(tot_len) * 0.0
		coherence_results2 = zeros(tot_len) * 0.0

		power_results = zeros(tot_len) * 0.0
		for (i, dd) in enumerate(dd_vals)

			#this coherence, 4 slits
			dat_4 = Data(exp_3_args.n, exp_3_args.zmax, num_slits, dd); 

			 #this coherence, 1 slit
			dat_1 = Data(exp_3_args.n, exp_3_args.zmax, 1, dd);

			sPow4, sCoh4, sInt4, ePow4, eCoh4, eInt4 = runExperiment(dat_4)
			sPow1, sCoh1, sInt1, ePow1, eCoh1, eInt1 = runExperiment(dat_1)
			coherence_results[i] = sCoh4
			coherence_results2[i] = sCoh1

			power_results[i] = abs(ePow4)/abs(ePow1)

		end

		p = plot(coherence_results, power_results, seriestype = :scatter)
		plot!(coherence_results, power_results)
		plot!(xlabel=L"\mu_g", ylabel=L"P^{(4)}/P^{(1)}")
		
		
		
		if(exp_3_title)
			plot!(title="Preliminary result showing how power gain\n for 4 slits decreases with a more coherent light source")
			saveFigToDir("figs/changing_coh","p4_over_p1_against_coh.pdf", p)

		else
			saveFigToDir("figs/changing_coh/no_title","p4_over_p1_against_coh_nt.pdf", p)

		end

		p;
	end
end

# ╔═╡ d3ddbe50-6158-11eb-142e-fdb83d14c5c2
md"""
We see that this form of the classical zeno effect seems to be mroe prominant for less coherent light sources.
"""

# ╔═╡ d9fe2a80-6156-11eb-003a-4b83419b6b49
md"""
# Increasing range of calculation parameters
We have now created the main set of graphs we wish to plot, and so we can extend the range of paramters we calculate, to give more informative plots.
"""

# ╔═╡ f6574b6e-6082-11eb-31f1-77523e3b4f6d
exp_final_args  = GlobalArgs(101, 0.5, 0)

# ╔═╡ bb62463e-6f6a-11eb-1020-4f7086afaaf8
@bind final_title CheckBox()

# ╔═╡ c2e0bae2-614d-11eb-2e3d-afb90e332436
@bind run_final CheckBox()

# ╔═╡ a5f281b8-6005-11eb-0bda-09ff565b6cd9
begin
	
	if(run_final)
		dd_values = [0.025,0.075, 0.1, 0.125, 0.15, 0.2, 0.3, 0.4,0.5, 0.6, 0.7, 0.8, 0.9, 0.98]
		
		nr_values = Array(1:20)


		results_final = zeros(size(dd_values)[1], size(nr_values)[1], 7) * 0.0

		Threads.@threads for i in 1:size(dd_values)[1]
			for (j, nr) in enumerate(nr_values)

				data_dd_nr = Data(exp_final_args.n, exp_final_args.zmax, nr, dd_values[i])
				sPow, sCoh, sInt, ePow, eCoh, eInt = runExperiment(data_dd_nr)

				results_final[i, j, 1] = nr
				results_final[i, j, 2] = sPow;
				results_final[i, j, 3] = sCoh;
				results_final[i, j, 4] = sInt;

				results_final[i, j, 5] = ePow;
				results_final[i, j, 6] = eCoh;
				results_final[i, j, 7] = eInt;
			end

		end

	end
	
end

# ╔═╡ b4cdba96-6007-11eb-25fa-532ffb5c372a
let
	if(run_final)
		
		
		
		
		pyplot()
		p = plot(xlabel=L"N", ylabel=L"P^{(N)}/P^{(1)}")

		
		#we select higher valu
		mu_index_to_plot = append!([1], Array(6:2:size(results_final)[1]))

		for i in mu_index_to_plot
			#the n axis, number of slits for each part
			n_axis = results_final[i,:,1]

			#The final power for each slit number
			P_n = results_final[i,:,5]
			#Divide by P^1 to get P^N/P^1
			P_n_scaled = P_n ./ P_n[1]

			start_coh = round(results_final[i,:, 3][1],digits=2)
			plot!(n_axis, P_n_scaled, lab=L"\mu_g="*"$start_coh")# $start_coh")

		end
		plot!(size=(600,500))
		
		if(final_title)
			plot!(title="Plot showing power gain as a function of\n interference slit count, shown for multiple \nstarting coherence values")
		saveFigToDir("figs/final_exp", "pN_over_p1_against_N.pdf", p)
		else
			
		saveFigToDir("figs/final_exp/no_title", "pN_over_p1_against_N_nt.pdf", p)
		end
		
		p	
	end
end

# ╔═╡ dc1bcd2c-609d-11eb-1457-172e8da11277
let
	if(run_final)
		pyplot()

		p = plot(xlabel=L"\mu_g", ylabel=L"P^{(N)}/P^{(1)}")
		

		#iterate all excluding P^1, as we scale by P^1
		for i in 2:2:size(results_final)[2]


			#the coherence axis, starting coherence
			coh_axis = results_final[:,1,3]

			power_1 = results_final[:, 1, 5]

			power_i = results_final[:, i, 5]

			power_i_scaled = power_i ./ power_1
			n_ = Int(results_final[1, i, 1])

			plot!(coh_axis, power_i_scaled, lab="N=$n_")# $start_coh")
			#
		end
		plot!(size=(600,500))
		
		if(final_title)
			plot!(title="Plot showing power gain for different starting coherence\n values, for multiple numbers of slits")
			saveFigToDir("figs/final_exp", "pN_over_p1_against_start_coh.pdf", p)
		else
			saveFigToDir("figs/final_exp/no_title", "pN_over_p1_against_start_coh_nt.pdf", p)

		end
		p	
	end
end

# ╔═╡ fe0b1ca2-609e-11eb-04b0-4d624a098f95
let
	if(run_final)
		pyplot()

		p = plot(xlabel=L"N", ylabel=L"\mu_g^{(N)}")
		

		mu_index_to_plot = append!([1], Array(6:2:size(results_final)[1]))

		for i in mu_index_to_plot

			#x axis is N, number of slits in the experiment
			n_axis = results_final[i,:,1]

			#the coherence axis, end coherence
			coh_axis = results_final[i,:,6]



			start_coh = round(results_final[i,:, 3][1],digits=2)
			plot!(n_axis, coh_axis, lab=L"\mu_g"*" = $start_coh")
			
		end
		plot!(size=(600,500))
		if(final_title)
			plot!(title="Plot showing global degree of coherence as we \nincrease number of slits, shown for different \nstarting values of coherence")
			saveFigToDir("figs/final_exp", "global_coh_against_N.pdf", p)
		else
			saveFigToDir("figs/final_exp/no_title", "global_coh_against_N_nt.pdf", p)
		end
		p	
	end
end

# ╔═╡ 75213890-6152-11eb-2174-01d67465f3f8
md"""
# Jump in coherence and power at each slit
The final aspect of this effect we wish to show, is how the coherence of light changes whenever a slit is encountered. To do this, we must modify our iterate over slits function to allow it to calculate the coherence after each slit. We can then return all these intermitent coherence values to be used for plotting. We shall also return the power at each slit.
"""

# ╔═╡ 093a9858-60b0-11eb-0ca1-b92defa3bb4d
"""
Runs the calculation, iterating over all slits within this experiment to
calculate the final result at the detector
"""
function iterateOverSlitsAndGetCoh(dat::Data, print_all=false)
	
	coh_results = zeros(dat.nr)
	power_results = zeros(dat.nr)
	#iterate all slits
	for m in 1:dat.nr
		
		#Fill the first octant of e
		for i in 1:convert(Int64, (dat.n+1)/2)
			for j in 1:i
				dat.e[i,j] = 0
				for k in 1:dat.n
					dat.c[k,j] = 0
					for l in 1:dat.n
						dat.c[k,j ] += dat.a[k,l] * dat.b[l,j]
					end
					dat.e[i,j] += dat.d[i,k] * dat.c[k,j]
				end
			end
		end
		 
		#Fill the octant below
		for i in convert(Int64, (dat.n+3)/2):dat.n
			for j in 1:(dat.n-i+1)
				dat.e[i,j] = 0
				for k in 1:dat.n
					dat.c[k,j] = 0
					for l in 1:dat.n
						dat.c[k,j] += dat.a[k,l] * dat.b[l,j]
					end
					dat.e[i,j] += dat.d[i,k]* dat.c[k,j]
				end
			end
		end
		
		#Fill the final quadrants
		for i in convert(Int64, (dat.n+3)/2):dat.n
			for j in (dat.n-i+2):i
				dat.e[i,j] = conj(dat.e[dat.n+1-j, dat.n+1-i])
			end
		end
		
		#finalising data	
		for i in 1:dat.n
			for j in 1:dat.n
				if( j<= i)
					dat.a[i,j] = dat.e[i,j]
				else
					dat.a[i,j] = conj(dat.e[j,i])
				end
			end
		end
		power = calculatePower(dat)
		coh_results[m] = calculateCoherence(dat, power)
		power_results[m] = power
	end
	
	return coh_results, power_results
end

# ╔═╡ 91c139aa-60b0-11eb-389b-9b490c4402fe
exp_jump_args = GlobalArgs(41, 0.5, 0.2)

# ╔═╡ 0864ff5e-6f6d-11eb-0982-2d626aa2b43a
@bind run_jump_title CheckBox()

# ╔═╡ 3da4cef6-6150-11eb-0e44-d74baa81dcd5
@bind run_jump_test CheckBox()

# ╔═╡ e957f802-6152-11eb-30a3-331d53e32933
"Tick above to run calculation"

# ╔═╡ 7adc8a16-60b0-11eb-022e-efa13742de0a

if (run_jump_test)
	
	n_test = [1, 2, 5, 10, 20, 50]
	
	jump_coh_results = Any[]
	jump_power_results = Any[]

	jump_x = Any[]
	
	#iterate all slit counts
	for nr in n_test
		#push the x coordinates
		push!(jump_x, Array(range(0, stop=1, length=nr+1)))
		
		data = Data(exp_jump_args.n, exp_jump_args.zmax, nr, exp_jump_args.dd)
		initExperimentValues(data)
		#get start coherence
		start_power = calculatePower(data)
		start_coh = calculateCoherence(data, start_power)
		
		jump_coh = [start_coh]
		jump_pow = [start_power]
		jump_it_result, jump_it_power_res = iterateOverSlitsAndGetCoh(data)
		append!(jump_coh, jump_it_result)
		
		append!(jump_pow, jump_it_power_res)

		push!(jump_coh_results, jump_coh)
		push!(jump_power_results, jump_pow)

	end
	"Coherence jumps measured succesfully"
end

# ╔═╡ 7368859e-60b2-11eb-329c-a18737645b13
let
	if(run_jump_test)
		p=plot(xlabel=L"z_D", ylabel=L"\mu_g")
		
		for (i, j_coh) in enumerate(jump_coh_results)

			final_plot_x = Any[]
			final_plot_y = Any[]

			#first convert our raw data to jump format
			for (j, x) in enumerate(jump_x[i])
				if(j == 1)
					append!(final_plot_x, x)
					append!(final_plot_y, j_coh[j])
				else
					append!(final_plot_x, x)
					append!(final_plot_y, j_coh[j-1])
					append!(final_plot_x, x)
					append!(final_plot_y, j_coh[j])
				end
			end

			append!(final_plot_x, 1.1)
			append!(final_plot_y, final_plot_y[size(final_plot_y)[1]])

			final_plot_x .*= exp_jump_args.zmax

			label_val = n_test[i]
			plot!(final_plot_x, final_plot_y, label="N=$label_val")
		end
		
		if(run_jump_title)
			plot!(title="Plot showing variation in global coherence\n after each slit, for different total slit numbers")		
			saveFigToDir("figs/final_exp", "coh_jump_each_slit.pdf", p)
		else
			saveFigToDir("figs/final_exp/no_title", "coh_jump_each_slit_nt.pdf", p)
		end
		p
	end
end

# ╔═╡ 59f1b756-6153-11eb-0685-3163bdc3debc

let
	if(run_jump_test)
		
		p = plot(xlabel=L"z_D", ylabel=L"P")
		for (i, j_coh) in enumerate(jump_power_results)

			final_plot_x = Any[]
			final_plot_y = Any[]

			#first convert our raw data to jump format
			for (j, x) in enumerate(jump_x[i])
				if(j == 1)
					append!(final_plot_x, x)
					append!(final_plot_y, j_coh[j])
					append!(final_plot_x, x)
					append!(final_plot_y, j_coh[j+1])
				elseif(j < size(j_coh)[1])
					append!(final_plot_x, x)
					append!(final_plot_y, j_coh[j])
					append!(final_plot_x, x)
					append!(final_plot_y, j_coh[j+1])
				 else
				 	append!(final_plot_x, x)
				 	append!(final_plot_y, j_coh[j])	
				end
			end

			append!(final_plot_x, 1.1)
			append!(final_plot_y, final_plot_y[size(final_plot_y)[1]])

			final_plot_x .*= exp_jump_args.zmax

			label_val = n_test[i]
			plot!(final_plot_x, final_plot_y, label="N=$label_val")
		end
		
		if(run_jump_title)
			plot!(title="Plot showing variation in power \n after each slit, for different total slit numbers")
			saveFigToDir("figs/final_exp", "pow_jump_each_slit.pdf", p)
		else
			saveFigToDir("figs/final_exp/no_title", "pow_jump_each_slit_nt.pdf", p)
		end
		p
	end
end

# ╔═╡ Cell order:
# ╟─6ef5c61e-5fb5-11eb-3c1a-f16d8288c036
# ╠═5b4c4458-5f38-11eb-0998-db9dfa182660
# ╠═cbd62ed8-5f3c-11eb-15bb-131bdf424007
# ╠═31566570-5f3d-11eb-2bed-91e88dad9bfd
# ╠═2c65ab14-5f3d-11eb-3f8d-67f64d0bc471
# ╠═ac2673b0-5f3f-11eb-21c0-9d6bce434a99
# ╠═5bfe1386-5f3d-11eb-03e8-75f1ba6e9db9
# ╠═c9463304-5f3a-11eb-30e6-e1882e01e062
# ╠═28168b1a-5f3e-11eb-0f4a-493b9b0399e8
# ╠═5a86cac8-5f3f-11eb-32c0-41ea6c56e911
# ╠═33803354-607c-11eb-15ec-d71ffaa07000
# ╟─49a99602-607c-11eb-24da-6900e24db0ae
# ╟─5a461e62-607d-11eb-310a-c3fac838e1c9
# ╠═bf547a4e-5f3b-11eb-17ed-f1685ba4801d
# ╟─29cb5eb4-6159-11eb-1fba-3164718f1193
# ╟─313fa49a-5f40-11eb-3e48-a90c4320c965
# ╟─8d995928-5fb8-11eb-28ec-e90a1a5bf8de
# ╟─3b6f627e-6081-11eb-2b90-4feedc83d256
# ╟─87a90dbc-5fb8-11eb-1ff7-15608ff9fb10
# ╟─292ab05a-5fbe-11eb-2419-51a6c705dec7
# ╟─80e20acc-5fbf-11eb-0960-21686dd875d8
# ╟─bd94c1d0-5ff1-11eb-2863-ef94e26d9cc8
# ╠═6cf6803c-607e-11eb-3bfc-69879c07310e
# ╠═f3c717b4-6233-11eb-3f44-b95ab2db7951
# ╟─0a5b3306-614f-11eb-08d6-07ac92e07c51
# ╠═8a0550d6-5ff2-11eb-33b5-79c6ac0b4c63
# ╟─ffa9ae64-6156-11eb-393b-2bf51a61773e
# ╟─9757caca-6083-11eb-1123-49788ed0e63d
# ╠═820aa19c-6f6a-11eb-0591-6d5de3b3587b
# ╠═267db182-6152-11eb-3284-f9f1ca65646a
# ╠═405c5032-5ff8-11eb-1423-3392b17d1945
# ╟─d3ddbe50-6158-11eb-142e-fdb83d14c5c2
# ╟─d9fe2a80-6156-11eb-003a-4b83419b6b49
# ╠═f6574b6e-6082-11eb-31f1-77523e3b4f6d
# ╠═bb62463e-6f6a-11eb-1020-4f7086afaaf8
# ╠═c2e0bae2-614d-11eb-2e3d-afb90e332436
# ╟─a5f281b8-6005-11eb-0bda-09ff565b6cd9
# ╟─b4cdba96-6007-11eb-25fa-532ffb5c372a
# ╟─dc1bcd2c-609d-11eb-1457-172e8da11277
# ╠═fe0b1ca2-609e-11eb-04b0-4d624a098f95
# ╟─75213890-6152-11eb-2174-01d67465f3f8
# ╟─093a9858-60b0-11eb-0ca1-b92defa3bb4d
# ╠═91c139aa-60b0-11eb-389b-9b490c4402fe
# ╠═0864ff5e-6f6d-11eb-0982-2d626aa2b43a
# ╠═3da4cef6-6150-11eb-0e44-d74baa81dcd5
# ╠═e957f802-6152-11eb-30a3-331d53e32933
# ╠═7adc8a16-60b0-11eb-022e-efa13742de0a
# ╠═7368859e-60b2-11eb-329c-a18737645b13
# ╠═59f1b756-6153-11eb-0685-3163bdc3debc
