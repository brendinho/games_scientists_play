#=
	Brendon Phillips
	Post-doctoral Visitor
	ABM-Lab
	Department of Mathsmatics and Statistics
	York University

	Julia simulation code to accompany the paper "Games scientists play ... "
	The model is fully described in the paper

	threaded code; start julia with -nthreads <number>
=#

using StatsBase
using DataFrames
using CSV
using Random
using Glob

M_A = 1 # funding available to acceptors
M_R = 0.34 # funding available to rejecters

#=
	One research question was to find out whether placing a cluster of acceptors (instead
	of distributing them randomly throughout the network) would influence eventual idea
	uptaske (depending on the interaction radius $r$). To run the simulation without this
	initial cluster, set M_vector=[0]. Else, initial $M\times M$ clusters will be placed
	on each network, with the other acceptors (total, proportion $p$ of the network) distributed
	randomly throughout the rest of the network.

	If $M^2>p.N$ (i.e., the cluster is larger than the desired number of initial acceptors),
	then the initial cluster will still be placed, and no acceptors will be placed on the network
	outside that cluster. However, the initial proportion in the file name and $p$ column of the
	output file will be incorrect. For analysis, a more appropriate guess of the initial proportion
	of acceptors would be $max(M^2, p.N)$.

	M is the linear dimension of the initial acceptor cluster, centred at (0,0), without loss of
	generality (due to the periodic boundary).
=#
M_vector = [0:10;]

D = 50 # linear dimension

time_limit = 10000 # number of time steps for the simulation

ensemble_size = 25 # number of trials per parameter tuple

neighbourhood_radius_vector = [1:20;] # 1:24 # Moore neighbourhood (interaction radius) size $r$
initial_acceptor_proportion_vector = [0.3,0.4] # initial proportion of accepting scientists

s_vector = [convert(Float64, x) for x in [0:10;]] # signal strength $s$
delta_vector = [convert(Float64, x) for x in [0:0.5:5;]] # acceptance payoff $delta$

#=
	get a list of the trials previously run and exclude them from the trials to be run
=#
previous_Parameters_hashes = [] # store the hashes of the pre-existing file names
list_of_previous_files = glob("Data_Files_Short/trial*csv") # all the trials that were previously run

#=
	for all the trials that were previously run, get a single row giving all the parameter values used
	and store its hash in a vector of hashes (checking hashes os much faster than comparing the tables
	themselves, or even rows)
=#
for data_file_name in list_of_previous_files
	data_file = CSV.read(data_file_name, DataFrame)
	append!(
		previous_Parameters_hashes,
		# hashing the entire table, not a row
		hash( unique(select(data_file, Not([:trial, :time, :prop_acceptors, :M_R, :M_A, :kappa]))) )
	)
end
previous_Parameters_hashes = unique(previous_Parameters_hashes)

#=
	for all the parameters we want to run, check if a results file with those parameter values in the name
	already exists, that is, the hash of the parameter tuple appears in the list of parameters already run.
	if so, skip; else, add to list of parameters to run
 =#

Parameters = DataFrame()
# for each possible combination of parameters
for r in neighbourhood_radius_vector, p in initial_acceptor_proportion_vector, s in s_vector, delta in delta_vector, M in M_vector
	# create a temp data frame with those parameter values for a check
	temp = DataFrame(newwork_linear_dimension=D, cluster_linear_dimension=M, initial_acceptor_proportion=p, s=s, interaction_radius=r, delta=delta)
	# if the hash is found, it's already been run. else, go ahead
	if !(hash(temp) in previous_Parameters_hashes)
		push!(Parameters, eachrow(temp)[1])
	end
end

#=
	get the radius $r$ Moore neighbood (i.e., possible queen's walks of max length $r$) centred at some node
	"academic" (from 1 to N). called "interaction neighbourhood" in the manuscript. node number is turned into
	square lattice coordinates, and neighbourhoods are found with periodic boundaries. lattice coordinates are
	then all returned as node numbers.

	i.e. with D=50,
	get_Moore_neighbourhood(5,2) -> 3, 4, 6, 7, 53, 54, 55, 56, 5, 10, 10, 10, 10, 107, 2403, 2404,
		2405, 2406, 2407, 2453, 2454, 2455, 2456, 2457
=#
function get_Moore_neighbourhood(academic::Int64, neighbourhood_radius::Int64)
	# get hte square lattice coordinates (starting at 1 in the upper left corner: 1-2-3-...)
	(y_coord, x_coord) = divrem(academic-1, D) .+ (1,1)
	# vector to store the neighbours of the node
	neighbours = Vector{Int64}()
	# temp function to turn latice coordinates into the corresponsing node number
	into_node(coords::Tuple{Int64, Int64}) = D*(coords[2]-1) + coords[1]
	# use successive hosizontal lines centred on the column containing the node
	for x in (x_coord - neighbourhood_radius):(x_coord + neighbourhood_radius)
		for y in (y_coord-neighbourhood_radius):(y_coord + neighbourhood_radius)
			push!(neighbours, into_node((
				if(x < 1) x + D elseif(x > D) x - D else x end, # +/-D for periodic boundaries
				if(y < 1) y + D elseif(y > D) y - D else y end # ... ^ ...
			)))
		end
	end
	 # sorted list of neighbour hodes (less the node itself)
	return sort(filter(x -> x != academic, unique(neighbours)))
end

Threads.@threads for combo in eachrow(Parameters)

	# random number generator
	rng = MersenneTwister(floor(Int, time()))

	# Parameter values for this trial
	M = combo.cluster_linear_dimension # dimension of the initial acceptor cluster
	p = combo.initial_acceptor_proportion # initial proportion of acceptors
	s = combo.s # social norm strrength
	r = combo.interaction_radius # dimension of the network total
	delta = combo.delta # utility of acceptance

	# #=
	# 	recording the comments at each time step is time-consuming, so we default just
	# 	record the final time step for each trial. to get the full time series, uncomment
	# 	every instance of Results_Here (the reporting in the loop and the file output at the end)
	# =#
	# Results_Here = DataFrame(fill((
	# 	M_A = -1.0,
	# 	M_R = -1.0,
	# 	kappa = -1.0,
	# 	linear_dimension = -1,
	# 	cluster_dimension = -1,
	# 	trial = -1,
	# 	initial_acceptor_proportion = -1.0,
	# 	s = -1.0,
	# 	interaction_radius = -1,
	# 	delta = -1.0,
	# 	time = -1,
	# 	prop_acceptors = -1.0
	# ), time_limit*ensemble_size))
	# #=
	# 	set the initial size rather than growing the table successively.
	# 	pre-allocation is always better
	# =#

	Results_Final = DataFrame(
		M_A = Float64[],
		M_R = Float64[],
		kappa = Float64[],
		network_linear_dimension = Int[],
		cluster_linear_dimension = Int[],
		trial = Int[],
		initial_acceptor_proportion = Float64[],
		s = Float64[],
		interaction_radius = Int[],
		delta = Float64[],
		time = Int[],
		prop_acceptors = Float64[]
	)

	# for each trial in the ensemble
	for trial in 1:ensemble_size

		# list of neighbours
		Neighbours = Vector{Vector{Int64}}(undef, D^2)
		for academic in 1:D^2
			Neighbours[academic] = get_Moore_neighbourhood(academic, r)
		end

		# gset the initial acceptance or rejection of the idea for all scientists
		Opinions = repeat([-1], D^2)
		if M == 0
			# if no cluster, choose at random to fill the quota
			Opinions = [sample([1,-1], Weights([p, 1-p])) for _ in 1:D^2]
		else
			# place the acceptor cluster
			number_acceptors = Int(floor(p*D^2))
			# get the initial_cluster of acceptors
			cluster = get_Moore_neighbourhood(1, M)
			Opinions[cluster] .= 1
			# set the rest of the randomly acceptors to acceptance
			Opinions[ shuffle(setdiff(1:D^2, cluster))[1:(number_acceptors-length(cluster))] ] .= 1
		end

		# payoff function for accepting  - fully described in the manuscript
		function payoff_for_Accepting(academic::Int64)
			#=
				depends on the social norm created by the acceptors in your interaction
				neighbourhood + the utility of accepting + the expected amount of funding
			=#
			num_Acceptors = count(x -> x == 1, Opinions)
			num_Accepting_neighbours = count(x -> x == 1, getindex(Opinions, Neighbours[academic]))

			first_fraction = M_A/(num_Acceptors + (Opinions[academic] == -1))
			second_fraction = num_Accepting_neighbours/((2*r+1)^2-1)

			return first_fraction + s*second_fraction + delta
		end

		# payoff function for rejecting the idea - see manuscript
		function payoff_for_Rejecting(academic::Int64)
			#=
				depends on the social norm created by the rejectors in your interaction
				neighbourhood + the expected amount of funding
			=#
			num_Rejectors = count(x -> x == -1, Opinions)
			num_Rejecting_neighbours = count(x -> x == -1, getindex(Opinions, Neighbours[academic]))

			first_fraction = M_R/(num_Rejectors + (Opinions[academic] == 1))
			second_fraction = num_Rejecting_neighbours/((2*r+1)^2-1)

			return first_fraction + s*second_fraction
		end

		for time in 1:time_limit

			# choose a random node
			random_academic = rand(rng, 1:D^2)

			# change to the opinion with higher payoff
			if payoff_for_Accepting(random_academic) == payoff_for_Rejecting(random_academic)
				continue
			elseif payoff_for_Accepting(random_academic) > payoff_for_Rejecting(random_academic)
				Opinions[random_academic] = 1
			else
				Opinions[random_academic] = -1
			end

			# # uncomment to record each time steps
			# Results_Here[(trial-1)*time_limit + time, :] = (
			# 	M_A, # acceptor resources
			# 	M_R, # rejector_resources
			# 	M_R/M_A, # kappa
			# 	D, # network linear dimension
			# 	M, # cluster linear dimension
			# 	trial, # trial number in the ensemble
			# 	p, # initial proportion of acceptors
			# 	s, # signal strength
			# 	r, # neighbourhood radius
			# 	delta, # delta
			# 	time, # time step
			# 	count(x -> x == 1, Opinions)/D^2 # final proportion of acceptors
			# )

		end

		push!(Results_Final,(
			M_A, # acceptor resources
			M_R, # rejector_resources
			M_R/M_A, # kappa
			D, # network linear dimension
			M, # cluster linear dimension
			trial, # trial number in the ensemble
			p, # initial proportion of acceptors
			s, # signal strength
			r, # neighbourhood radius
			delta, # delta
			time_limit, # time step
			count(x -> x == 1, Opinions)/D^2 # final proportion of acceptors
		))

	end

	CSV.write("Data_Files_Short/trial_p_$(p)_s_$(s)_r_$(r)_delta_$(delta)_D_$(D)_M_$(M).csv", Results_Final)
	# # uncomment to record the per-time-step data
	# CSV.write("Data_Files_Short/times_p_$(p)_s_$(s)_r_$(r)_delta_$(delta)_D_$(D)_M_$(M).csv", Results_Here)

end
