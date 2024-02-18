using Catalyst, DifferentialEquations, Plots, Random, Distributions, DataFrames, CSV, Printf, ModelingToolkit
# Function to initialize the reaction network
function initialize_rn()
    # Define the parameters, variables, and species of the reaction network
    @parameters kOn kOff k kT deg_R deg_G;
    @variables t;
    @species A(t) DNA(t) A_DNA(t) RNA(t) GFP(t);

    # Define the reaction network topology
    rxs = [
        (@reaction kOn, A + DNA --> A_DNA), # binding of Transcription factor complex to DNA
        (@reaction kOff, A_DNA --> A + DNA), # unbinding of Transcription factor complex to DNA
        (@reaction k, A_DNA --> A_DNA + RNA), # Transcription of DNA to mRNA
        (@reaction kT, RNA --> RNA + GFP), # Translation of RNA to reporter protein (GFP)
        (@reaction deg_R, RNA --> 0), # mRNA Degradation
        (@reaction deg_G, GFP --> 0) # Reporter Degradation
    ]

    # Define the reaction system
    @named rn = ReactionSystem(rxs, t, [A, DNA, A_DNA, RNA, GFP], [kOn, kOff, k, kT, deg_R, deg_G])

    return rn
end

# Call reaction network initilization function
rn = initialize_rn();

#use this function to generate the normally distributed parameters
function generate_parameter_samples(means::Tuple, std_devs::Tuple, num_samples::Int)
    param_values = [rand(Normal(mean, std), num_samples) for (mean, std) in zip(means, std_devs)]
    return param_values
end

#use this function to display the distributions as histogram
function display_histograms(param_values::Vector{Vector{Float64}})
    histograms = []
    param_labels = ["kOn", "kOff", "k", "kT", "deg_R", "deg_G"]
    for (i, values) in enumerate(param_values)
        hist = histogram(values, label=param_labels[i], xlabel="Parameter Value", ylabel="Frequency", bins=20)
        push!(histograms, hist)
    end

    # Display all histograms
    plot(histograms..., layout=(2, 3), legend=:topright, size=(1200, 800), spacing=10)  # Adjust left margin here
end

# Values used for generate_parameter_samples function
num_samples = 1000; #of samples to draw from for each parameter
means = (10.0, 10.0, 5.0, 2.0, 0.08, 0.05); # mean value of parameter distribution (kOn, kOff, k, kT, deg_R, deg_G)
std_devs = (1.0, 1.0, 0.5, 0.035, 0.001, 0.01); # Standard deviation of parameter distribution

#generate the distributions of parameters (priors) to be selected from randomly at each iteration of the simulation loop
param_values = generate_parameter_samples(means, std_devs, num_samples);
# Display distributrions
display_histograms(param_values)

#Function to call to run the simulations
function run_simulations(num_simulations::Int, param_values::Vector{Vector{Float64}}, rn::ReactionSystem, u0::Vector{Pair{Num, Int64}}, tspan::Tuple{Float64, Float64}, saveat::StepRangeLen{Float64, Float64})
    all_species_values = []
    all_time_points = []
    param_sets = []

    Random.seed!(123) # Set a seed for reproducibility

    for i in 1:num_simulations
        # Randomly choose values from the pre-generated arrays for each parameter
        random_param_values = [values[rand(1:length(values))] for values in param_values]
        push!(param_sets, random_param_values) # Store the parameter set for this simulation

        # Create a single DiscreteProblem and JumpProblem
        p = Dict(Symbol("kOn") => random_param_values[1],
                 Symbol("kOff") => random_param_values[2],
                 Symbol("k") => random_param_values[3],
                 Symbol("kT") => random_param_values[4],
                 Symbol("deg_R") => random_param_values[5],
                 Symbol("deg_G") => random_param_values[6])

        dprob = DiscreteProblem(rn, u0, tspan, p)
        jprob = JumpProblem(rn, dprob, Direct())

        # Solve the problem using SSAStepper
        sol = solve(jprob, SSAStepper(), saveat=saveat)
        push!(all_species_values, sol(tspan[1]:saveat[2]:tspan[2])) # Store all species values at each time point
        push!(all_time_points, tspan[1]:saveat[2]:tspan[2]) # Store all time points
    end

    return all_species_values, all_time_points, param_sets
end

#set the number of simulations to be run
num_simulations = 5
# Set span of time (in hours - parameter rates calculated in hr^-1)
tspan = (0,1000);
saveat = 0:1:1000
#  Set initial conditions of species in the reaction network
u0 = [A => 1, DNA => 1, A_DNA => 0, RNA => 0, GFP => 0];
#call function to run simulations
solutions, all_time_points, param_sets = run_simulations(num_simulations, param_values, rn, u0);

#Store each species values for all time points as individual arrays
Sim_A = []
Sim_DNA = []
Sim_ADNA = []
Sim_RNA = []
Sim_GFP = []
for sim in 1:length(solutions)
    species_values_GFP = []
    species_values_A = []
    species_values_DNA = []
    species_values_ADNA = []
    species_values_RNA = []
    for t in 1:length(solutions[sim])
        push!(species_values_GFP, solutions[sim][t][5])
        push!(species_values_A, solutions[sim][t][1])
        push!(species_values_DNA, solutions[sim][t][2])
        push!(species_values_ADNA, solutions[sim][t][3])
        push!(species_values_RNA, solutions[sim][t][4])
    end
    push!(Sim_GFP, species_values_GFP)
    push!(Sim_A, species_values_A)
    push!(Sim_DNA, species_values_DNA)
    push!(Sim_ADNA, species_values_ADNA)
    push!(Sim_RNA, species_values_RNA)
end

# Initialize an empty array to store the lengths
lengths = []
# Iterate over each array in Sim_GFP
for sim_values in Sim_GFP
    # Get the length of the current array and store it
    push!(lengths, length(sim_values))
end
# Print the lengths
println(lengths)
# Define the time point from which you want to select values
SS_time = 250
# Initialize an empty array to store the selected values
Sim_GFP_SS = []
# Iterate over each simulation
for sim_values in Sim_GFP
    # Slice the values to select only those after the start_time
    selected_values = sim_values[SS_time:end]
    # Append the selected values to the new array
    push!(Sim_GFP_SS, selected_values)
end
plot(Sim_GFP_SS)