using Catalyst, DifferentialEquations, Plots, Random, Distributions, DataFrames, CSV

@parameters kOn kOff k kT deg_R deg_G;
@variables t;
@species A(t) DNA(t) A_DNA(t) RNA(t) GFP(t);
# Define the reaction network
rxs = [
    (@reaction kOn, A + DNA --> A_DNA), # binding of Transcription factor complex to DNA
    (@reaction kOff, A_DNA --> A + DNA), # unbinding of Transcription factor complex to DNA
    (@reaction k, A_DNA --> A_DNA + RNA), # Transcription of DNA to mRNA
    (@reaction kT, RNA --> RNA + GFP), # Translation of RNA to reporter protein (GFP)
    (@reaction deg_R, RNA --> 0), # mNRA Degradation
    (@reaction deg_G, GFP --> 0) # Reporter Degradation
];

tspan = (0.0, 100); # reaction time span
u0 = [A => 1, DNA => 1, A_DNA => 0, RNA => 0, GFP => 0];  # starting conditions

@named rn = ReactionSystem(rxs, t, [A, DNA, A_DNA, RNA, GFP], [kOn, kOff, k, kT, deg_R, deg_G]);

# Generate arrays of ~1000 values for each parameter using normal distribution
num_samples = 1000;

means = (10.0, 10.0, 5.0, 2.0, 0.08, 0.03); # mean value of parameter distribution (kOn, kOff, k, kT, deg_R, deg_G)
std_devs = (1.0, 1.0, 0.5, 0.35, 0.001, 0.001); # Standard deviation of parameter distribution

# Create arrays of ~1000 values for each parameter
param_values = [rand(Normal(mean, std), num_samples) for (mean, std) in zip(means, std_devs)];

# Plot histograms for each parameter
histograms = [];
param_labels = ["kOn", "kOff", "k", "kT", "deg_R", "deg_G"];
for (i, values) in enumerate(param_values);
    hist = histogram(values, label=param_labels[i], xlabel="Parameter Value", ylabel="Frequency", bins=20);
    push!(histograms, hist);
end

# Display all histograms
plot(histograms..., layout=(2, 3), legend=:topright)

# Arrays to store simulated GFP values and parameters for each simulation
num_simulations = 5;
solutions = [];
param_sets = [];

Random.seed!(123); # Set a seed for reproducibility
for i in 1:num_simulations;
    # Randomly choose values from the pre-generated arrays for each parameter
    random_params = [values[rand(1:num_samples)] for values in param_values];
    push!(param_sets, random_params); # Store the parameter set for this simulation

    # Create a single DiscreteProblem and JumpProblem
    p = Dict(Symbol("kOn") => random_params[1],
             Symbol("kOff") => random_params[2],
             Symbol("k") => random_params[3],
             Symbol("kT") => random_params[4],
             Symbol("deg_R") => random_params[5],
             Symbol("deg_G") => random_params[6])

    dprob = DiscreteProblem(rn, u0, tspan, p)
    jprob = JumpProblem(rn, dprob, Direct())

    # Solve the problem using SSAStepper
    sol = solve(jprob, SSAStepper())
    push!(solutions, sol)
end

# Preallocate array to store GFP values for each simulation
num_steps = length(tspan[1]:0.1:tspan[2]);
gfp_values = zeros(num_steps, num_simulations);

# Store GFP values from simulations
for i in 1:num_simulations
    sol = solutions[i]
    for j in 1:num_steps
        t = tspan[1] + (j - 1) * 0.1
        gfp_values[j, i] = sol(t)[5]
    end
end

# Create a dictionary with column names and GFP values for each simulation
column_names = ["Simulation_$i" for i in 1:num_simulations];
gfp_dict = Dict(Symbol(column_names[i]) => gfp_values[:, i] for i in 1:num_simulations);


# Initialize an empty DataFrame with the parameter names as the first column
param_df = DataFrame(Parameter = param_labels);

# Add parameter values to the DataFrame for each simulation run
for i in 1:num_simulations
    sim_col = Symbol("Simulation_$i");
    param_df[!, sim_col] = [param_sets[i][j] for j in 1:length(param_labels)];
end
# Save the parameters used for each simulation to a CSV file
CSV.write("C://Users//Strey Lab//Desktop//Github Reps//Outputs//12_14_23//parameters.csv", param_df);

# Create a DataFrame with all simulated GFP values
df_gfp_values = DataFrame(gfp_dict);
# Save the simulated GFP time traces to a CSV file
CSV.write("C://Users//Strey Lab//Desktop//Github Reps//Outputs//12_14_23//simulated_traces.csv", df_gfp_values, bufsize=100000000);
rename!(df_gfp_values, Symbol.(column_names));

# Create a plot object
plot_obj = plot(size=(1000,800),
    tickfont=font(12),  # Adjust the font size of ticks
    xlabelfont=font(14),  # Adjust the font size of x-axis label
    ylabelfont=font(14),  # Adjust the font size of y-axis label
    titlefont=font(16)    # Adjust the font size of title
);
# Plot the results for each simulation with parameter values in the legend
for i in 1:num_simulations
    # Extract GFP values from the solution and plot
    sol = solutions[i]
    plot!(sol, vars=:GFP, label="Simulation $i")
end
# Customize plot labels and title, add legend manually
xlabel!("Time");
ylabel!("GFP");
title!("GFP over Time for Multiple Simulations");
plot!(legend=:topleft);
# Display the final plot
display(plot_obj)
