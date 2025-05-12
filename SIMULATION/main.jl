using Dates
using Printf
using DelimitedFiles

# Load all modules
include("constants.jl")
include("boundary_conditions.jl")
include("forces.jl")
include("monte_carlo.jl")
include("visualization.jl")

using .Constants, .BoundaryConditions, .Forces, .MonteCarlo, .Visualization

function main(;
    N=50,                  # Number of particles
    box_size=1e-6,         # Simulation box size (m)
    pH=7.0,                # System pH
    ionic_strength=0.01,   # Ionic strength (mol/L)
    T=298.15,              # Temperature (K)
    H_ext=[0.0, 0.0, 1e4], # External field (A/m)
    n_steps=20_000,        # MC steps
    output_dir="results"   # Output directory
)
    try
        # Create output directory with timestamp
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
        output_path = joinpath(output_dir, "simulation_$timestamp")
        mkpath(output_path)

        @info "Starting simulation with $N particles at pH $pH"
        @info "Output will be saved to: $output_path"

        # 1. Set up system parameters
        sys_params = BoundaryConditions.SystemParameters(
            N=N,
            box_size=box_size,
            pH=pH,
            ionic_strength=ionic_strength,
            T=T,
            H_ext=H_ext,
            eta=1e-3,
            h_surfactant=2e-9
        )

        # 2. Initialize Monte Carlo parameters
        mc_params = MonteCarlo.MCSimulation(
            n_steps=n_steps,
            step_size=0.2*Constants.a
        )

        # 3. Run simulation
        @info "Running Monte Carlo simulation..."
        start_time = time()
        results = MonteCarlo.run_simulation(sys_params, mc_params)
        elapsed_time = time() - start_time

        @info "Simulation completed in $(@sprintf("%.2f", elapsed_time)) seconds"
        
        # 4. Save parameters to file
        params_path = joinpath(output_path, "simulation_parameters.txt")
        open(params_path, "w") do f
            println(f, "=== Ferrofluid Simulation Parameters ===")
            println(f, "Timestamp: ", timestamp)
            println(f, "\n[System Parameters]")
            println(f, "Particles (N): ", N)
            println(f, "Box size (m): ", box_size)
            println(f, "pH: ", pH)
            println(f, "Ionic strength (M): ", ionic_strength)
            println(f, "Temperature (K): ", T)
            println(f, "External field (A/m): ", H_ext)
            println(f, "Viscosity (Pa·s): ", 1e-3)
            println(f, "Surfactant thickness (m): ", 2e-9)
            println(f, "Particle radius (m): ", Constants.a)
            println(f, "Saturation magnetization (A/m): ", Constants.M_s)
            
            println(f, "\n[Monte Carlo Parameters]")
            println(f, "Total steps: ", n_steps)
            println(f, "Equilibration steps: ", mc_params.equilibration_steps)
            println(f, "Step size (m): ", mc_params.step_size)
            println(f, "Target acceptance: ", mc_params.target_acceptance)
            println(f, "Sampling interval: ", mc_params.sampling_interval)
            
            println(f, "\n[Results Summary]")
            println(f, "Acceptance rate: ", round(results.acceptance_rate*100, digits=2), "%")
            println(f, "Final energy (J): ", results.energy_history[end])
            println(f, "Energy min/max (J): ", extrema(results.energy_history))
        end
        @info "Parameters saved to: $params_path"

        # 5. Save numerical results
        data_path = joinpath(output_path, "energy_history.csv")
        writedlm(data_path, results.energy_history, ',')
        @info "Energy history saved to: $data_path"

        # 6. Generate visualizations
        viz_paths = Visualization.save_visualizations(
            results,
            sys_params,
            directory=output_path
        )

        @info """
        Visualization outputs:
        - Configuration plot: $(viz_paths.config_path)
        - Energy plot: $(viz_paths.energy_path)
        - Summary plot: $(viz_paths.summary_path)
        """

        return (parameters=params_path, data=data_path, visualizations=viz_paths)

    catch e
        @error "Simulation failed!" exception=(e, catch_backtrace())
        return nothing
    end
end

# Run with default parameters if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
