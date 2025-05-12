module Visualization
using Plots
using FileIO
using Dates
using Measures  # Add this import for mm unit
using ..BoundaryConditions
using ..Constants

export plot_configuration, plot_energy_history, save_visualizations

"""
    plot_configuration(positions, moments; title="", savepath=nothing)
Create 3D plot of particle configuration with optional saving.
"""
function plot_configuration(positions, moments; title="Ferrofluid Configuration", savepath=nothing)
    plt = plot3d(
        legend=false,
        xlabel="x (m)", ylabel="y (m)", zlabel="z (m)",
        title=title,
        camera=(30, 30),
        size=(800, 600),
        marker=:circle,
        markersize=3
    )
    
    # Plot particles with color-coded magnetic moments
    for (pos, mom) in zip(positions, moments)
        # Normalize moment for color mapping
        color = RGB((mom[1]+1)/2, (mom[2]+1)/2, (mom[3]+1)/2)
        scatter!(
            plt,
            [pos[1]], [pos[2]], [pos[3]],
            markercolor=color,
            markersize=4,
            markerstrokewidth=0
        )
    end
    
    # Add colorbar for magnetic moments
    plot!(plt, colorbar=true, colorbartitle="Magnetic Moment")
    
    if !isnothing(savepath)
        mkpath(dirname(savepath))
        savefig(plt, savepath)
        @info "Configuration plot saved to $savepath"
    end
    
    return plt
end

"""
    plot_energy_history(energy_history; title="", savepath=nothing)
Plot energy convergence with optional saving.
"""
function plot_energy_history(energy_history; title="Energy Convergence", savepath=nothing)
    plt = plot(
        energy_history,
        xlabel="MC Step",
        ylabel="Energy (J)",
        title=title,
        label="System Energy",
        linewidth=2,
        size=(800, 400),
        margin=5mm  # Now properly using the imported mm unit
    )
    
    if !isnothing(savepath)
        mkpath(dirname(savepath))
        savefig(plt, savepath)
        @info "Energy plot saved to $savepath"
    end
    
    return plt
end

"""
    save_visualizations(results, sys_params; prefix="ferrofluid_", directory="results")
Generate and save all visualizations for a simulation.
"""
function save_visualizations(results, sys_params; prefix="ferrofluid_", directory="results")
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    basepath = joinpath(directory, prefix * timestamp)
    
    # Create output directory
    mkpath(directory)
    
    # 1. Configuration plot
    config_path = basepath * "_config.png"
    plot_configuration(
        results.final_config.positions,
        results.final_config.moments,
        title="Final Configuration (N=$(sys_params.N))",
        savepath=config_path
    )
    
    # 2. Energy convergence plot
    energy_path = basepath * "_energy.png"
    plot_energy_history(
        results.energy_history,
        title="Energy Convergence (Acceptance: $(round(results.acceptance_rate*100, digits=1))%)",
        savepath=energy_path
    )
    
    # 3. Combined summary plot
    summary_path = basepath * "_summary.png"
    plt_summary = plot(
        plot_configuration(results.final_config.positions, results.final_config.moments),
        plot_energy_history(results.energy_history),
        layout=(1,2),
        size=(1200, 500),
        plot_title="Ferrofluid Simulation Summary - $(timestamp)"
    )
    savefig(plt_summary, summary_path)
    @info "Summary plot saved to $summary_path"
    
    return (config_path=config_path, energy_path=energy_path, summary_path=summary_path)
end

end  # module Visualization
