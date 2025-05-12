module MonteCarlo
using ..Constants
using ..BoundaryConditions
using ..Forces
using Distributions
using LinearAlgebra
using Random
using Statistics

export MCSimulation, run_simulation, generate_initial_configuration

"""
    MCSimulation
Parameters for Monte Carlo simulation.
"""
struct MCSimulation
    n_steps::Int
    equilibration_steps::Int
    step_size::Float64
    target_acceptance::Float64
    sampling_interval::Int

    # Add keyword constructor
    function MCSimulation(;
        n_steps::Int,
        equilibration_steps::Int = div(n_steps, 5),
        step_size::Float64,
        target_acceptance::Float64 = 0.5,
        sampling_interval::Int = 10
    )
        new(n_steps, equilibration_steps, step_size, target_acceptance, sampling_interval)
    end
end

"""
    generate_initial_configuration(sys_params::SystemParameters)
Create initial particle configuration.
"""
function generate_initial_configuration(sys_params::SystemParameters)
    positions = [sys_params.box_size * rand(3) for _ in 1:sys_params.N]
    moments = [sys_params.m_s * normalize(randn(3)) for _ in 1:sys_params.N]
    
    # Ensure no overlaps
    for i in 1:sys_params.N, j in i+1:sys_params.N
        while norm(positions[i] - positions[j]) < 2.2*Constants.a
            positions[j] = sys_params.box_size * rand(3)
        end
    end
    
    return (positions=positions, moments=moments)
end

"""
    calculate_system_energy(positions, moments, sys_params, derived_params)
Compute total energy of the system.
"""
function calculate_system_energy(positions, moments, sys_params, derived_params)
    energy = 0.0
    for i in 1:sys_params.N, j in i+1:sys_params.N
        r = BoundaryConditions.minimum_image_distance(positions[i], positions[j], sys_params.box_size)
        d = norm(r) - 2*a
        r_hat = r / norm(r)
        
        # Steric
        if d < sys_params.h_surfactant
            energy += derived_params.U_rep_steric * exp(-d/sys_params.h_surfactant)
        end
        
        # DLVO
        if d > 0
            kappa = 1/derived_params.lambda_D
            energy += 64*pi*epsilon_0*epsilon_r*a*(derived_params.k_B_T/e)^2 *
                     tanh(e*derived_params.zeta/(4*derived_params.k_B_T))^2 *
                     exp(-kappa*d)
            energy += -derived_params.A_H*a/(12*d)
        end
        
        # Magnetic
        energy += (mu_0/(4*pi)) * (
            dot(moments[i], moments[j])/norm(r)^3 -
            3*dot(moments[i], r_hat)*dot(moments[j], r_hat)/norm(r)^3
        )
    end
    
    # External field
    for mom in moments
        energy += -dot(mom, sys_params.H_ext)
    end
    
    return energy
end

"""
    run_simulation(sys_params::SystemParameters, mc_params::MCSimulation)
Run Metropolis Monte Carlo simulation.
"""
function run_simulation(sys_params::SystemParameters, mc_params::MCSimulation)
    derived_params = BoundaryConditions.initialize_system(sys_params)
    config = generate_initial_configuration(sys_params)
    positions, moments = deepcopy(config.positions), deepcopy(config.moments)
    
    # Trackers
    current_energy = calculate_system_energy(positions, moments, sys_params, derived_params)
    energy_history = Float64[]
    accepted = 0
    step_size = mc_params.step_size
    
    for step in 1:mc_params.n_steps
        # Random particle move
        p_idx = rand(1:sys_params.N)
        new_positions = deepcopy(positions)
        new_moments = deepcopy(moments)
        
        # Position or orientation move
        if rand() < 0.5
            # Position move (Gaussian)
            new_positions[p_idx] += step_size * randn(3)
            BoundaryConditions.apply_periodic_boundary!([new_positions[p_idx]], sys_params.box_size)
        else
            # Orientation move (Von Mises-Fisher)
            axis = normalize(randn(3))
            angle = step_size * randn()
            new_moments[p_idx] = cos(angle)*moments[p_idx] + 
                                sin(angle)*cross(axis, moments[p_idx]) +
                                (1-cos(angle))*dot(axis, moments[p_idx])*axis
            new_moments[p_idx] *= sys_params.m_s / norm(new_moments[p_idx])
        end
        
        # Energy change
        new_energy = calculate_system_energy(new_positions, new_moments, sys_params, derived_params)
        ΔE = new_energy - current_energy
        
        # Metropolis criterion
        if ΔE < 0 || rand() < exp(-ΔE/derived_params.k_B_T)
            positions, moments = new_positions, new_moments
            current_energy = new_energy
            accepted += 1
        end
        
        # Adaptive step size
        if step % 100 == 0
            current_acceptance = accepted / 100
            step_size *= (current_acceptance / mc_params.target_acceptance)
            accepted = 0
        end
        
        # Sample after equilibration
        if step > mc_params.equilibration_steps && step % mc_params.sampling_interval == 0
            push!(energy_history, current_energy)
        end
    end
    
    return (
        final_config = (positions=positions, moments=moments),
        energy_history = energy_history,
        acceptance_rate = accepted/mc_params.n_steps
    )
end

end  # module MonteCarlo
