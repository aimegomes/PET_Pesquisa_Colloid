module BoundaryConditions
using ..Constants

export SystemParameters, initialize_system, calculate_zeta, calculate_debye_length,
       apply_periodic_boundary!, minimum_image_distance

"""
    SystemParameters
Struct containing all system boundary conditions and parameters.
"""
struct SystemParameters
    N::Int                     # Number of particles
    box_size::Float64          # Simulation box side length (m)
    pH::Float64                # System pH
    ionic_strength::Float64    # Ionic strength (mol/L)
    T::Float64                 # Temperature (K)
    H_ext::Vector{Float64}     # External magnetic field (A/m)
    eta::Float64               # Dynamic viscosity (Pa·s)
    h_surfactant::Float64      # Surfactant thickness (m)
    m_s::Float64               # Saturation magnetization (A·m²)

    # Keyword constructor
    function SystemParameters(;
        N::Int,
        box_size::Float64,
        pH::Float64,
        ionic_strength::Float64,
        T::Float64,
        H_ext::Vector{Float64},
        eta::Float64,
        h_surfactant::Float64,
        m_s::Float64 = M_s * (4/3)*pi*a^3  # Default calculation
    )
        new(N, box_size, pH, ionic_strength, T, H_ext, eta, h_surfactant, m_s)
    end
end

"""
    initialize_system(params::SystemParameters)
Initialize system with derived parameters based on boundary conditions.
Returns a NamedTuple of calculated properties.
"""
function initialize_system(params::SystemParameters)
    # Calculate particle volume
    V_p = (4/3)*pi*a^3
    
    # Calculate thermal energy
    k_B_T = k_B * params.T
    
    # Calculate zeta potential and Debye length
    zeta, lambda_D = calculate_zeta(params.pH, params.ionic_strength, params.T)
    
    # Calculate Hamaker constant (adjusted for medium)
    A_H = 1e-19  # Default value, can be made temperature-dependent
    
    # Calculate steric repulsion energy
    U_rep_steric = 10*k_B_T
    
    return (
        V_p = V_p,
        particle_mass = 5.3e3 * V_p,  # Cobalt ferrite density
        k_B_T = k_B_T,
        zeta = zeta,
        lambda_D = lambda_D,
        A_H = A_H,
        U_rep_steric = U_rep_steric
    )
end

"""
    calculate_debye_length(ionic_strength, T)
Calculate Debye screening length (in meters) for given ionic strength and temperature.
"""
function calculate_debye_length(ionic_strength, T)
    sqrt(epsilon_0 * epsilon_r * k_B * T / (2 * e^2 * ionic_strength * 1000 * N_A))
end

"""
    calculate_zeta(pH, ionic_strength, T)
Calculate zeta potential (in volts) using simplified surface complexation model.
"""
function calculate_zeta(pH, ionic_strength, T)
    # Surface charge density (C/m²)
    sigma_0 = surface_site_density * e * tanh((pHzpc - pH)/2)
    
    # Debye length and kappa
    lambda_D = calculate_debye_length(ionic_strength, T)
    kappa = 1/lambda_D
    
    # Zeta potential approximation (simplified)
    zeta = sigma_0 / (epsilon_0 * epsilon_r * kappa * (1 + kappa*a))
    
    # Ensure reasonable values
    zeta = clamp(zeta, -0.1, 0.1)  # Limit to ±100 mV
    
    return zeta, lambda_D
end

"""
    apply_periodic_boundary!(positions, box_size)
Apply periodic boundary conditions to particle positions.
"""
function apply_periodic_boundary!(positions, box_size)
    for pos in positions
        pos .= mod.(pos, box_size)
    end
    return positions
end

"""
    minimum_image_distance(r_i, r_j, box_size)
Calculate minimum image distance between two particles in periodic boundary conditions.
"""
function minimum_image_distance(r_i, r_j, box_size)
    dr = r_i - r_j
    dr .= mod.(dr .+ box_size/2, box_size) .- box_size/2
    return dr
end

end  # module BoundaryConditions
