module Forces
using ..Constants
using ..BoundaryConditions
using LinearAlgebra
using Random

export total_force, calculate_steric_force, calculate_dlvo_force, calculate_magnetic_force

"""
    calculate_steric_force(d, r_hat, h_surfactant, U_rep_steric)
Calculate steric repulsion force between surfactant-coated particles.
"""
function calculate_steric_force(d, r_hat, h_surfactant, U_rep_steric)
    if d < h_surfactant
        return (U_rep_steric/h_surfactant) * exp(-d/h_surfactant) * r_hat
    else
        return zeros(3)
    end
end

"""
    calculate_dlvo_force(d, r_hat, zeta, lambda_D, a, T)
Calculate DLVO (Derjaguin-Landau-Verwey-Overbeek) force between particles.
"""
function calculate_dlvo_force(d, r_hat, zeta, lambda_D, a, T)
    if d <= 0
        return 1e6 * r_hat  # Strong repulsion when particles overlap
    end
    
    kappa = 1/lambda_D
    k_B_T = k_B * T
    
    # Electrostatic repulsion (linearized Poisson-Boltzmann)
    U_elec = 64 * pi * epsilon_0 * epsilon_r * a * (k_B_T/e)^2 * 
             tanh(e*zeta/(4*k_B_T))^2 * exp(-kappa*d)
    F_elec = -kappa * U_elec * r_hat
    
    # Van der Waals attraction
    A_H = 1e-19  # Hamaker constant (J)
    U_vdw = -A_H * a / (12*d)
    F_vdw = -(A_H * a / (12*d^2)) * r_hat
    
    return F_elec + F_vdw
end

"""
    calculate_magnetic_force(r, m1, m2)
Calculate magnetic dipole-dipole force between particles.
"""
function calculate_magnetic_force(r, m1, m2)
    r_norm = norm(r)
    r_hat = r / r_norm
    
    term1 = dot(m1, r_hat) * m2
    term2 = dot(m2, r_hat) * m1
    term3 = dot(m1, m2) * r_hat
    term4 = 5 * dot(m1, r_hat) * dot(m2, r_hat) * r_hat
    
    return (3*mu_0/(4*pi)) * (term1 + term2 + term3 - term4) / r_norm^4
end

"""
    calculate_brownian_force(k_B_T, eta, a)
Calculate Brownian (thermal) random force.
"""
function calculate_brownian_force(k_B_T, eta, a)
    sqrt(2 * k_B_T * eta) * randn(3)
end

"""
    total_force(r_i, r_j, m_i, m_j, sys_params, derived_params)
Calculate total force between two particles including:
- Steric repulsion
- DLVO forces (electrostatic + van der Waals)
- Magnetic dipole-dipole
"""
function total_force(r_i, r_j, m_i, m_j, sys_params, derived_params)
    r = BoundaryConditions.minimum_image_distance(r_i, r_j, sys_params.box_size)
    d = norm(r) - 2*a
    r_hat = r / norm(r)
    
    # Steric force
    F_steric = calculate_steric_force(d, r_hat, 
                  sys_params.h_surfactant, 
                  derived_params.U_rep_steric)
    
    # DLVO forces
    F_dlvo = calculate_dlvo_force(d, r_hat,
                derived_params.zeta,
                derived_params.lambda_D,
                a, sys_params.T)
    
    # Magnetic force
    F_mag = calculate_magnetic_force(r, m_i, m_j)
    
    return F_steric + F_dlvo + F_mag
end

"""
    external_forces(pos, mom, sys_params, derived_params)
Calculate external forces (magnetic field + Brownian)
"""
function external_forces(pos, mom, sys_params, derived_params)
    # Magnetic torque
    torque = -0.1 * cross(mom, sys_params.H_ext) / derived_params.m_s
    
    # Brownian force
    F_brownian = calculate_brownian_force(
        derived_params.k_B_T,
        sys_params.eta,
        a
    )
    
    # Drag force (would need velocity)
    F_drag = zeros(3)  # Will be calculated in dynamics
    
    return (torque=torque, F_brownian=F_brownian, F_drag=F_drag)
end

end  # module Forces
