module Constants
    # Fundamental constants
    export k_B, mu_0, epsilon_0, e, N_A, epsilon_r
    const k_B = 1.380649e-23      # Boltzmann constant
    const mu_0 = 4 * pi * 1e-7    # Vacuum permeability
    const epsilon_0 = 8.8541878128e-12
    const e = 1.602176634e-19     # Elementary charge
    const N_A = 6.02214076e23     # Avogadro's number
    const epsilon_r = 78.5        # Relative permittivity

    # Material properties
    export a, M_s, pHzpc, surface_site_density
    const a = 10e-9               # Particle radius
    const M_s = 4.8e5             # Saturation magnetization
    const pHzpc = 6.5             # pH at zero point charge
    const surface_site_density = 5e18 # Surface sites/m²
end
