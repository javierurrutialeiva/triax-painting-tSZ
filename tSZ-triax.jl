struct Battaglia16ThermalSZProfile{T,C} <: AbstractGNFW{T}
    f_b::T
    cosmo::C
end

struct RotationEulerMatrix{T} <: AbstractRotationMatrix{T}
    ϕ::T
    θ::T
    ψ::T
    R::Matrix{T} 
end

function Battaglia16ThermalSZProfile(; Omega_c::T=0.2589, Omega_b::T=0.0486, h::T=0.6774) where {T <: Real}
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(T, h=h, OmegaM=OmegaM)
    return Battaglia16ThermalSZProfile(f_b, cosmo)
end

function get_params(::AbstractGNFW{T}, M_200, z) where T
	z₁ = z + 1
	m = M_200 / (1e14M_sun)
	P₀ = 18.1 * m^0.154 * z₁^-0.758
	xc = 0.497 * m^-0.00865 * z₁^0.731
	β = 4.35 * m^0.0393 * z₁^0.415
	α = 1
    γ = -0.3
    β = γ - α * β  # Sigurd's conversion from Battaglia to standard NFW
    return (xc=T(xc), α=T(α), β=T(β), γ=T(γ), P₀=T(P₀))
end

function generalized_nfw(x, xc, α, β, γ)
    x̄ = x / xc
    return x̄^γ * (1 + x̄^α)^((β - γ) / α)
end

function _generalized_scaled_nfw(x̄, α, β, γ)
    return x̄^γ * (1 + x̄^α)^((β - γ) / α)
end

function _nfw_profile_los_quadrature(x, xc, α, β, γ; zmax=1e5, rtol=eps(), order=9)
    x² = x^2
    scale = 1e9
    integral, err = quadgk(y -> scale * generalized_nfw(√(y^2 + x²), xc, α, β, γ),
                      0.0, zmax, rtol=rtol, order=order)
    return 2integral / scale
end

function pressure_profile(model::Battaglia16ThermalSZProfile{T}, r, M, z;  
        ϵ₁ = 0, ϵ₂ = 0,  ϕ = 0., θ = 0., ψ = 0., 
        z2r::Union{DataInterpolations.LinearInterpolation, Nothing} = nothing,
        ΔZ = 1e-6) where T
    δ = α = r
    cosmo = model.cosmo
    R_rot = RotationEulerMatrix(ϕ, θ, ψ)
    R, X, Y, Zlos = get_coords(δ, α, z, cosmo, R_rot, ϵ₁ = ϵ₁, ϵ₂ = ϵ₂, z2r = z2r, ΔZ = 1e-6, full = true)
    ymap = ycompton(model, R, M, z; Zlos = Zlos, return_3D = false)
    return ymap
end

function ycompton(model::Battaglia16ThermalSZProfile{T}, θ, M, z; θ_physical = true, 
        Zlos::Any = nothing, return_3D = false) where T
    factor = constants.G * M * 200 * ρ_crit(model, z) * model.f_b / 2 * 0.5176 * P_e_factor
    R_200 = R_Δ(model, M, z, 200)
    par = get_params(model, M, z)
    θ = θ_physical == true ? angular_size(model, θ, z) : θ
    Zlos = θ_physical == true ? angular_size(model, Zlos, z) : Zlos
    x = θ_physical == true ? θ / angular_size(model ,R_200, z) : θ / R_200
    if ndims(θ) == 1
        P = @. par.P₀ * _nfw_profile_los_quadrature(x, par.xc, par.α, par.β, par.γ)
        return P
    elseif ndims(θ) == 3 && Zlos !== nothing
        P₀, xc, α, β, γ = par.P₀, par.xc, par.α, par.β, par.γ 
        P3D = @. P₀ * ((x/xc)^γ) * (1  + (x/xc)^α)^((β - γ) / α)
        projected_P = trapz((:,:,Zlos), P3D)
        if return_3D == false
            return projected_P
        else
            return P3D
        end
    end
end


function profile_grid(model::AbstractGNFW{T}; Nr = 128, NlogM = 128, Nz = 128, zmin = 1e-3, zmax = 5.,
    logMmin = 11., logMmax = 16., minR = -10, maxR = 10., Nϵ = 128, Nrot = 128)where T
    
    Rs = LinRange(minR, maxR, length = Nr)
    logMs = LinRange(logMmin, logMmax, length = NlogM)
    redshifts = LinRange(z_min, z_max, N_z)
    ϵ = LinRange(0.01, 0.99, Nϵ)
    ϕs = LinRange(0, 2π, Nrot)
    θs = LinRange(0, 2π, Nrot)
    ψs = LinRange(0, 2π, Nrot)
    return profile_grid(model, Rs, logMs, redshifts, ϵ, ϕs, θs, ψs)
    
end

function profile_grid(model::AbstractGNFW{T}, Rs, logMs, redshifts, ϵ, ϕs, θs, ψs) where T

    Nr, Nz, NlogM, Nϵ, Nrot = length(Rs), length(redshifts), length(logMs), length(ϵ), length(ϕs)
    
    A = zeros(T, (Nr, Nr, Nz, NlogM, Nϵ, Nϵ, Nrot, Nrot, Nrot))
    z2r_interpolator = build_z2r_interpolator(1e-5, 5., model.cosmo);
    @threads for iz in 1:Nz
      z_val = redshifts[iz]

      for iM in 1:NM
        M_val = logMs[iM]

        for ie1 in 1:Nϵ, ie2 in 1:Nϵ
          eps1 = ϵs[ie1]
          eps2 = ϵs[ie2]

          for iϕ in 1:Nϕ, iθ in 1:Nθ, iψ in 1:Nψ
            φ_val = ϕs[iϕ]
            θ_val = θs[iθ]
            ψ_val = ψs[iψ]

            block = pressure_profile(model, Rs, M_val, z_val;
                             ϵ₁ =  eps1, ϵ₂ = eps2, ϕ = φ_val, 
                             θ = θ_val, ψ = ψ_val,
                             z2r = z2r_interpolato)
                    
            @inbounds A[:, :, iz, iM, ie1, ie2, iϕ, iθ, iψ] = block
          end
        end

      end
    end

    return Rs, redshifts, logMs, ϵs, ϕs, θs, ψs, A
end


