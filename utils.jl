abstract type AbstractProfileWorkspace{T} end
abstract type AbstractProfile{T} end
abstract type AbstractGNFW{T} <: AbstractProfile{T} end
abstract type AbstractInterpolatorProfile{T} <: AbstractProfile{T} end


struct LogInterpolatorProfile{T, P <: AbstractProfile{T}, I1, C} <: AbstractInterpolatorProfile{T}
    model::P
    itp::I1
    cosmo::C
end

function LogInterpolatorProfile(model::AbstractProfile, itp)
    return LogInterpolatorProfile(model, itp, model.cosmo)  # use wrapped cosmology
end


struct Battaglia16ThermalSZProfile{T,C} <: AbstractGNFW{T}
    f_b::T
    cosmo::C
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

function get_cosmology(::Type{T}; h=0.69,
                   Neff=3.04,
                   OmegaK=0.0,
                   OmegaM=0.29,
                   OmegaR=nothing,
                   Tcmb=2.7255,
                   w0=-1,
                   wa=0) where T

    if OmegaR === nothing
        OmegaG = 4.48131e-7*Tcmb^4/h^2
        OmegaN = Neff*OmegaG*(7/8)*(4/11)^(4/3)
        OmegaR = OmegaG + OmegaN
    end

    OmegaL = 1 - OmegaK - OmegaM - OmegaR

    if !(w0 == -1 && wa == 0)
        return Cosmology.WCDM{T}(h, OmegaK, OmegaL, OmegaM, OmegaR, w0, wa)
    end

    if OmegaK < 0
        return Cosmology.ClosedLCDM{T}(h, OmegaK, OmegaL, OmegaM, OmegaR)
    elseif OmegaK > 0
        return Cosmology.OpenLCDM{T}(h, OmegaK, OmegaL, OmegaM, OmegaR)
    else
        return Cosmology.FlatLCDM{T}(h, OmegaL, OmegaM, OmegaR)
    end
end
get_cosmology(; h=0.69, Neff=3.04, OmegaK=0.0, OmegaM=0.29, OmegaR=nothing, Tcmb=2.7255, 
    w0=-1, wa=0) = get_cosmology(Float32; h=h, Neff=Neff, OmegaK=OmegaK, 
        OmegaM=OmegaM, OmegaR=OmegaR, Tcmb=Tcmb, w0=w0, wa=wa)


function ρ_crit(model, z)
    H_z = H(model.cosmo, z)
    return uconvert(u"kg/m^3", 3H_z^2 / (8π * constants.G))
end

function R_Δ(model, M_Δ, z, Δ=200)
    return ∛(M_Δ / (4π/3 * Δ * ρ_crit(model, z)))
end

function angular_size(model::AbstractProfile{T}, physical_size, z) where T
    d_A = angular_diameter_dist(model.cosmo, z)

    # convert both to the same units and strip units for atan
    phys_siz_unitless = @. T.(ustrip.(uconvert.(unit(d_A), physical_size)))
    d_A_unitless = T(ustrip(d_A))
    return atan.(phys_siz_unitless, d_A_unitless)
end

function build_z2r_interpolator(min_z::T, max_z::T,
    cosmo::Cosmology.AbstractCosmology; n_bins=2000) where T

    zrange = LinRange(min_z, max_z, n_bins)
    rrange = zero(zrange)
    for i in 1:n_bins
        rrange[i] = ustrip(T, u"Mpc",
            Cosmology.comoving_radial_dist(u"Mpc", cosmo, zrange[i]))
    end
    z2r = DataInterpolations.LinearInterpolation(rrange, zrange);
    return z2r
end

function euler_rotation_matrix(θ, φ, ψ)
    Rzϕ = @SMatrix[
        cos(φ)  -sin(φ)   0
        sin(φ)   cos(φ)   0
           0        0     1
    ]
    Rxθ = @SMatrix[
           1        0         0
           0   cos(θ)   -sin(θ)
           0   sin(θ)    cos(θ)
    ]
    Rzψ = @SMatrix[
        cos(ψ)  -sin(ψ)   0
        sin(ψ)   cos(ψ)   0
           0        0     1
    ]

    return Rzψ * Rxθ * Rzϕ
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
