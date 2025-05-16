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

function δαz2gnom(cosmo::Cosmology.AbstractCosmology, z2r::DataInterpolations.LinearInterpolation,
    δ, α, z; δ₀ = 0., α₀ = 0., ΔZ = 5e-6, NZ::Int = 500)
    χ₀ = z2r(z)
    Dₐ = angular_diameter_dist(cosmo,z)
    χ₁, χ₂ = z2r(z - ΔZ), z2r(z + ΔZ)
    Z̃ = collect(range(χ₁ -  χ₀, χ₂ - χ₀, length=length(δ)))u"Mpc"
    C̃ = @.cos(δ₀)*cos(δ)*cos(α - α₀) + sin(δ)*sin(δ₀)
    X̃ = @. Dₐ*(cos(δ)*sin(α - α₀))/C̃
    Ỹ = @. Dₐ*(sin(δ₀)*cos(δ)*cos(α - α₀) - cos(δ₀)*sin(δ))/C̃
    X̃ = uconvert.(u"m", X̃)
    Ỹ = uconvert.(u"m", Ỹ)
    Z̃ = uconvert.(u"m", Z̃)
    return X̃, Ỹ, Z̃
end

function RotationEulerMatrix(ϕ::T, θ::T, ψ::T) where T<:Real
    R_z1 = [ cos(ϕ)  -sin(ϕ)   0;
                       sin(ϕ)   cos(ϕ)   0;
                         0         0     1 ]
    R_x  = [ 1      0          0;
                      0  cos(θ)  -sin(θ);
                      0  sin(θ)   cos(θ) ]
    R_z2 = [ cos(ψ)  -sin(ψ)  0;
                      sin(ψ)   cos(ψ)  0;
                        0         0    1 ]
    R = R_z1 * R_x * R_z2
    return RotationEulerMatrix{T}(ϕ, θ, ψ, Matrix{T}(R))
end

rotation_matrix(m::RotationEulerMatrix) = m.R


function get_coords(δ, α, z, δ₀, α₀, 
        cosmo::Cosmology.AbstractCosmology; 
        ϕ = 0, θ = 0, ψ = 0, ϵ₁ = 0, ϵ₂ = 0, ΔZ = 1e-6, 
        z2r::Union{DataInterpolations.LinearInterpolation, Nothing} = nothing, full = false)
    R_rot = rotation_matrix(RotationEuelerMatrix(ϕ, θ, ψ))
    if z2r === nothing
        z2r = build_z2r_interpolator(1e-10, 5., cosmo)
    end
    x,y,z = δαz2gnom(cosmo, z2r, δ, α, z; δ₀ =  δ₀, α₀ = α₀, ΔZ = ΔZ)
    xs, ys, zs = vec(x), vec(y), vec(z)
    v⃗ = [xs'; ys'; zs']
    v⃗′ = R_rot * v⃗
    x′ = v⃗′[1, :]
    y′ = v⃗′[2, :]
    z′ = v⃗′[3, :]

    X,Y,Z = meshgrid(x′, y′, z′)  
    print(minimum(Z), "\n")
    print(maximum(Z), "\n")
    R = @. sqrt(X^2 + Y^2/((1 - ϵ₁)^2) + Z^2/((1 - ϵ₂)^2))
    if full == true
        return R, x′, y′, z′
    else
        return R
    end
end

function get_coords(δ, α, z, cosmo::Cosmology.AbstractCosmology, R_rot::RotationEulerMatrix{T}; 
        ϵ₁ = 0, ϵ₂ = 0, ΔZ = 1e-6, z2r::Union{DataInterpolations.LinearInterpolation, Nothing} = nothing, 
        full = false) where T
    R_rot = rotation_matrix(R_rot)
    if z2r === nothing
        z2r = build_z2r_interpolator(1e-10, 5., cosmo)
    end
    x,y,z = δαz2gnom(cosmo, z2r, δ, α, z; δ₀ =  0, α₀ = 0, ΔZ = ΔZ)
    xs, ys, zs = vec(x), vec(y), vec(z)
    v⃗ = [xs'; ys'; zs']
    v⃗′ = R_rot * v⃗
    x′ = v⃗′[1, :]
    y′ = v⃗′[2, :]
    z′ = v⃗′[3, :]

    X,Y,Z = meshgrid(x′, y′, z′)  
    R = @. sqrt(X^2 + Y^2/((1 - ϵ₁)^2) + Z^2/((1 - ϵ₂)^2))
    if full == true
        return R, x′, y′, z′
    else
        return R
    end
end
