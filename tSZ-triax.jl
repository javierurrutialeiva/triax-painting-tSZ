function δαz2gnom(cosmo::Cosmology.AbstractCosmology, z2r::DataInterpolations.LinearInterpolation,
    δ, α, z; δ₀ = 0., α₀ = 0., ΔZ = 5e-6, NZ::Int = 500, match_z = true)
    
    Dₐ = angular_diameter_dist(cosmo,z)
    
    C̃ = @.cos(δ₀)*cos(δ)*cos(α - α₀) + sin(δ)*sin(δ₀)
    X̃ = @. Dₐ*(cos(δ)*sin(α - α₀))/C̃
    Ỹ = @. Dₐ*(sin(δ₀)*cos(δ)*cos(α - α₀) - cos(δ₀)*sin(δ))/C̃
  
    if match_z == false
        χ₀ = z2r(z)
        χ₁, χ₂ = z2r(z - ΔZ), z2r(z + ΔZ)
        Z̃ = collect(range(χ₁ -  χ₀, χ₂ - χ₀, length=length(δ)))u"Mpc"
    else
        vmin, vmax = (minimum(X̃) + minimum(Ỹ))/2 , (maximum(X̃) + maximum(Ỹ))/2
        Z̃ = collect(range(vmin, vmax, length = length(δ)))
    end
    X̃ = uconvert.(u"m", X̃)
    Ỹ = uconvert.(u"m", Ỹ)
    Z̃ = uconvert.(u"m", Z̃)
    return X̃, Ỹ, Z̃
end



function get_coords(δ, α, z, δ₀, α₀, 
        cosmo::Cosmology.AbstractCosmology; 
        θ = 0, φ = 0, ψ = 0, ϵ₁ = 0, ϵ₂ = 0, ΔZ = 1e-6, 
        z2r::Union{DataInterpolations.LinearInterpolation, Nothing} = nothing, full = false)
    
    q₁, q₂ = 1 - ϵ₁, 1 - ϵ₂
    a, b, c = 1.0, q₁, q₂
    R = euler_rotation_matrix(θ, φ, ψ)
    T = Diagonal([1/a^2, 1/b^2, 1/c^2])
    T_rot = R*T*transpose(R)
    T_proj = T_rot[1:2, 1:2]
    λ = sort(eigvals(T_proj))
    
    lₛ = 1 / sqrt(λ[1])
    lₘ = 1 / sqrt(λ[2])
    qₚ = lₘ / lₛ
    
    cθ, cφ, cψ, = cos(θ), cos(φ), cos(ψ)
    c2θ, c2φ, c2ψ, = cos(2*θ), cos(2*φ), cos(2*ψ)
    sθ, sφ, sψ = sin(θ), sin(φ), sin(ψ)
    s2θ, s2φ, s2ψ = sin(2*θ), sin(2*φ), sin(2*ψ)

    j = 1/2*(1/q₁^2 + 1/q₂^2) - ((sθ^2*cψ^2)*(q₁^2 + q₂^2 - 2))/(q₁^2*q₂^2) + (1/q₁^2 - 1/q₂^2)*(c2φ*(cθ^2*cψ^2 - sψ^2) - cθ*s2φ*s2ψ)
    k = 1/(4*q₁^2*q₂^2)*(2*cθ*(q₁^2 - q₂^2)*c2φ*c2ψ + (sθ^2 * (q₁^2 + q₂^2 -2 ) + (1+cθ^2)*(q₁^2 - q₂^2)*c2φ)*s2ψ)
    l = 1/2*(1/q₁^2 + 1/q₂^2) - ((sθ^2*sψ^2)*(q₁^2 + q₂^2 - 2))/(q₁^2*q₂^2) + (1/q₁^2 - 1/q₂^2)*(c2φ*(cθ^2*sψ^2 - cψ^2) - cθ*s2φ*s2ψ)
        
    l = clamp(l, 1, Inf)
    k = clamp(k, 1, Inf)
    j = clamp(j, 1, Inf)
        
    if z2r === nothing
        z2r = build_z2r_interpolator(1e-10, 5., cosmo)
    end
    
    x₁,x₂,x₃ = δαz2gnom(cosmo, z2r, δ, α, z; δ₀ =  δ₀, α₀ = α₀, ΔZ = ΔZ)
    x⃗₁, x⃗₂, x⃗₃ = vec(x₁), vec(x₂), vec(x₃)
    r⃗ = SMatrix{3,Int(length(x⃗₁))}([x⃗₁'; x⃗₂'; x⃗₃';]...)
    r⃗′ = R*r⃗
    x₁′ = r⃗′[1, :]
    x₂′ = r⃗′[2, :]
    x₃′ = r⃗′[3, :]
    
    X₁,X₂,X₃ = meshgrid(x₁′, x₂′, x₃′)
    
    ζ = @. sqrt((X₁/q₁)^2 + (X₂/q₂)^2 + X₃^2)
    if full == true
        q_p = sqrt( (j + l - sqrt((j - l)^2 + 4*k^2)) / 
            (j + l + sqrt((j - l)^2 + 4*k^2)))
        
        θϵ = atan((l - j + sqrt((j - l)^2 + 4*k^2))/(2*k))
        f = sθ^2 * ((sφ/q₁)^2 + (cφ/q₂)^2) + cθ^2
        ϵ_proj = sqrt(qₚ / (q₁ * q₂)) * f^(-3/4)
        lₚ = lₛ/(ϵ_proj * sqrt(f))
        l_los = lₛ/sqrt(f)
        X̃₁, X̃₂ = meshgrid(x₁′, x₂′)
        ξ = @. sqrt(X̃₁^2 + (X̃₂^2 / qₚ^2)*(lₛ/lₚ)^2)
    end
    if full == true
        return ζ, ξ, qₚ, θϵ, ϵ_proj, lₛ, lₚ, l_los, x₁′, x₂′, x₃′
    else
        return ζ
    end
end

function generalized_nfw(x, xc, α, β, γ)
    x̄ = x / xc
    return x̄^γ * (1 + x̄^α)^((β - γ) / α)
end


function ycompton(model::Battaglia16ThermalSZProfile{T}, δ, α, δ₀, α₀, M, z; θ_physical = true, 
        Zlos::Any = nothing, return_3D = false, full = false, θ = 0, φ = 0, ψ = 0, ϵ₁ = 0, ϵ₂ = 0, 
        z2r::Union{DataInterpolations.LinearInterpolation, Nothing} = nothing, 
        match_z = true, Δz = 1e-3) where T
    
    ζ, ξ, qₚ, θϵ, ϵ_proj, lₛ, lₚ, l_los, x, y, zlos  = get_coords(δ, α, z, δ₀, α₀,cosmo; θ = θ, φ = φ, ψ = ψ, 
                                                                z2r = z2r, ϵ₁ = ϵ₁, ϵ₂ = ϵ₂, full = true, 
                                                                ΔZ = Δz);
    factor = constants.G * M * 200 * ρ_crit(model, z) * model.f_b / 2 * 0.5176
    R_200 = R_Δ(model, M, z, 200)
    par = get_params(model, M, z)
    xζ, xξ = ustrip(ζ./(lₛ*R_200)), ustrip(ξ./(lₚ*R_200))
    
    P₀, xc, α, β, γ = par.P₀, par.xc, par.α, par.β, par.γ 
    if return_3D == false
        scale = 1e9
        x² = xξ
        zmax=1e5
        rtol=eps()
        order=9
        F₂D = [ quadgk(
                    y -> scale * generalized_nfw( sqrt(y^2 + xi), xc, α, β, γ ),
                    0.0, zmax;
                    rtol = rtol,
                    order = order
                )[1]
                for xi in x²
              ]
        if full == true
            return x, y, zlos, 2*lₚ*ϵ_proj*F₂D/scale
        end
            return 2*lₚ*ϵ_proj*F₂D/scale
    end
    F₃D = @. P₀ * ((xζ/xc)^γ) * (1  + (xζ/xc)^α)^((β - γ) / α)
    return F₃D
end


function profile_grid(model::AbstractGNFW{T}; Nr = 128, NlogM = 128, Nz = 128, zmin = 1e-3, zmax = 5.,
    logMmin = 11., logMmax = 16., minR = -10, maxR = 10., Nϵ₁ = 128, Nϵ₂ = 128, Nrot = 128)where T
    
    Rs = LinRange(minR, maxR, length = Nr)
    logMs = LinRange(logMmin, logMmax, length = NlogM)
    redshifts = LinRange(z_min, z_max, N_z)
    ϵ₁ = LinRange(0.01, 0.99, Nϵ₁)
    ϵ₂ = LinRange(0.01, 0.99, Nϵ₂)
    ϕs = LinRange(0, 2π, Nrot)
    θs = LinRange(0, 2π, Nrot)
    ψs = LinRange(0, 2π, Nrot)
    return profile_grid(model, Rs, logMs, redshifts, ϵ₁, ϵ₂, ϕs, θs, ψs)
    
end

function profile_grid(model::AbstractGNFW{T}, Rs, logMs, redshifts, ϵ₁, ϵ₂, ϕs, θs, ψs) where T

    Nr, Nz, NlogM = length(Rs), length(redshifts), length(logMs)
    Nϵ₁, Nϵ₂, Nrot = length(ϵ₁), length(ϵ₂), length(ϕs)
    
    A = zeros(T, (Nr, Nr, Nz, NlogM, Nϵ, Nϵ, Nrot, Nrot, Nrot))
    z2r_interpolator = build_z2r_interpolator(1e-5, 5., model.cosmo);
    Threads.@threads for iz in 1:Nz
      z_val = redshifts[iz]

      for iM in 1:NM
        M_val = logMs[iM]

        for ie1 in 1:Nϵ₁, ie2 in 1:Nϵ₂
          eps1 = ϵ₁[ie1]
          eps2 = ϵ₂[ie2]

          for iϕ in 1:Nϕ, iθ in 1:Nθ, iψ in 1:Nψ
            φ_val = ϕs[iϕ]
            θ_val = θs[iθ]
            ψ_val = ψs[iψ]

            block = ycompton(model, Rs, Rs, 0, 0, M_val, z_val;
                             ϵ₁ =  eps1, ϵ₂ = eps2, φ = φ_val, 
                             θ = θ_val, ψ = ψ_val,
                             z2r = z2r_interpolato)
                    
            @inbounds A[:, :, iz, iM, ie1, ie2,iθ, iψ, iϕ] = block
          end
        end

      end
    end

    return Rs, redshifts, logMs, ϵ₁, ϵ₂, ϕs, θs, ψs, A
end

function build_interpolator(model::AbstractProfile; cache_file::String="", 
                            Nr=256, pad=256, overwrite=true, verbose=true)

    if overwrite || (isfile(cache_file) == false)
        verbose && print("Building new interpolator from model.\n")
        rs, z, logM, ϵ₁, ϵ₂, ϕ, θ, ψ, prof_y  = profile_grid(model; Nr = Nr)
        if length(cache_file) > 0
            verbose && print("Saving new interpolator to $(cache_file).\n")
            save(cache_file, 
                Dict("R"=>rs, "z"=>z, "logM"=>logM, 
                     "ϵ₁"=>ϵ₁, "ϵ₂"=>ϵ₂, "ϕ"=>ϕ, "θ"=>θ, "ψ"=>ψ,
                     "prof_y"=>prof_y ))

        end
    else
        print("Found cached Battaglia profile model. Loading from disk.\n")
        model_grid = load(cache_file)
        prof_logθs, prof_redshift, prof_logMs, prof_y = model_grid["prof_logθs"], 
            model_grid["prof_redshift"], model_grid["prof_logMs"], model_grid["prof_y"]
    end

    itp = Interpolations.interpolate(log.(prof_y), BSpline(Cubic(Line(OnGrid()))))
    interp_model = scale(itp, rs, z, logM, ϵ₁, ϵ₂, ϕ, θ, ψ)
    return LogInterpolatorProfile(model, interp_model)
end


