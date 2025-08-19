function δαz2gnom(cosmo::Cosmology.AbstractCosmology,
                  z2r::DataInterpolations.LinearInterpolation,
                  δ::AbstractVector, α::AbstractVector, z;
                  δ₀=0., α₀=0., ΔZ=5e-4, NZ::Int=100,
                  match_z=true, in_degrees=true, RotationMatrix=nothing)

    # --- angle units (radians internally) ---
    if in_degrees
        δ_unit = deg2rad.(δ)
        α_unit = deg2rad.(α)
        δ₀_unit = deg2rad(δ₀)
        α₀_unit = deg2rad(α₀)
    else
        δ_unit = δ
        α_unit = α
        δ₀_unit = δ₀
        α₀_unit = α₀
    end

    # --- redshift / comoving distance grid ---
    if match_z
        zgrid = collect(range(z-ΔZ, z+ΔZ, length=length(δ)))
    else
        zgrid = collect(range(z-ΔZ, z+ΔZ, length=NZ))
    end
    χ_grid = z2r.(zgrid) |> collect
    χ0 = z2r(z)                 # reference comoving distance
    χ_offsets = χ_grid .- χ0    # this will be used for LOS displacement

    # ensure strictly increasing χ_grid
    if any(diff(χ_grid) .<= 0)
        perm = sortperm(χ_grid)
        tmp = χ_grid[perm]
        keep = [1; findall(d -> d > 0, diff(tmp)) .+ 1]
        χ_grid = tmp[keep]
        χ_offsets = χ_grid .- χ0
    end
    Nz = length(χ_grid)

    # --- build angular grid (Ny x Nx) ---
    Ny = length(δ_unit)
    Nx = length(α_unit)
    αgrid = repeat(reshape(α_unit, 1, Nx), Ny, 1)   # (Ny, Nx)
    δgrid = repeat(reshape(δ_unit, Ny, 1), 1, Nx)   # (Ny, Nx)

    # --- unit vectors on unit sphere for every angular pixel (Ny*Nx x 3) ---
    ux = vec(cos.(δgrid) .* cos.(αgrid))   # (Nvec,)
    uy = vec(cos.(δgrid) .* sin.(αgrid))
    uz = vec(sin.(δgrid))
    Nvec = length(ux)
    U = hcat(ux, uy, uz)                   # (Nvec, 3)

    # --- center unit vector (u0) and tangent basis at center ---
    u0 = [cos(δ₀_unit)*cos(α₀_unit), cos(δ₀_unit)*sin(α₀_unit), sin(δ₀_unit)]  # (3,)
    # unit vector in Dec direction (∂/∂δ, evaluated at center) -> already unit length
    e_dec = [-sin(δ₀_unit)*cos(α₀_unit), -sin(δ₀_unit)*sin(α₀_unit), cos(δ₀_unit)]
    # unit vector in RA direction (∂/∂α)/|∂/∂α| evaluated at center:
    # ∂/∂α = [-cosδ sinα, cosδ cosα, 0] ; its norm = cosδ0 -> unit vector:
    e_ra = [-sin(α₀_unit), cos(α₀_unit), 0.0]  # already unit length

    # sanity: ensure e_dec and e_ra are orthonormal with u0 (they should be)
    # (optional debug)
    # println("dots: u0·e_dec = ", dot(u0,e_dec), " u0·e_ra = ", dot(u0,e_ra), " e_dec·e_ra = ", dot(e_dec,e_ra))

    # Preallocate output arrays (Nz,Ny,Nx)
    Xbox = Array{Float64}(undef, Nz, Ny, Nx)  # dec direction
    Ybox = Array{Float64}(undef, Nz, Ny, Nx)  # ra direction
    Zbox = Array{Float64}(undef, Nz, Ny, Nx)  # los displacement

    # center position vector (3,)
    r_center = χ0 .* u0   # 3-vector

    # For each χ slice compute r = χ * U (Nvec x 3), then d = r - r_center
    # project transverse part onto e_dec and e_ra, and compute LOS displacement
    for (k, χk) in enumerate(χ_grid)
        # r_k: (Nvec, 3)
        r_k = χk .* U                     # broadcasting scalar times each row
        # difference to center
        d = r_k .- r_center'              # r_center' makes (1,3) so broadcast -> (Nvec,3)
        # remove radial component w.r.t. center direction (project out along u0)
        d_dot_u0 = d * u0                 # (Nvec,)  (matrix-vector)
        # d_t = d - (d·u0) * u0  -> (Nvec,3)
        d_t = d .- d_dot_u0 .* (u0')      # broadcast; u0' is (1,3)
        # projections
        s_dec = d_t * e_dec               # (Nvec,)
        s_ra  = d_t * e_ra                # (Nvec,)
        # LOS displacement along central LOS (signed)
        los = (r_k * u0) .- χ0            # (Nvec,)  (dot(r_k,u0) - χ0)
        # reshape back to (Ny,Nx) and store in slice k
        Xbox[k, :, :] .= reshape(s_dec, Ny, Nx)
        Ybox[k, :, :] .= reshape(s_ra, Ny, Nx)
        Zbox[k, :, :] .= reshape(los, Ny, Nx)
    end

    # optional rotation of the 3D vectors (apply same rotation to the vector components)
    if RotationMatrix !== nothing
        # RotationMatrix expected 3x3; we rotate each vector (s_dec,s_ra,los) written in (dec,ra,los) basis
        # But RotationMatrix should be expressed in same coordinate basis; if intended for (X,Y,Z) physical axes,
        # apply rotation to the 3-vector [s_dec, s_ra, los] for each element.
        # We'll apply it elementwise:
        for k in 1:Nz, j in 1:Ny, i in 1:Nx
            v = RotationMatrix * [Xbox[k,j,i], Ybox[k,j,i], Zbox[k,j,i]]
            Xbox[k,j,i] = v[1]; Ybox[k,j,i] = v[2]; Zbox[k,j,i] = v[3]
        end
    end

    # convert Mpc -> meters and attach unit
    Mpc_to_m = 3.085677581491367e22
    X̃ = (Xbox .* Mpc_to_m) .* u"m"
    Ỹ = (Ybox .* Mpc_to_m) .* u"m"
    Z̃ = (Zbox .* Mpc_to_m) .* u"m"

    # diagnostics (print min/max of physical boxes)
    println("X (dec) range [m]: ", X̃[1], " -> ", X̃[end])
    println("Y (ra)  range [m]: ", Ỹ[1], " -> ", Ỹ[end])
    println("Z (los) range [m]: ", Z̃[1], " -> ", Z̃[end])

    return X̃, Ỹ, Z̃
end



function get_coords(δ, α, z, δ₀, α₀, 
        cosmo::Cosmology.AbstractCosmology; 
        θ = 0, φ = 0, ψ = 0, ϵ₁ = 0, ϵ₂ = 0, ΔZ = 5e-4, 
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
    
    X₁,X₂,X₃ = δαz2gnom(cosmo, z2r, δ, α, z; δ₀ =  δ₀, α₀ = α₀, ΔZ = ΔZ, RotationMatrix = R)
    
    Nz, Ny, Nx = size(X₁)
    
    ζ = @. sqrt((X₁/q₁)^2 + (X₂/q₂)^2 + X₃^2)
    
    kmid = fld(Nz + 1, 2)
    X̃₁ = Array(X₁[kmid, :, :])
    X̃₂ = Array(X₂[kmid, :, :])
    x₁ = vec(X₁[kmid, :, 1])
    x₂ = vec(X₂[kmid, 1, :]) 
    x₃ = vec(X₃[:, 1, 1])
    if full == true
        q_p = sqrt( (j + l - sqrt((j - l)^2 + 4*k^2)) / 
            (j + l + sqrt((j - l)^2 + 4*k^2)))
        
        θϵ = atan((l - j + sqrt((j - l)^2 + 4*k^2))/(2*k))
        f = sθ^2 * ((sφ/q₁)^2 + (cφ/q₂)^2) + cθ^2
        ϵ_proj = sqrt(qₚ / (q₁ * q₂)) * f^(-3/4)
        lₚ = lₛ/(ϵ_proj * sqrt(f))
        l_los = lₛ/sqrt(f)
        ξ = @. sqrt(X̃₁^2 + (X̃₂^2 / qₚ^2)*(lₛ/lₚ)^2)
    end
    if full == true
        return ζ, ξ, qₚ, θϵ, ϵ_proj, f, lₛ, lₚ, l_los, x₁, x₂, x₃
    else
        return ζ
    end
end

function generalized_nfw(x, xc, α, β, γ, R_200)
    x̄ = x / xc
    return @. x̄^γ * (1 + x̄^α)^((β - γ) / α)
end

function ycompton(model::Battaglia16ThermalSZProfile{T}, δ, α, δ₀, α₀, M, z; θ_physical = true, 
        Zlos::Any = nothing, return_3D = false, full = false, θ = 0, φ = 0, ψ = 0, ϵ₁ = 0, ϵ₂ = 0, 
        z2r::Union{DataInterpolations.LinearInterpolation, Nothing} = nothing, 
        match_z = true, Δz = 5e-4) where T
    
    ζ, ξ, qₚ, θϵ, ϵ_proj, f, lₛ, lₚ, l_los, x, y, zlos  = get_coords(δ, α, z, δ₀, α₀,cosmo; θ = θ, φ = φ, ψ = ψ, 
                                                                z2r = z2r, ϵ₁ = ϵ₁, ϵ₂ = ϵ₂, full = true, 
                                                                ΔZ = Δz);
    R_200 = R_Δ(model, M, z, 200)
    factor = constants.G * M * 200 * ρ_crit(model, z) * model.f_b / 2 * 0.5176 * P_e_factor
    par = get_params(model, M, z)
    P₀, xc, α, β, γ = par.P₀, par.xc, par.α, par.β, par.γ 
    if return_3D == false
        scale = 1e9
        x² = ustrip.((ξ/R_200).^2)
        zmax=1e5
        rtol=eps()
        order=9
        F₂D = [ quadgk(
                    ζ -> scale * P₀ * generalized_nfw(√(ζ^2 +xi), xc, α, β, γ),
                    0, zmax;
                    rtol = rtol,
                    order = order
                )[1]
                for xi in x²
              ]
        F₂D *= factor 
        if full == true
            return x, y, zlos, 2/sqrt(f) * F₂D/scale
        end
            return 2*factor/sqrt(f)*F₂D/scale
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


