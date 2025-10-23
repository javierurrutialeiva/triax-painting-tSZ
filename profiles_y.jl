function δαz2gnom(cosmo::Cosmology.AbstractCosmology,
                  z2r::DataInterpolations.LinearInterpolation,
                  δ::AbstractVector, α::AbstractVector, z;
                  δ₀=0., α₀=0., ΔZ=5e-4, Nz::Int=100,
                  match_z=true, in_degrees=true, RotationMatrix=nothing)
    
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
    # --- the χ are made based on the maximum physical offset along RA or DEC direction
    Δα_max = maximum(abs.(α_unit .- α₀_unit))
    Δδ_max = maximum(abs.(δ_unit .- δ₀_unit))
    χ0 = z2r(z)   # comoving distance at central z
    Δx_max = χ0 * cos(δ₀_unit) * Δα_max
    Δy_max = χ0 * Δδ_max
    max_size = max(Δx_max, Δy_max)   # maximum transverse scale in Mpc
    
    # --- define LOS comoving distance interval ---
    χ_min = χ0 - max_size
    χ_max = χ0 + max_size

    function f_inv(zz)
        zz = zz + 1e-5
        if zz === 0
            return NaN
        end
        return z2r(zz)
    end
    z_min = find_zero(zz -> f_inv(zz) - χ_min, z - 0.01)
    z_max = find_zero(zz -> f_inv(zz) - χ_max, z + 0.01)
    Nz = length(δ_unit)
    # --- define LOS redshift grid ---
    zgrid = collect(range(z_min, z_max, length=Nz))
    χ_grid = z2r.(zgrid) |> collect
    χ0 = z2r(z)     
    
    # if any(diff(χ_grid) .<= 0)
    #     perm = sortperm(χ_grid)
    #     tmp = χ_grid[perm]
    #     keep = [1; findall(d -> d > 0, diff(tmp)) .+ 1]
    #     χ_grid = tmp[keep]
    #     χ_offsets = χ_grid .- χ0
    # end
    
    Nz = length(χ_grid)

    # --- build angular grid (Ny x Nx) ---
    Ny = length(δ_unit)
    Nx = length(α_unit)
    αgrid = repeat(reshape(α_unit, 1, Nx), Ny, 1)   # (Ny, Nx)
    δgrid = repeat(reshape(δ_unit, Ny, 1), 1, Nx)   # (Ny, Nx)

    # -- It compute a grid of each angular offset
    
    # --- unit vectors on unit sphere for every angular pixel (Ny*Nx x 3) ---
    ux = vec(cos.(δgrid) .* cos.(αgrid))   # (Nvec,) unit vector pointing x
    uy = vec(cos.(δgrid) .* sin.(αgrid))   # unit vector pointing y
    uz = vec(sin.(δgrid))                  # unit vector pointing z
    Nvec = length(ux)
    
    U = hcat(ux, uy, uz)                   # (Nvec, 3), set of vectors pointing to each direction
    
    u0 = [cos(δ₀_unit)*cos(α₀_unit), cos(δ₀_unit)*sin(α₀_unit), sin(δ₀_unit)]  # (3,)
    e_dec = [-sin(δ₀_unit)*cos(α₀_unit), -sin(δ₀_unit)*sin(α₀_unit), cos(δ₀_unit)]
    e_ra = [-sin(α₀_unit), cos(α₀_unit), 0.0]  # already unit length

    Xbox = Array{Float64}(undef, Nz, Ny, Nx)  # dec direction
    Ybox = Array{Float64}(undef, Nz, Ny, Nx)  # ra direction
    Zbox = Array{Float64}(undef, Nz, Ny, Nx)  # los displacement

    r_center = χ0 .* u0   # 3-vector

    # For each χ slice compute r = χ * U (Nvec x 3), then d = r - r_center
    # project transverse part onto e_dec and e_ra, and compute LOS displacement
    for (k, χk) in enumerate(χ_grid)
        # r_k: (Nvec, 3)
        r_k = χk .* U                     # broadcasting scalar times each row
        # difference to center
        d = r_k .- r_center'              # r_center' makes (1,3) so broadcast -> (Nvec,3)
        d_dot_u0 = d * u0                 # (Nvec,)  (matrix-vector)
        # d_t = d - (d·u0) * u0  -> (Nvec,3)
        d_t = d .- d_dot_u0 .* (u0')      # broadcast; u0' is (1,3)
        # projections
        s_dec = d_t * e_dec               # (Nvec,)
        s_ra  = d_t * e_ra                # (Nvec,)
        los = (r_k * u0) .- χ0            # (Nvec,)  (dot(r_k,u0) - χ0)
        Xbox[k, :, :] .= reshape(s_dec, Ny, Nx)
        Ybox[k, :, :] .= reshape(s_ra, Ny, Nx)
        Zbox[k, :, :] .= reshape(los, Ny, Nx)
    end
    # optional rotation of the 3D vectors (apply same rotation to the vector components)
    if RotationMatrix !== nothing
        
        coords = hcat(vec(Xbox), vec(Ybox), vec(Zbox))'
        rot_coords = RotationMatrix * coords
        Xbox = reshape(rot_coords[1,:], Nx, Ny, Nz)
        Ybox = reshape(rot_coords[2,:], Nx, Ny, Nz)
        Zbox = reshape(rot_coords[3,:], Nx, Ny, Nz)
    end

    # convert Mpc -> meters and attach unit
    Mpc_to_m = 3.085677581491367e22
    X̃ = (Xbox .* Mpc_to_m) .* u"m"
    Ỹ = (Ybox .* Mpc_to_m) .* u"m"
    Z̃ = (Zbox .* Mpc_to_m) .* u"m"

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
        z2r = build_z2r_interpolator(1e-5, 10., cosmo)
    end
    
    X₁,X₂,X₃ = δαz2gnom(cosmo, z2r, δ, α, z; δ₀ =  δ₀, α₀ = α₀, ΔZ = ΔZ, RotationMatrix = R)
    
    Nz, Ny, Nx = size(X₁)
    
    ζ = @. sqrt((X₁/q₁)^2 + (X₂/q₂)^2 + X₃^2)
    
    kmid = fld(Nz + 1, 2)
    X̃₁ = Array(X₁[kmid, :, :])
    X̃₂ = Array(X₂[kmid, :, :])

    x₁ = vec(X₁[kmid, :, 1])
    x₂ = vec(X₂[kmid, 1, :]) 
    x₃ = vec(X₃[:, 1, kmid])
    if full == true
        q_jl = sqrt((j - l)^2 + 4*k^2)
        dem = clamp(j + l - q_jl, 0, Inf)
        num = clamp(j + l + q_jl, 0, Inf)
        q_p = sqrt(dem)/sqrt(num)
        
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
    
    ζ, ξ, qₚ, θϵ, ϵ_proj, f, lₛ, lₚ, l_los, x, y, zlos  = get_coords(δ, α, z, δ₀, α₀, model.cosmo; θ = θ, φ = φ, ψ = ψ, 
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
    else
        x̃ = ζ / R_200
        F₃D = @. P₀ * ((x̃/xc)^γ) * (1  + (x̃/xc)^α)^((β - γ) / α)
        if full == true
            return x, y, zlos, F₃D * factor
        else
            return F₃D * factor
        end
    end
end

function plot_example(model::Battaglia16ThermalSZProfile{T}, δ, α, δ₀, α₀, M, z, output; θ_physical = true, 
        Zlos::Any = nothing, return_3D = false, full = false, θ = 0, φ = 0, ψ = 0, ϵ₁ = 0, ϵ₂ = 0, 
        z2r::Union{DataInterpolations.LinearInterpolation, Nothing} = nothing, 
        match_z = true, Δz = 5e-4) where T
    
    x,y, zlos, Y = ycompton(model, δ, α, δ₀, α₀, M, z; θ_physical = θ_physical, Zlos = Zlos, return_3D = false,
        full = true, θ = θ, φ = φ, ψ = ψ, ϵ₁ = ϵ₁, ϵ₂ = ϵ₂, match_z = true, Δz = Δz)
    
    p = heatmap(sort(δ) .*60, sort(α) .* 60, log10.(Y))
    xlabel!("δ [arcmin]")
    ylabel!("α [arcmin]")
    
    # Build annotation text including angles
    param_text = @sprintf(
        "M = %.2e\nz = %.2f\nϵ₁ = %.2f\nϵ₂ = %.2f\nθ = %.2f\nφ = %.2f\nψ = %.2f",
        M, z, ϵ₁, ϵ₂, θ, φ, ψ
    )

    annotate!(p, minimum(δ) * 60 *0.95, maximum(α)*60 *0.6, text(param_text, :left, 10, "white"))
    
    savefig(p, output)
end



function profile_grid(model::AbstractGNFW{T}; Nr = 64, NlogM = 58, Nz = 58, zmin = 1e-2, zmax = 4.,
    logMmin = 15, logMmax = 15.7, minR = -0.05, maxR = 0.05, Nϵ₁ = 16, Nϵ₂ = 16, Nrot = 8,
    params = ["ϵ₁","ϵ₂","ϕ","θ","ψ"]) where T
    
    Rs = collect(LinRange(minR, maxR, Nr))
    logMs = collect(LinRange(logMmin, logMmax, NlogM))
    redshifts = collect(LinRange(zmin, zmax, Nz))
    ϵ₁, ϵ₂, ϕs, θs, ψs = 0., 0., 0., 0., 0.
    for i in 1:length(params)
        pname = params[i]
        if pname == "ϵ₁"
            ϵ₁ = collect(LinRange(0.01, 0.99, Nϵ₁))
        elseif pname == "ϵ₂"
            ϵ₂ = collect(LinRange(0.01, 0.99, Nϵ₂))
        elseif pname == "ϕ"
            ϕs = collect(LinRange(0, 2π, Nrot))
        elseif pname == "θ"
            θs = collect(LinRange(0, 2π, Nrot))
        elseif panem == "ψ"
            ψs = collect(LinRange(0, 2π, Nrot))
        end
    end
    return profile_grid(model, Rs, logMs, redshifts; ϵ₁ = ϵ₁, ϵ₂ = ϵ₂, 
        ϕs = ϕs, θs = θs, ψs = ψs, params = params)
    
end

function profile_grid(model::AbstractGNFW{T}, Rs, logMs, redshifts; 
        ϵ₁ = [0.0], ϵ₂ = [0.0], ϕs = [0.0], θs = [0.0], ψs = [0.0],
        params = ["ϵ₁","ϵ₂","ϕ","θ","ψ"]) where T

    all_params = ["ϵ₁","ϵ₂","ϕ","θ","ψ"]
    if !("ϵ₁" in params); ϵ₁ = [first(ϵ₁)]; end
    if !("ϵ₂" in params); ϵ₂ = [first(ϵ₂)]; end
    if !("ϕ" in params);  ϕs = [first(ϕs)]; end
    if !("θ" in params);  θs = [first(θs)]; end
    if !("ψ" in params);  ψs = [first(ψs)]; end
    
    Nr, Nz, NM = length(Rs), length(redshifts), length(logMs)
    pvals = [ϵ₁, ϵ₂, ϕs, θs,  ψs]

    plength = 1
    varying = Dict()
    fixed   = Dict()
    original_shape = [Nr, Nr, Nz, NM]
    vparams = []
    for (pname, pval) in zip(all_params, pvals)
        if pname in params && length(pval) > 1
            varying[pname] = pval
            plength*=length(pval)
            push!(original_shape, length(pval))
            push!(vparams, pval)
        else
            fixed[pname] = first(pval)
        end
    end

    dims = (Nr, Nr, Nz, NM, plength)
    println("interpolator dim = ", dims)
    println("original shape = ", original_shape)
    println("N elements = ",prod(BigInt.(dims)))
    A = zeros(Float32, dims)
    grids = collect(Iterators.product(values(varying)...))
    mat = permutedims(hcat([collect(g) for g in grids]...))
    z2r_interpolator = build_z2r_interpolator(1e-30, 10., model.cosmo);
    ncounts = 0;
    nevals = prod(dims)/(Nr*Nr)
    Threads.@threads for iz in 1:Nz
      z_val = redshifts[iz]
      for iM in 1:NM
        M_val = 10^logMs[iM]*M_sun

        for i in 1:size(mat)[1]
            pi = mat[i,:]
            p = zeros(length(all_params))
            for j in 1:length(p)
                pname = all_params[j]
                if pname in params
                        p[j] = pi[j]
                else
                        p[j] = fixed[pname]
                end
            end
            eps1, eps2, φ_val, θ_val, ψ_val = p
            block = ycompton(model, Rs, Rs, 0, 0, M_val, z_val;
                             ϵ₁ =  eps1, ϵ₂ = eps2, φ = φ_val,
                             θ = θ_val, ψ = ψ_val,
                             z2r = z2r_interpolator)  
            A[:, :, iz, iM, i] = block
            ncounts += 1
            print("\r $(ncounts) / $(nevals)")
            flush(stdout)
          end
        end
    end
    A = reshape(A, original_shape...)
    return A, Rs, redshifts, logMs, vparams
end

function build_interpolator(model::AbstractProfile; cache_file::String="", 
                            Nr=256, pad=256, overwrite=true, verbose=true,
                            params = ["ϵ₁"])

    if overwrite || (isfile(cache_file) == false)
        verbose && print("Building new interpolator from model.\n")
        prof_y, rs, z, logM, vparams = profile_grid(model; Nr = Nr, params = params)
        if length(cache_file) > 0
            verbose && print("Saving new interpolator to $(cache_file).\n")
            out = Dict("R"=>rs, "z"=>z, "logM"=>logM, "prof_y"=>prof_y)
            for i in 1:length(params)
                pname = params[i]
                vals = vparams[i]
                out[pname] = vals
            end
            save(cache_file, out)
        end
    else
        print("Found cached Battaglia profile model. Loading from disk.\n")
        model_grid = load(cache_file)
        rs = model_grid["R"]
        z = model_grid["z"]
        logM = model_grid["logM"]
        prof_y = model_grid["prof_y"]
        vparams = [model_grid[p] for p in params]
    end
     # Collect all axes into one tuple
    axes_all = (rs, rs, z, logM, vparams...)

    # Helper: check if axis is regular (uniform spacing)
    is_regular_axis(ax) = length(ax) > 1 && all(isapprox.(diff(ax), diff(ax)[1]))

    # If all axes are regular -> use scale()
    if all(is_regular_axis(ax) for ax in axes_all)
        verbose && println("Using scaled interpolation (regular axes).")
        # Convert to ranges
        axes_ranges = map(ax -> range(first(ax), stop=last(ax), length=length(ax)), axes_all)
        itp_raw = interpolate(log10.(prof_y), BSpline(Linear()))
        itp = Interpolations.scale(itp_raw, axes_ranges...)
    else
        verbose && println("Using gridded interpolation (irregular axes).")
        itp = interpolate((axes_all...), log10.(prof_y), Gridded(Linear()))
    end
    return LogInterpolatorProfile(model, itp)
end

# function create_animation(model::Battaglia16ThermalSZProfile{T}, δ, α, δ₀, α₀, M, z; 
#             pmin = 0., pmax = 0.95, N = 10, output = "anim.gif", param = "ϵ₁" , θ = 0, 
#             φ = 0, ψ = 0, ϵ₁ = 0, ϵ₂ = 0, plot_slices = false, vmin = 1e-8, vmax = 1e-5) where T
#     println("Starting animation with model : $model over param : $param")
#     pars = collect(LinRange(pmin, pmax, N))
#     gr(fmt = :png)
#     anim = @animate for val in pars
#         θ′, φ′, ψ′, ϵ₁′, ϵ₂′ = θ, φ, ψ, ϵ₁, ϵ₂
#         if param == "θ"
#             θ′ = val
#         elseif param == "φ"
#             φ′ = val
#         elseif param == "ψ"
#             ψ′ = val
#         elseif param == "ϵ₁"
#             ϵ₁′ = val
#         elseif param == "ϵ₂"
#             ϵ₂′ = val
#         else
#             error("Unknown parameter: $param")
#         end
#         println("$param : $val")
#         x, y, zlos, Y = ycompton(
#             model, δ, α, δ₀, α₀, M, z;
#             θ_physical = true,
#             full = true, θ = θ′, φ = φ′, ψ = ψ′, 
#             ϵ₁ = ϵ₁′, ϵ₂ = ϵ₂′)
#         if plot_slices == false
#             p = heatmap(sort(δ) .*60, sort(α) .* 60, log10.(Y), title = "Projected")
#             param_text = @sprintf(
#                 "M = %.2e\nz = %.2f\nϵ₁ = %.2f\nϵ₂ = %.2f\nθ = %.2f\nφ = %.2f\nψ = %.2f",
#                 M, z, ϵ₁′, ϵ₂′, θ′, φ′, ψ′
#             )
#             annotate!(p, minimum(δ) * 60 *0.95, maximum(α)*60 *0.6, text(param_text, :left, 10, "white"))
#             frame = plot(p, size = (600, 600))       
#         elseif plot_slices == true
#             x, y, zlos, Y3D = ycompton(
#                 model, δ, α, δ₀, α₀, M, z;
#                 θ_physical = true, return_3D = true,
#                 full = true, θ = θ′, φ = φ′, ψ = ψ′, 
#                 ϵ₁ = ϵ₁′, ϵ₂ = ϵ₂′)
#             m_y = div(length(δ)+1, 2)
#             m_x = div(length(α)+1, 2)
#             m_z = div(length(zlos)+1,2)
#             p1 = heatmap(log10.(Y3D[:,:, m_z]), title = "XY", xlabel = "δ [arcmin]", ylabel = "α [arcmin]"
#                 , c = :thermal, clim=(vmin, vmax))
#             p2 = heatmap(log10.(Y3D[:, m_y, :]), title = "XZ", ylabel = "los", xlabel = "x"
#                 , c = :thermal, clim=(vmin, vmax))
#             p3 = heatmap(log10.(Y3D[m_x, :, :]), title = "YZ", ylabel = "los", xlabel = "y"
#                 , c = :thermal, clim=(vmin, vmax))
        
#             layout = @layout [a b c]
#             param_text = @sprintf(
#                 "M = %.2e\nz = %.2f\nϵ₁ = %.2f\nϵ₂ = %.2f\nθ = %.2f\nφ = %.2f\nψ = %.2f",
#                 M, z, ϵ₁′, ϵ₂′, θ′, φ′, ψ′
#             )
#             annotate!(p1, minimum(δ) * 60 *0.85, maximum(α)*60 *0.75, 
#                 text(param_text, :left, 10, "white"))
#             frame = plot(p1, p2, p3; layout, size = (1200, 400))
#             end
#         end
#     println("Saving animation to $output")
#     gif(anim, output, fps=2)
# end

function profile_paint_triax!(m::Enmap{T, 2, Matrix{T}, CarClenshawCurtis{T}},
                        workspace::CarClenshawCurtisProfileWorkspace, model, Mh, z, α₀, δ₀, 
                        θmax; normalization=1, params = [0.0]) where T
    i1, j1 = sky2pix(m, deg2rad(α₀) - θmax, deg2rad(δ₀) - θmax)
    i2, j2 = sky2pix(m, deg2rad(α₀) + θmax, deg2rad(δ₀) + θmax)
    i_start = floor(Int, max(min(i1, i2), 1))
    i_stop = ceil(Int, min(max(i1, i2), size(m, 1)))
    j_start = floor(Int, max(min(j1, j2), 1))
    j_stop = ceil(Int, min(max(j1, j2), size(m, 2)))
    θmin = compute_θmin(model)
    
    ni = abs(i_stop - i_start) + 1
    nj = abs(j_stop - j_start) + 1

    α = collect(LinRange(-θmax, θmax, ni))
    δ = collect(LinRange(-θmax, θmax, nj))
    idx_j = 1
    @inbounds for j in j_start:j_stop
        idx_i = 1
        for i in i_start:i_stop
            αi, δj = α[idx_i], δ[idx_j]
            if params !== nothing
                m[i,j] += normalization * 10 .^model.itp(αi, δj, z, log10.(Mh), params...)
            else
                m[i,j] += normalization * 10 .^model.itp(αi, δj, z, log10.(Mh))
            end
            idx_i+=1
        end
        idx_j+=1
    end
end

function paint_catalog!(m::Enmap{T, 2, Matrix{T}, CarClenshawCurtis{T}},
                        workspace::CarClenshawCurtisProfileWorkspace, model, Mhs, zs, αs, δs, 
                        params = nothing) where T
    for i in 1:length(Mhs)
        Mh = Mhs[i]
        z = zs[i]
        α = αs[i]
        δ = δs[i]
        if params !== nothing
            p = params[i]
        else
            p = nothing
        end
        θmax = compute_θmax(model, Mh*M_sun, z)
        profile_paint_triax!(m, workspace, model, Mh, z, α, δ, θmax; params = p)
    end
end
