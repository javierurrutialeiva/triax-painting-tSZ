using Plots 
using Interpolations
using QuadGK
using Unitful, UnitfulAstro
using Random
using SpecialFunctions
using Roots
using Trapz
using PhysicalConstants
using Cosmology
using MeshGrid
import DataInterpolations 
using Base.Threads


import PhysicalConstants.CODATA2018 as constants
const M_sun = 1.98847e30u"kg"
const T_cmb =  2.725 * u"K"
const P_e_factor = constants.σ_e / (constants.m_e * constants.c_0^2)

abstract type AbstractProfileWorkspace{T} end
abstract type AbstractRotationMatrix{T} end
abstract type AbstractProfile{T} end
abstract type AbstractGNFW{T} <: AbstractProfile{T} end
abstract type AbstractInterpolatorProfile{T} <: AbstractProfile{T} end


include("./tSZ-triax.jl")
include("./utils.jl")
