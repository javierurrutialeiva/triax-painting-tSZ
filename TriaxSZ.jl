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
using LinearAlgebra
using Statistics
using StaticArrays
using SkyCoords
using Printf
using JLD2, FileIO
using Pixell
using XGPaint
import PhysicalConstants.CODATA2018 as constants
using Roots

const M_sun = 1.98847e30u"kg"
const T_cmb =  2.725 * u"K"
const P_e_factor = constants.σ_e / (constants.m_e * constants.c_0^2)
const deg = u"°"

abstract type AbstractProfileWorkspace{T} end
abstract type AbstractRotationMatrix{T} end
abstract type AbstractProfile{T} end
abstract type AbstractGNFW{T} <: AbstractProfile{T} end
abstract type AbstractInterpolatorProfile{T} <: AbstractProfile{T} end

include("./utils.jl")
include("./profiles_y.jl")

