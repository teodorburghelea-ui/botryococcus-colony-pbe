#!/usr/bin/env julia

# Numerical-integrity checks for the structured population-balance equation
# in manuscript Eq. (3).  This file intentionally uses only the Julia
# standard library so that the tests can be rerun without a parameter fit or
# a package environment.

using Printf
using Statistics

const DEFAULT_PHI = 0.86
const BIO_FRACTIONS = (0.42, 0.58)
const HYD_FRACTIONS = (0.22, 0.78)

struct XGrid
    edges::Vector{Float64}
    centers::Vector{Float64}
    widths::Vector{Float64}
end

struct AGrid
    edges::Vector{Float64}
    centers::Vector{Float64}
    widths::Vector{Float64}
end

"""A log-spaced sectional grid for colony material x."""
function x_grid_from_diameter(dmin::Float64, dmax::Float64, nx::Int;
                              phi::Float64 = DEFAULT_PHI)
    xmin = pi * phi * dmin^3 / 6
    xmax = pi * phi * dmax^3 / 6
    edges = exp.(collect(range(log(xmin), log(xmax), length = nx + 1)))
    centers = sqrt.(edges[1:end-1] .* edges[2:end])
    return XGrid(edges, centers, diff(edges))
end

"""A uniform finite-volume grid for the bounded internal coordinate a."""
function a_grid(na::Int)
    edges = collect(range(0.0, 1.0, length = na + 1))
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
    return AGrid(edges, centers, diff(edges))
end

diameter(x::Float64; phi::Float64 = DEFAULT_PHI) = (6x / (pi * phi))^(1 / 3)

Base.@kwdef struct ModelParameters
    gamma::Float64 = 0.0                 # G_x = gamma x
    state_relaxation::Float64 = 0.0       # G_a = relaxation (a_eq - a)
    state_equilibrium::Float64 = 0.50

    bio_rate::Float64 = 0.0
    bio_dcrit::Float64 = 2.0
    bio_dwidth::Float64 = 0.30
    bio_acrit::Float64 = 0.50
    bio_awidth::Float64 = 0.08

    hyd_rate::Float64 = 0.0
    stress_prefactor::Float64 = 1.0
    cohesion0::Float64 = 1.0
    cohesion_state_slope::Float64 = 0.0
    hyd_exponent::Float64 = 1.0

    death_rate::Float64 = 0.0
    death_size_slope::Float64 = 0.0
    dilution::Float64 = 0.0
    shell_loss::Float64 = 0.0
    aggregation_coefficient::Float64 = 0.0
end

logistic(z::Float64) = z >= 0 ? 1 / (1 + exp(-z)) : exp(z) / (1 + exp(z))

function beta_bio(x::Float64, a::Float64, p::ModelParameters)
    d = diameter(x)
    return p.bio_rate * logistic((d - p.bio_dcrit) / p.bio_dwidth) *
           logistic((a - p.bio_acrit) / p.bio_awidth)
end

function beta_hyd(x::Float64, a::Float64, p::ModelParameters)
    d = diameter(x)
    cohesion = p.cohesion0 * (1 + p.cohesion_state_slope * a)
    excess = max(0.0, p.stress_prefactor * d / cohesion - 1)
    return p.hyd_rate * excess^p.hyd_exponent
end

death_coefficient(x::Float64, p::ModelParameters) =
    p.death_rate * (1 + p.death_size_slope * diameter(x) / (1 + diameter(x)))

growth_x(x::Float64, p::ModelParameters) = p.gamma * x
growth_a(a::Float64, p::ModelParameters) = p.state_relaxation * (p.state_equilibrium - a)

"""Return the two bracketing sectional centers and linear weights.

The x interpolation preserves both daughter number and the first material
moment exactly.  Leaving the chosen grid is an explicit error rather than a
silent numerical material loss.
"""
function brackets_and_weights(c::Vector{Float64}, value::Float64;
                              quantity::String = "coordinate")
    tolerance = 128 * eps(max(abs(value), maximum(abs, c)))
    if value < c[1] - tolerance || value > c[end] + tolerance
        error("$quantity target $value lies outside the sectional grid")
    end
    value = clamp(value, c[1], c[end])
    if value <= c[1]
        return 1, 1, 1.0, 0.0
    elseif value >= c[end]
        n = length(c)
        return n, n, 1.0, 0.0
    end
    left = searchsortedlast(c, value)
    right = left + 1
    wright = (value - c[left]) / (c[right] - c[left])
    return left, right, 1 - wright, wright
end

inside_sectional_domain(c::Vector{Float64}, value::Float64) =
    c[1] <= value <= c[end]

"""Deposit a sectional source in x and a using positive linear weights."""
function deposit!(du::Matrix{Float64}, xg::XGrid, ag::AGrid,
                  x::Float64, a::Float64, source::Float64)
    source == 0 && return
    il, ir, wl, wr = brackets_and_weights(xg.centers, x; quantity = "material")
    kl, kr, vl, vr = brackets_and_weights(ag.centers, a; quantity = "state")
    du[il, kl] += source * wl * vl
    du[il, kr] += source * wl * vr
    du[ir, kl] += source * wr * vl
    du[ir, kr] += source * wr * vr
    return nothing
end

"""First-order upwind finite-volume discretisation of ∂_x(G_x n)."""
function x_transport_rhs!(du::Matrix{Float64}, u::Matrix{Float64},
                          xg::XGrid, p::ModelParameters)
    nx, na = size(u)
    material_rate = 0.0
    # The two outer face fluxes are zero: a sufficiently wide grid makes the
    # no-flux numerical boundary independent of the test solution.
    for i in 1:(nx - 1)
        gface = growth_x(xg.edges[i + 1], p)
        for k in 1:na
            flux = gface >= 0 ? gface * u[i, k] / xg.widths[i] :
                                gface * u[i + 1, k] / xg.widths[i + 1]
            du[i, k] -= flux
            du[i + 1, k] += flux
            material_rate += flux * (xg.centers[i + 1] - xg.centers[i])
        end
    end
    return material_rate
end

"""First-order upwind finite-volume discretisation of ∂_a(G_a n)."""
function a_transport_rhs!(du::Matrix{Float64}, u::Matrix{Float64},
                          ag::AGrid, p::ModelParameters)
    nx, na = size(u)
    # Both a-boundaries are no-flux.  The state transport therefore leaves M1
    # invariant by construction, as required by Appendix B.
    for k in 1:(na - 1)
        gface = growth_a(ag.edges[k + 1], p)
        for i in 1:nx
            flux = gface >= 0 ? gface * u[i, k] / ag.widths[k] :
                                gface * u[i, k + 1] / ag.widths[k + 1]
            du[i, k] -= flux
            du[i, k + 1] += flux
        end
    end
    return nothing
end

bio_child_state(a::Float64) = 0.20 + 0.25 * a
hyd_child_state(a::Float64) = a

"""Biological and hydrodynamic binary daughter operators.

The biological fractions lose `shell_loss` of parent material; the
hydrodynamic fractions conserve it.  Daughter mapping is performed in x,
then interpreted as a diameter transformation only for reporting.
"""
function fragmentation_rhs!(du::Matrix{Float64}, u::Matrix{Float64},
                            xg::XGrid, ag::AGrid, p::ModelParameters)
    nx, na = size(u)
    bio_material_rate = 0.0
    for i in 1:nx, k in 1:na
        parent = u[i, k]
        parent == 0 && continue
        x = xg.centers[i]
        a = ag.centers[k]

        rbio = beta_bio(x, a, p) * parent
        biological_targets = all(inside_sectional_domain(xg.centers, fraction * x)
                                 for fraction in BIO_FRACTIONS) &&
                             inside_sectional_domain(ag.centers, bio_child_state(a))
        if rbio != 0 && biological_targets
            du[i, k] -= rbio
            daughter_state = bio_child_state(a)
            for fraction in BIO_FRACTIONS
                deposit!(du, xg, ag, fraction * x, daughter_state,
                         (1 - p.shell_loss) * rbio)
            end
            bio_material_rate -= p.shell_loss * rbio * x
        end

        rhyd = beta_hyd(x, a, p) * parent
        hydrodynamic_targets = all(inside_sectional_domain(xg.centers, fraction * x)
                                   for fraction in HYD_FRACTIONS) &&
                               inside_sectional_domain(ag.centers, hyd_child_state(a))
        if rhyd != 0 && hydrodynamic_targets
            du[i, k] -= rhyd
            daughter_state = hyd_child_state(a)
            for fraction in HYD_FRACTIONS
                deposit!(du, xg, ag, fraction * x, daughter_state, rhyd)
            end
        end
    end
    return bio_material_rate
end

"""Conservative constant-kernel aggregation operator.

For an unordered pair, two parents are removed and a single aggregate with
the summed material and material-weighted state is deposited.  This is a
sectional representation of the optional aggregation operator A[n].
"""
function aggregation_rhs!(du::Matrix{Float64}, u::Matrix{Float64},
                          xg::XGrid, ag::AGrid, p::ModelParameters)
    K = p.aggregation_coefficient
    K == 0 && return nothing
    nx, na = size(u)
    ncells = nx * na
    for first in 1:ncells
        i1 = ((first - 1) % nx) + 1
        k1 = ((first - 1) ÷ nx) + 1
        n1 = u[i1, k1]
        n1 == 0 && continue
        for second in first:ncells
            i2 = ((second - 1) % nx) + 1
            k2 = ((second - 1) ÷ nx) + 1
            n2 = u[i2, k2]
            n2 == 0 && continue
            event_rate = first == second ? 0.5 * K * n1 * n2 : K * n1 * n2
            event_rate == 0 && continue
            xnew = xg.centers[i1] + xg.centers[i2]
            anew = (xg.centers[i1] * ag.centers[k1] +
                    xg.centers[i2] * ag.centers[k2]) / xnew
            # Closed material boundary: an aggregation event that would
            # leave the represented domain is not applied.  This avoids a
            # hidden material sink at either material-grid boundary.
            inside_sectional_domain(xg.centers, xnew) || continue
            inside_sectional_domain(ag.centers, anew) || continue
            if first == second
                du[i1, k1] -= 2 * event_rate
            else
                du[i1, k1] -= event_rate
                du[i2, k2] -= event_rate
            end
            deposit!(du, xg, ag, xnew, anew, event_rate)
        end
    end
    return nothing
end

function material_moment(u::Matrix{Float64}, xg::XGrid)
    total = 0.0
    for i in axes(u, 1), k in axes(u, 2)
        total += xg.centers[i] * u[i, k]
    end
    return total
end

number_moment(u::Matrix{Float64}) = sum(u)

"""Return the complete semi-discrete right-hand side of manuscript Eq. (3).

`u[i,k]` stores an integrated cell count, rather than a point density.  The
returned named tuple is an independently accumulated M1 ledger for the terms
in Appendix B.
"""
function rhs!(du::Matrix{Float64}, u::Matrix{Float64}, xg::XGrid, ag::AGrid,
              p::ModelParameters, inlet::Matrix{Float64})
    fill!(du, 0.0)
    growth = x_transport_rhs!(du, u, xg, p)
    a_transport_rhs!(du, u, ag, p)
    biological = fragmentation_rhs!(du, u, xg, ag, p)
    aggregation_rhs!(du, u, xg, ag, p)

    linear = 0.0
    for i in axes(u, 1), k in axes(u, 2)
        x = xg.centers[i]
        loss = p.dilution + death_coefficient(x, p)
        du[i, k] += -loss * u[i, k] + p.dilution * inlet[i, k]
        linear += x * (-loss * u[i, k] + p.dilution * inlet[i, k])
    end
    # The conservative hydrodynamic and aggregation operators, and closed
    # internal-state transport, have zero M1 contribution.  They are tested
    # separately below rather than assumed away in the implementation.
    return (growth = growth, biological = biological, linear = linear)
end

"""CFL bound for non-negative explicit finite-volume updates."""
function positivity_timestep(u::Matrix{Float64}, xg::XGrid, ag::AGrid,
                             p::ModelParameters; cfl::Float64 = 0.75)
    nx, na = size(u)
    m0 = number_moment(u)
    max_outflow = 0.0
    for i in 1:nx, k in 1:na
        outflow = p.dilution + death_coefficient(xg.centers[i], p) +
                  beta_bio(xg.centers[i], ag.centers[k], p) +
                  beta_hyd(xg.centers[i], ag.centers[k], p) +
                  p.aggregation_coefficient * m0
        if i > 1
            outflow += max(0.0, -growth_x(xg.edges[i], p)) / xg.widths[i]
        end
        if i < nx
            outflow += max(0.0, growth_x(xg.edges[i + 1], p)) / xg.widths[i]
        end
        if k > 1
            outflow += max(0.0, -growth_a(ag.edges[k], p)) / ag.widths[k]
        end
        if k < na
            outflow += max(0.0, growth_a(ag.edges[k + 1], p)) / ag.widths[k]
        end
        max_outflow = max(max_outflow, outflow)
    end
    return max_outflow == 0 ? Inf : cfl / max_outflow
end

"""Advance with positivity-preserving forward Euler finite-volume steps."""
function evolve(u0::Matrix{Float64}, xg::XGrid, ag::AGrid,
                p::ModelParameters, inlet::Matrix{Float64};
                tend::Float64, dtmax::Float64 = 0.002)
    u = copy(u0)
    du = similar(u)
    m1_initial = material_moment(u, xg)
    m1_ledger = m1_initial
    min_cell = minimum(u)
    time = 0.0
    steps = 0
    while time < tend - 32eps(tend)
        dt = min(positivity_timestep(u, xg, ag, p), dtmax, tend - time)
        isfinite(dt) && dt > 0 || error("could not determine a positive time step")
        rates = rhs!(du, u, xg, ag, p, inlet)
        unew = u .+ dt .* du
        local_minimum = minimum(unew)
        local_minimum < -5e-13 && error("positivity violation: $local_minimum")
        # Only remove possible round-off values after the strict positivity
        # check.  A genuine negative value stops the test above.
        if local_minimum < 0
            unew .= max.(unew, 0.0)
        end
        m1_ledger += dt * (rates.growth + rates.biological + rates.linear)
        u = unew
        min_cell = min(min_cell, minimum(u))
        time += dt
        steps += 1
    end
    return (u = u, m1_initial = m1_initial, m1_final = material_moment(u, xg),
            m1_ledger = m1_ledger, min_cell = min_cell, steps = steps)
end

"""Smooth, deliberately uncalibrated initial or inlet sectional population."""
function lognormal_population(xg::XGrid, ag::AGrid;
                              dmedian::Float64 = 3.0,
                              sigma_logd::Float64 = 0.28,
                              amean::Float64 = 0.62,
                              sigma_a::Float64 = 0.12,
                              number::Float64 = 1.0)
    u = zeros(length(xg.centers), length(ag.centers))
    for i in axes(u, 1), k in axes(u, 2)
        x = xg.centers[i]
        d = diameter(x)
        p_d = exp(-0.5 * ((log(d) - log(dmedian)) / sigma_logd)^2) /
              (d * sigma_logd * sqrt(2pi))
        # p_x = p_d (dD/dx), and dD/dx = D/(3x).
        p_x = p_d * d / (3x)
        p_a = exp(-0.5 * ((ag.centers[k] - amean) / sigma_a)^2)
        u[i, k] = p_x * xg.widths[i] * p_a * ag.widths[k]
    end
    u .*= number / sum(u)
    return u
end

relative_residual(value::Float64, reference::Float64) =
    abs(value - reference) / max(abs(reference), 1e-30)

function test_row(name::String, result; analytic::Union{Nothing, Float64} = nothing)
    balance = relative_residual(result.m1_final, result.m1_ledger)
    analytic_error = analytic === nothing ? NaN : relative_residual(result.m1_final, analytic)
    return (test = name, m1_initial = result.m1_initial, m1_final = result.m1_final,
            m1_ledger = result.m1_ledger, balance_relative = balance,
            analytic_relative = analytic_error, min_cell = result.min_cell,
            steps = result.steps)
end

function write_moment_csv(path::String, rows)
    open(path, "w") do io
        println(io, "test,m1_initial,m1_final,m1_ledger,relative_balance_residual,relative_analytic_error,min_cell,steps")
        for row in rows
            analytic = isnan(row.analytic_relative) ? "" : @sprintf("%.8e", row.analytic_relative)
            @printf(io, "%s,%.12e,%.12e,%.12e,%.8e,%s,%.8e,%d\n",
                    row.test, row.m1_initial, row.m1_final, row.m1_ledger,
                    row.balance_relative, analytic, row.min_cell, row.steps)
        end
    end
end

function run_moment_tests()
    rows = NamedTuple[]
    xg = x_grid_from_diameter(0.03, 100.0, 48)
    ag = a_grid(8)
    u0 = lognormal_population(xg, ag)
    zero_inlet = zeros(size(u0))

    bio = ModelParameters(bio_rate = 0.60, bio_dcrit = 0.10, bio_dwidth = 0.10,
                          bio_acrit = 0.05, bio_awidth = 0.05)
    res = evolve(u0, xg, ag, bio, zero_inlet; tend = 0.45)
    push!(rows, test_row("biological_daughter_material", res))

    hyd = ModelParameters(hyd_rate = 0.42, stress_prefactor = 1_000.0,
                          cohesion0 = 1.0, hyd_exponent = 0.0)
    res = evolve(u0, xg, ag, hyd, zero_inlet; tend = 0.45)
    push!(rows, test_row("hydrodynamic_daughter_material", res))

    growth = ModelParameters(gamma = 0.16)
    res = evolve(u0, xg, ag, growth, zero_inlet; tend = 0.35)
    push!(rows, test_row("diameter_transport", res))

    state = ModelParameters(state_relaxation = 0.55, state_equilibrium = 0.44)
    res = evolve(u0, xg, ag, state, zero_inlet; tend = 0.45)
    push!(rows, test_row("state_transport", res))

    xga = x_grid_from_diameter(0.03, 100.0, 18)
    aga = a_grid(4)
    ua = lognormal_population(xga, aga)
    aggregation = ModelParameters(aggregation_coefficient = 0.12)
    res = evolve(ua, xga, aga, aggregation, zeros(size(ua)); tend = 0.35)
    push!(rows, test_row("aggregation_material", res))

    inlet = lognormal_population(xg, ag; dmedian = 2.0, amean = 0.38, number = 0.65)
    exchange = ModelParameters(dilution = 0.30)
    res = evolve(u0, xg, ag, exchange, inlet; tend = 0.70)
    analytic = material_moment(inlet, xg) +
               (material_moment(u0, xg) - material_moment(inlet, xg)) * exp(-exchange.dilution * 0.70)
    push!(rows, test_row("chemostat_exchange", res; analytic = analytic))

    shell = ModelParameters(bio_rate = 0.48, bio_dcrit = 0.10, bio_dwidth = 0.10,
                            bio_acrit = 0.05, bio_awidth = 0.05, shell_loss = 0.11)
    res = evolve(u0, xg, ag, shell, zero_inlet; tend = 0.45)
    push!(rows, test_row("biological_shell_loss", res))

    # A small grid keeps the explicit O(N^2) aggregation check inexpensive;
    # it nevertheless exercises every term in Eq. (3) simultaneously.
    xgf = x_grid_from_diameter(0.03, 100.0, 32)
    agf = a_grid(6)
    uf = lognormal_population(xgf, agf)
    inf = lognormal_population(xgf, agf; dmedian = 2.1, amean = 0.35, number = 0.70)
    full = ModelParameters(gamma = 0.12, state_relaxation = 0.24,
                           state_equilibrium = 0.45, bio_rate = 0.36,
                           bio_dcrit = 2.0, bio_dwidth = 0.35,
                           bio_acrit = 0.46, bio_awidth = 0.10,
                           hyd_rate = 0.09, stress_prefactor = 0.85,
                           cohesion0 = 0.80, cohesion_state_slope = 0.25,
                           hyd_exponent = 1.2, death_rate = 0.025,
                           death_size_slope = 0.20, dilution = 0.07,
                           shell_loss = 0.06, aggregation_coefficient = 0.045)
    res = evolve(uf, xgf, agf, full, inf; tend = 0.30)
    push!(rows, test_row("full_eq3", res))

    max_residual = maximum(row.balance_relative for row in rows)
    min_population = minimum(row.min_cell for row in rows)
    max_residual < 2e-11 || error("material-moment verification failed: $max_residual")
    min_population >= -5e-13 || error("positivity verification failed: $min_population")
    return rows
end

function standard_normal_cdf(z::Float64)
    # `erf` is supplied by the system C math library.  Calling it directly
    # keeps this verification script free of a Julia package dependency.
    erf_value = ccall(:erf, Cdouble, (Cdouble,), z / sqrt(2))
    return 0.5 * (1 + erf_value)
end

function diameter_transformation_convergence()
    dmin, dmax = 0.20, 40.0
    dmedian, sigma = 3.0, 0.35
    ag = a_grid(10)
    rows = NamedTuple[]
    previous_cdf = NaN
    previous_dv = NaN
    normalizer = standard_normal_cdf((log(dmax) - log(dmedian)) / sigma) -
                 standard_normal_cdf((log(dmin) - log(dmedian)) / sigma)
    # The numerical population is normalised on [dmin,dmax], so compare it
    # with the correspondingly truncated analytic material-weighted diameter,
    # E[D^4]/E[D^3], rather than with the infinite-domain lognormal moment.
    function truncated_lognormal_moment(order::Int)
        shifted_upper = (log(dmax) - log(dmedian) - order * sigma^2) / sigma
        shifted_lower = (log(dmin) - log(dmedian) - order * sigma^2) / sigma
        return exp(order * log(dmedian) + 0.5 * order^2 * sigma^2) *
               (standard_normal_cdf(shifted_upper) - standard_normal_cdf(shifted_lower))
    end
    exact_dv = truncated_lognormal_moment(4) / truncated_lognormal_moment(3)
    for nx in (32, 64, 128, 256)
        xg = x_grid_from_diameter(dmin, dmax, nx)
        u = lognormal_population(xg, ag; dmedian = dmedian, sigma_logd = sigma)
        masses = vec(sum(u, dims = 2))
        discrete_cdf = cumsum(masses)
        exact_cdf = [(standard_normal_cdf((log(diameter(xg.edges[j + 1])) - log(dmedian)) / sigma) -
                      standard_normal_cdf((log(dmin) - log(dmedian)) / sigma)) / normalizer
                     for j in 1:nx]
        cdf_error = maximum(abs.(discrete_cdf .- exact_cdf))
        dv = sum(xg.centers[i] * diameter(xg.centers[i]) * masses[i] for i in 1:nx) /
             sum(xg.centers[i] * masses[i] for i in 1:nx)
        dv_error = relative_residual(dv, exact_dv)
        cdf_order = isnan(previous_cdf) ? NaN : log(previous_cdf / cdf_error) / log(2)
        dv_order = isnan(previous_dv) ? NaN : log(previous_dv / dv_error) / log(2)
        push!(rows, (nx = nx, cdf_linf = cdf_error, volume_weighted_diameter_relative = dv_error,
                     cdf_order = cdf_order, diameter_order = dv_order))
        previous_cdf, previous_dv = cdf_error, dv_error
    end
    rows[end].cdf_linf < rows[1].cdf_linf || error("diameter transformation did not converge")
    rows[end].volume_weighted_diameter_relative < rows[1].volume_weighted_diameter_relative ||
        error("diameter transformation moment did not converge")
    return rows
end

function write_diameter_csv(path::String, rows)
    open(path, "w") do io
        println(io, "nx,cdf_linf,volume_weighted_diameter_relative,cdf_observed_order,diameter_observed_order")
        for row in rows
            cdf_order = isnan(row.cdf_order) ? "" : @sprintf("%.6f", row.cdf_order)
            d_order = isnan(row.diameter_order) ? "" : @sprintf("%.6f", row.diameter_order)
            @printf(io, "%d,%.8e,%.8e,%s,%s\n", row.nx, row.cdf_linf,
                    row.volume_weighted_diameter_relative, cdf_order, d_order)
        end
    end
end

function normalized_diameter_cdf(u::Matrix{Float64}, xg::XGrid,
                                 dquery::Vector{Float64})
    masses = vec(sum(u, dims = 2))
    masses ./= sum(masses)
    dcenters = diameter.(xg.centers)
    cumulative = cumsum(masses)
    result = zeros(length(dquery))
    for q in eachindex(dquery)
        index = searchsortedlast(dcenters, dquery[q])
        result[q] = index == 0 ? 0.0 : cumulative[index]
    end
    return result
end

function kernel_parameters(kind::Symbol)
    if kind == :biological
        return ModelParameters(bio_rate = 0.62, bio_dcrit = 0.10, bio_dwidth = 0.10,
                               bio_acrit = 0.05, bio_awidth = 0.05)
    elseif kind == :hydrodynamic
        return ModelParameters(hyd_rate = 0.62, stress_prefactor = 1_000.0,
                               cohesion0 = 1.0, hyd_exponent = 0.0)
    end
    error("unknown daughter kernel $kind")
end

function kernel_solution(kind::Symbol, nx::Int)
    xg = x_grid_from_diameter(0.03, 100.0, nx)
    ag = a_grid(10)
    u0 = lognormal_population(xg, ag; dmedian = 3.0, sigma_logd = 0.23,
                              amean = 0.63, sigma_a = 0.10)
    result = evolve(u0, xg, ag, kernel_parameters(kind), zeros(size(u0)); tend = 0.38)
    return result.u, xg
end

function daughter_kernel_convergence()
    rows = NamedTuple[]
    dquery = exp.(collect(range(log(0.035), log(80.0), length = 401)))
    for kind in (:biological, :hydrodynamic)
        uref, xref = kernel_solution(kind, 512)
        cdfref = normalized_diameter_cdf(uref, xref, dquery)
        m1ref = material_moment(uref, xref)
        previous = NaN
        first_error = NaN
        for nx in (32, 64, 128, 256)
            u, xg = kernel_solution(kind, nx)
            cdf = normalized_diameter_cdf(u, xg, dquery)
            cdf_l1 = mean(abs.(cdf .- cdfref))
            cdf_linf = maximum(abs.(cdf .- cdfref))
            m1error = relative_residual(material_moment(u, xg), m1ref)
            order = isnan(previous) ? NaN : log(previous / cdf_linf) / log(2)
            push!(rows, (kernel = String(kind), nx = nx, cdf_l1 = cdf_l1,
                         cdf_linf = cdf_linf, material_relative = m1error,
                         observed_order = order))
            isnan(first_error) && (first_error = cdf_linf)
            previous = cdf_linf
        end
        previous < first_error || error("$kind daughter-kernel calculation did not converge")
    end
    return rows
end

function write_kernel_csv(path::String, rows)
    open(path, "w") do io
        println(io, "kernel,nx,diameter_cdf_l1,diameter_cdf_linf,material_relative_to_512,observed_order")
        for row in rows
            order = isnan(row.observed_order) ? "" : @sprintf("%.6f", row.observed_order)
            @printf(io, "%s,%d,%.8e,%.8e,%.8e,%s\n", row.kernel, row.nx,
                    row.cdf_l1, row.cdf_linf, row.material_relative, order)
        end
    end
end

function write_summary(path::String, moment_rows, diameter_rows, kernel_rows)
    max_residual = maximum(row.balance_relative for row in moment_rows)
    min_population = minimum(row.min_cell for row in moment_rows)
    diameter_first, diameter_last = diameter_rows[1], diameter_rows[end]
    bio = filter(row -> row.kernel == "biological", kernel_rows)
    hyd = filter(row -> row.kernel == "hydrodynamic", kernel_rows)
    open(path, "w") do io
        println(io, "# Numerical-integrity results")
        println(io)
        println(io, "All values use fixed, deliberately uncalibrated numerical test parameters.")
        @printf(io, "The largest independent M1-ledger residual across the eight Eq. (3) tests is %.3e; the smallest cell population is %.3e.\n", max_residual, min_population)
        @printf(io, "For the diameter transformation, refinement from Nx=%d to Nx=%d reduces the CDF L∞ error from %.3e to %.3e and the volume-weighted-diameter error from %.3e to %.3e.\n",
                diameter_first.nx, diameter_last.nx, diameter_first.cdf_linf,
                diameter_last.cdf_linf, diameter_first.volume_weighted_diameter_relative,
                diameter_last.volume_weighted_diameter_relative)
        @printf(io, "For the biological daughter kernel, the CDF L∞ error relative to Nx=512 falls from %.3e (Nx=32) to %.3e (Nx=256).\n", bio[1].cdf_linf, bio[end].cdf_linf)
        @printf(io, "For the hydrodynamic daughter kernel, the CDF L∞ error relative to Nx=512 falls from %.3e (Nx=32) to %.3e (Nx=256).\n", hyd[1].cdf_linf, hyd[end].cdf_linf)
        println(io)
        println(io, "The test is intentionally a numerical-verification artifact, not a calibrated biological prediction.")
    end
end

function main()
    resultdir = joinpath(@__DIR__, "results")
    mkpath(resultdir)
    moment_rows = run_moment_tests()
    diameter_rows = diameter_transformation_convergence()
    kernel_rows = daughter_kernel_convergence()
    write_moment_csv(joinpath(resultdir, "moment_balance.csv"), moment_rows)
    write_diameter_csv(joinpath(resultdir, "diameter_transform_convergence.csv"), diameter_rows)
    write_kernel_csv(joinpath(resultdir, "daughter_kernel_convergence.csv"), kernel_rows)
    write_summary(joinpath(resultdir, "numerical_integrity_summary.md"), moment_rows,
                  diameter_rows, kernel_rows)
    println("Numerical integrity checks passed. Results written to $(resultdir)")
end

main()
