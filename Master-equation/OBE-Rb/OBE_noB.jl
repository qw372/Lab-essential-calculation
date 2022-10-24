using DifferentialEquations
using Plots: plot, plot!, png, display
using WignerSymbols: wigner3j, wigner6j

# follow https://iopscience.iop.org/article/10.1088/1367-2630/18/12/123017
function solve_obe_zeeman(Ω, Δ, Fprime, F, u0; tspan=(0, 20), stepsize=0.02)
    # all frequencies are in unit of Γ, and time is in unit of 1/Γ
    # Ω = 1
    # δ = -1

    function ε(q)
        q == 1 && return 1
        q == -1 && return 0
        q == 0 && return 0

        error("undefined component of electric field.")
    end

    function wigner3j_modified(j1, j2, j3, m1, m2, m3)
        abs(m1) > j1 && return 0
        abs(m3) > j3 && return 0
        return wigner3j(Float64, j1, j2, j3, m1, m2, m3)
    end


    function Ωpartial(q, n)
        return im*Ω/2*ε(-q)*((-1)^(F-n))*wigner3j_modified(F, 1, Fprime, -n, -q, n+q)
    end

    function Γpartial(q, m, n)
        return (2Fprime+1)*((-1)^(-m-n))*wigner3j_modified(F, 1, Fprime, -m, -q, m+q)*wigner3j_modified(F, 1, Fprime, -n, -q, n+q)
    end

    function getu(i1, i2, i3, i4, i5, i6, u)
        uu = try
            u[i1, i2, i3, i4, i5, i6]
        catch ex
            if isa(ex, BoundsError)
                0
            else
                error(ex)
            end
        end

        return uu
    end

    # elements in du and u are in order of ρ_gg, ρ_ee, ρ_ge, ρ_eg
    function obe!(du, u, p, t)
        # u[1, 2, 3, 4, 5, 6]
        # first  and fourth index: 1 or 2, indicate ground or excited states
        # second and fifth index: 1 to (max) number of hyperfine manifolds in ground or excited states
        # third and sixth index: 1 to (max) number of Zeeman states in a hyperfine manifold

        for m in -F:F, n in -F:F
            du[1, 1, (F+1)+m, 1, 1, (F+1)+n] = 0

            for q in -1:1
                du[1, 1, (F+1)+m, 1, 1, (F+1)+n] += (-conj(Ωpartial(q, n))*getu(1, 1, (F+1)+m, 2, 1, (Fprime+1)+n+q, u)
                                                    - Ωpartial(q, m)*getu(2, 1, (Fprime+1)+m+q, 1, 1, (F+1)+n, u)
                                                    + Γpartial(q, m, n)*getu(2, 1, (Fprime+1)+m+q, 2, 1, (Fprime+1)+n+q, u))
            end
        end

        for m in -Fprime:Fprime, n in -Fprime:Fprime
            du[2, 1, (Fprime+1)+m, 2, 1, (Fprime+1)+n] = 0

            for q in -1:1
                du[2, 1, (Fprime+1)+m, 2, 1, (Fprime+1)+n] += (conj(Ωpartial(q, m-q))*getu(1, 1, (F+1)+m-q, 2, 1, (Fprime+1)+n, u)
                                                                + Ωpartial(q, n-q)*getu(2, 1, (Fprime+1)+m, 1, 1, (F+1)+n-q, u))
            end
            du[2, 1, (Fprime+1)+m, 2, 1, (Fprime+1)+n] += -u[2, 1, (Fprime+1)+m, 2, 1, (Fprime+1)+n]
        end

        for m in -F:F, n in -Fprime:Fprime
            du[1, 1, (F+1)+m, 2, 1, (Fprime+1)+n] = 0

            for q in -1:1
                du[1, 1, (F+1)+m, 2, 1, (Fprime+1)+n] += (Ωpartial(q, n-q)*getu(1, 1, (F+1)+m, 1, 1, (F+1)+n-q, u)
                                                            - Ωpartial(q, m)*getu(2, 1, (Fprime+1)+m+q, 2, 1, (Fprime+1)+n, u))
            end
            du[1, 1, (F+1)+m, 2, 1, (Fprime+1)+n] += (im*Δ-1/2)*u[1, 1, (F+1)+m, 2, 1, (Fprime+1)+n]
        end

        for m in -Fprime:Fprime, n in -F:F
            du[2, 1, (Fprime+1)+m, 1, 1, (F+1)+n] = conj(du[1, 1, (F+1)+n, 2, 1, (Fprime+1)+m])
        end

    end

    prob = ODEProblem(obe!, u0, tspan, saveat=stepsize)
    sol = solve(prob, Tsit5())

    return sol
end

begin
    Ω = 5
    Δ = -1.0
    F = 2
    Fprime = 3
    tspan = (0.0, 600.0)
    u0 = zeros(Complex{Float64}, (2, 1, 2*max(F, Fprime)+1, 2, 1, 2*max(F, Fprime)+1))
    for m in 1:2F+1
        u0[1, 1, m, 1, 1, m] = 1/(2F+1)+0im
    end
    @time sol = solve_obe_zeeman(Ω, Δ, Fprime, F, u0, tspan=tspan, stepsize=0.02)
    matele = [[real(sol.u[i][1, 1, m, 1, 1, m]) for i in 1:length(sol.t)] for m in 1:2F+1]
    groundpopulation = [sum([real(sol.u[i][1, 1, m, 1, 1, m]) for m in 1:2F+1]) for i in 1:length(sol.t)]
    fig = plot(sol.t, matele, labels=reshape(["ground mF=$i" for i in -F:F], (1, 2F+1)), legend=:right)
    # plot!(fig, sol.t, groundpopulation, label="population in all ground states")
    display(fig)
    png("latest.png")
end
