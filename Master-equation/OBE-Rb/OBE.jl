using DifferentialEquations
using Plots: plot, plot!, png, display
using WignerSymbols: wigner3j, wigner6j
using Printf

# follow https://iopscience.iop.org/article/10.1088/1367-2630/18/12/123017
function solve_obe_zeeman(Ω, Δ, Fprime, F, u0; tspan=(0, 20), savestep=0.02)
    # all frequencies are in unit of Γ, and time is in unit of 1/Γ

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

    function μpartial(q, j, n, type)
        return -im*μ(j, type)*(-1)^(j-n+q)*sqrt(j*(j+1)*(2*j+1))*wigner3j_modified(j, 1, j, -n, -q, n+q)
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

            bb = 0
            for q in -1:1
                bb += 1/2*(-B(q)*(conj(μpartial(q, F, n-q, "ground"))*getu(1, 1, (F+1)+m, 1, 1, (F+1)+n-q, u)
                                  -conj(μpartial(q, F, m, "ground"))*getu(1, 1, (F+1)+m+q, 1, 1, (F+1)+n, u))
                            +conj(B(q))*(μpartial(q, F, n, "ground")*getu(1, 1, (F+1)+m, 1, 1, (F+1)+n+q, u)
                                         -μpartial(q, F, m-q, "ground")*getu(1, 1, (F+1)+m-q, 1, 1, (F+1)+n, u)))
            end

            du[1, 1, (F+1)+m, 1, 1, (F+1)+n] += bb
        end

        for m in -Fprime:Fprime, n in -Fprime:Fprime
            du[2, 1, (Fprime+1)+m, 2, 1, (Fprime+1)+n] = 0

            for q in -1:1
                du[2, 1, (Fprime+1)+m, 2, 1, (Fprime+1)+n] += (conj(Ωpartial(q, m-q))*getu(1, 1, (F+1)+m-q, 2, 1, (Fprime+1)+n, u)
                                                                + Ωpartial(q, n-q)*getu(2, 1, (Fprime+1)+m, 1, 1, (F+1)+n-q, u))
            end

            du[2, 1, (Fprime+1)+m, 2, 1, (Fprime+1)+n] += -u[2, 1, (Fprime+1)+m, 2, 1, (Fprime+1)+n]

            bb = 0
            for q in -1:1
                bb += 1/2*(-B(q)*(conj(μpartial(q, Fprime, n-q, "excited"))*getu(2, 1, (Fprime+1)+m, 2, 1, (Fprime+1)+n-q, u)
                                  -conj(μpartial(q, Fprime, m, "excited"))*getu(2, 1, (Fprime+1)+m+q, 2, 1, (Fprime+1)+n, u))
                            +conj(B(q))*(μpartial(q, Fprime, n, "excited")*getu(2, 1, (Fprime+1)+m, 2, 1, (Fprime+1)+n+q, u)
                                         -μpartial(q, Fprime, m-q, "excited")*getu(2, 1, (Fprime+1)+m-q, 2, 1, (Fprime+1)+n, u)))
            end

            du[2, 1, (Fprime+1)+m, 2, 1, (Fprime+1)+n] += bb
        end

        for m in -F:F, n in -Fprime:Fprime
            du[1, 1, (F+1)+m, 2, 1, (Fprime+1)+n] = 0

            for q in -1:1
                du[1, 1, (F+1)+m, 2, 1, (Fprime+1)+n] += (Ωpartial(q, n-q)*getu(1, 1, (F+1)+m, 1, 1, (F+1)+n-q, u)
                                                            - Ωpartial(q, m)*getu(2, 1, (Fprime+1)+m+q, 2, 1, (Fprime+1)+n, u))
            end

            du[1, 1, (F+1)+m, 2, 1, (Fprime+1)+n] += (im*Δ-1/2)*u[1, 1, (F+1)+m, 2, 1, (Fprime+1)+n]

            bb = 0
            for q in -1:1
                bb += 1/2*(-B(q)*(conj(μpartial(q, Fprime, n-q, "excited"))*getu(1, 1, (F+1)+m, 2, 1, (Fprime+1)+n-q, u)
                                  -conj(μpartial(q, F, m, "ground"))*getu(1, 1, (F+1)+m+q, 2, 1, (Fprime+1)+n, u))
                            +conj(B(q))*(μpartial(q, Fprime, n, "excited")*getu(1, 1, (F+1)+m, 2, 1, (Fprime+1)+n+q, u)
                                         -μpartial(q, F, m-q, "ground")*getu(1, 1, (F+1)+m-q, 2, 1, (Fprime+1)+n, u)))
            end

            du[1, 1, (F+1)+m, 2, 1, (Fprime+1)+n] += bb
        end

        for m in -Fprime:Fprime, n in -F:F
            du[2, 1, (Fprime+1)+m, 1, 1, (F+1)+n] = conj(du[1, 1, (F+1)+n, 2, 1, (Fprime+1)+m])
        end

    end

    prob = ODEProblem(obe!, u0, tspan, saveat=savestep)
    sol = solve(prob, Tsit5())

    return sol
end

begin
    # magneton including g-factor
    # https://steck.us/alkalidata/rubidium87numbers.1.6.pdf
    function μ(j, type)
        if type == "ground"
            j == 2 && return 0.70/6.065 # (MHz/Gauss)/Γ
            j == 1 && return -0.70/6.065 # (MHz/Gauss)/Γ
            error("hyperfine level F=$j doesn't exist in ground states")
        elseif type == "excited"
            j == 3 && return 0.93/6.065 # (MHz/Gauss)/Γ
            j == 2 && return 0.93/6.065 # (MHz/Gauss)/Γ
            j == 1 && return 0.93/6.065 # (MHz/Gauss)/Γ
            j == 0 && return 0.93/6.065 # (MHz/Gauss)/Γ
            error("hyperfine level F=$j doesn't exist in excited states")
        else
            error("state type not supported")
        end
    end

    # laser polarization, sum of squares has to be 1
    function ε(q)
        q == 1 && return 0
        q == -1 && return 1
        q == 0 && return 0

        error("undefined component of electric field.")
    end

    # magnetic field
    bx = 1.9
    by = -0.14
    bz = -0.99

    b1 = -1/sqrt(2)*(bx+by*im)
    bm1 = 1/sqrt(2)*(bx-by*im)

    function B(q)
        # in Gauss
        # return 0

        q == 1 && return b1
        q == -1 && return bm1
        q == 0 && return bz

        error("undefined component of magnetic field.")
    end

    F = 2
    Fprime = 3

    # laserintensity = 10e-6/(π*0.01^2) # in W/m^2, laser intensity
    # ϵ = 8.8541878128e-12 # in SI units, vacuum permittivity
    # c = 299792458
    # Efield = sqrt(2*laserintensity/c/ϵ)
    # # println(Efield)
    # ħ = 1.0545718e-34
    # electron = 1.60217662e-19 # in C, electron charge
    # Γ = 2*π*6.065e6
    # red_mat_elem_F = 4.18e-10 # in m, reduced matrix element about r
    # rabifreq = Efield*red_mat_elem_F*electron/(ħ*Γ)
    #
    # Ω = Efield*red_mat_elem_F*electron/(ħ*Γ)

    # from the line below Eq. 16 of https://iopscience.iop.org/article/10.1088/1367-2630/18/12/123017
    satintensity = 1.67 # mW/cm^2
    laserintensity = 10e-3/(π*1^2) # in mW/cm^2, laser intensity

    # from Eq. 16 of https://iopscience.iop.org/article/10.1088/1367-2630/18/12/123017
    Ω = sqrt((2Fprime+1)/2*laserintensity/satintensity)
    # print(Ω)
    Δ = -26/6.065
    tspan = (0.0, 300.0)
    u0 = zeros(Complex{Float64}, (2, 1, 2*max(F, Fprime)+1, 2, 1, 2*max(F, Fprime)+1))
    for m in 1:2F+1
        u0[1, 1, m, 1, 1, m] = 1/(2F+1)+0im
    end
    @time sol = solve_obe_zeeman(Ω, Δ, Fprime, F, u0, tspan=tspan, savestep=0.1)
    # matele = [[real(sol.u[i][1, 1, m, 1, 1, m]) for i in 1:length(sol.t)] for m in 1:2F+1]
    matele = [[real(sol.u[i][1, 1, m, 1, 1, m]) for i in 1:length(sol.t)] for m in 1:2F+1]
    groundpopulation = [sum([real(sol.u[i][1, 1, m, 1, 1, m]) for m in 1:2F+1]) for i in 1:length(sol.t)]
    expopulation = [sum([real(sol.u[i][2, 1, m, 2, 1, m]) for m in 1:2Fprime+1]) for i in 1:length(sol.t)]
    # fig = plot(sol.t, matele, labels=reshape(["ground mF=$i" for i in -F:F], (1, 2F+1)), legend=:right)
    # plot!(fig, sol.t, groundpopulation, label="population in all ground states")
    fig = plot(sol.t, expopulation, label="total population in excited states")
    display(fig)
    png("latest.png")

    expopulation_ave = sum(expopulation)/length(sol.t)
    @sprintf("Total excited state population: %.3E", expopulation_ave)
end
