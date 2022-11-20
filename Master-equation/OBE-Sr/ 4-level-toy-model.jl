using DifferentialEquations
using Plots: plot, plot!, png, display
using BenchmarkTools

function Lindblad_21(den)
    Gamma_21 = 2*pi*32 # MHz, decay rate

    return Gamma_21*[den[2, 2]    -den[1, 2]/2  0             0
                    -den[2, 1]/2  -den[2, 2]    -den[2, 3]/2  -den[2, 4]/2
                    0             -den[3, 2]/2  0             0
                    0             -den[4, 2]/2  0             0]
end

function Lindblad_32(den)
    Gamma_32 = 2*pi*3 # MHz, decay rate

    return Gamma_32*[0            0             -den[1, 3]/2   0
                    0             den[3, 3]     -den[2, 3]/2   0
                    -den[3, 1]/2  -den[3, 2]/2  -den[3, 3]     -den[3, 4]/2
                    0             0             -den[4, 3]/2   0]
end

function Lindblad_34(den)
    Gamma_34 = 2*pi*0.5 # MHz, decay rate

    return Gamma_34*[0            0             -den[1, 3]/2  0
                    0             0             -den[2, 3]/2  0
                    -den[3, 1]/2  -den[3, 2]/2  -den[3, 3]    -den[3, 4]/2
                    0             0             -den[4, 3]/2  -den[3, 3]]
end

function Hamiltonian()
    Omega_12 = 10 # MHz, Rabi freq
    Omega_23 = 3 # MHz, Rabi freq

    Delta_12 = 0 # MHz, detuning
    Delta_23 = 0 # MHz, detuning

    return [0           Omega_12/2  0                   0
            Omega_12/2  -Delta_12   Omega_23/2          0
            0           Omega_23/2  -Delta_12-Delta_23  0
            0           0           0                   0]
end

function MasterEquation!(dden, den, p, t)
    d = -im * (Hamiltonian()*den - den*Hamiltonian())
    d += Lindblad_21(den)
    d += Lindblad_32(den)
    d += Lindblad_34(den)

    for i in 1:4
        for j in 1:4
            dden[i, j] = d[i, j]
        end
    end

end

function SolveME(;tspan=(0, 1000.0), savestep=5)
    den0 = [1.0+0*im 0 0 0
            0      0 0 0
            0      0 0 0
            0      0 0 0]
    prob = ODEProblem(MasterEquation!, den0, tspan, saveat=savestep)
    sol = solve(prob)

    return sol
end

begin
    @time sol = SolveME(tspan=(0, 10000))

    pop_1 = [abs(sol.u[i][1, 1])^2 for i in 1:length(sol.t)]
    pop_2 = [abs(sol.u[i][2, 2])^2 for i in 1:length(sol.t)]
    pop_3 = [abs(sol.u[i][3, 3])^2 for i in 1:length(sol.t)]
    pop_4 = [abs(sol.u[i][4, 4])^2 for i in 1:length(sol.t)]

    fig = plot(sol.t, pop_1, label="state 1 population")
    fig = plot!(sol.t, pop_2, label="state 2 population")
    fig = plot!(sol.t, pop_3, label="state 3 population")
    fig = plot!(sol.t, pop_4, label="state 4 population")

    display(fig)
end