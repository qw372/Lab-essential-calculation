using DifferentialEquations
using Plots: heatmap, png, display
using BenchmarkTools


function MasterEquation!(dden, den, p, t)
    Ham = p[1]
    Gamma_21 = p[2]
    Gamma_32 = p[3]
    Gamma_34 = p[4]

    d = -im * (Ham*den - den*Ham)

    d += Gamma_21*[ den[2, 2]    -den[1, 2]/2  0             0
                    -den[2, 1]/2  -den[2, 2]    -den[2, 3]/2  -den[2, 4]/2
                    0             -den[3, 2]/2  0             0
                    0             -den[4, 2]/2  0             0]

    d += Gamma_32*[ 0            0             -den[1, 3]/2   0
                    0             den[3, 3]     -den[2, 3]/2   0
                    -den[3, 1]/2  -den[3, 2]/2  -den[3, 3]     -den[3, 4]/2
                    0             0             -den[4, 3]/2   0]

    d += Gamma_34*[ 0            0             -den[1, 3]/2  0
                    0             0             -den[2, 3]/2  0
                    -den[3, 1]/2  -den[3, 2]/2  -den[3, 3]    -den[3, 4]/2
                    0             0             -den[4, 3]/2  -den[3, 3]]

    for i in 1:4
        for j in 1:4
            dden[i, j] = d[i, j]
        end
    end

end

function SolveME(Omega_12, Delta_12; tspan=(0, 1000.0), savestep=5)
    Gamma_21 = 2*pi*32 # MHz, decay rate
    Gamma_32 = 2*pi*3 # MHz, decay rate
    Gamma_34 = 2*pi*0.5 # MHz, decay rate

    Omega_23 = 20-Omega_12
    Omega_12 < 0 && return
    Omega_23 < 0 && return

    Delta_23 = -Delta_12

    Ham = [ 0           Omega_12/2  0                   0
            Omega_12/2  -Delta_12   Omega_23/2          0
            0           Omega_23/2  -Delta_12-Delta_23  0
            0           0           0                   0]

    den0 = [1.0+0*im 0 0 0
            0      0 0 0
            0      0 0 0
            0      0 0 0]

    prob = ODEProblem(MasterEquation!, den0, tspan, (Ham, Gamma_21, Gamma_32, Gamma_34), saveat=savestep)
    sol = solve(prob)

    # pop_1 = [abs(sol.u[i][1, 1])^2 for i in 1:length(sol.t)]
    # pop_2 = [abs(sol.u[i][2, 2])^2 for i in 1:length(sol.t)]
    # pop_3 = [abs(sol.u[i][3, 3])^2 for i in 1:length(sol.t)]
    pop_4 = [abs(sol.u[i][4, 4])^2 for i in 1:length(sol.t)]

    return pop_4[end]
end

function SolveMEBayesian()
    # not implemented yet
end

begin
    Omega_12_list = LinRange(0, 20, 20)
    Delta_12_list = LinRange(-10, 10, 20)
    pop_4_list = zeros(Float64, length(Omega_12_list), length(Delta_12_list))

    @time for (i, Omega_12) in enumerate(Omega_12_list)
        for (j, Delta_12) in enumerate(Delta_12_list)
            pop_4 = SolveME(Omega_12, Delta_12, tspan=(0, 1000))
            pop_4_list[i, j] = pop_4
        end
    end
    fig = heatmap(Delta_12_list, Omega_12_list, pop_4_list, xlabel="Delta_12", ylabel="Omega_12", title="julia code")
    png("2D-map-julia.png")
    display(fig)
end