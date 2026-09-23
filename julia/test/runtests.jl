using Test
using CalciumBoneModels
using SciMLBase

@testset "Lemaire" begin
    p = LemaireParams()
    sol = simulate_lemaire(p; tspan=(0.0, 10.0))
    @test SciMLBase.successful_retcode(sol)
    @test all(sol[1, :] .> 0)
    @test all(sol[2, :] .> 0)
    @test all(sol[3, :] .> 0)
end

@testset "Kroll" begin
    p = KrollParams()
    @test kroll_aux([10.0, 500.0, 10.0, 200.0, 1.0e-6], p, 1.0).pth_input == p.input_rate
    @test kroll_aux([10.0, 500.0, 10.0, 200.0, 1.0e-6], p, 7.0).pth_input == 0.0
    @test kroll_aux([10.0, 500.0, 10.0, 200.0, 1.0e-6],
                    KrollParams(mode=:continuous), 7.0).pth_input == 10.0

    sol = simulate_kroll(p; tspan=(0.0, 48.0))
    @test SciMLBase.successful_retcode(sol)
    @test all(sol[1, :] .> 0)
    @test all(sol[2, :] .> 0)
    @test all(sol[3, :] .> 0)
    @test all(sol[4, :] .> 0)
    @test all(sol[5, :] .>= 0)
end

@testset "Raposo" begin
    p = RaposoParams()
    u₀ = raposo_initial_state(p)
    a = raposo_aux(u₀, p)
    @test length(u₀) == 11
    @test a.J_bp_P ≈ p.Stoic_Ca_P * p.J_bp_Ca
    @test a.J_pb_Ca ≈ p.k_Ca_b * u₀[5]
    @test isfinite(a.Ca_T)

    # The public CellML defaults are not initialized at a physiological
    # equilibrium: J_pb_Ca starts at 100 while J_bp_Ca is 3.3. Preserve that
    # as-coded behavior and test it explicitly rather than hiding it.
    sol = simulate_raposo(p; tspan=(0.0, 1.0))
    @test SciMLBase.successful_retcode(sol)
    @test all(isfinite, sol)
    @test all(sol[5, :] .> 0)
    @test minimum(sol[1, :]) < 0
end

@testset "Peterson/Khurana integrated model" begin
    p = PetersonKhuranaParams()
    u₀ = peterson_khurana_initial_state(p)
    @test length(u₀) == 32
    du = zeros(32)
    p_eff = CalciumBoneModels.peterson_khurana_derived_params(p, u₀)
    peterson_khurana_rhs!(du, u₀, p_eff, 0.0)
    @test all(isfinite, du)

    sol = simulate_peterson_khurana(p; tspan=(0.0, 48.0), with_events=false)
    @test SciMLBase.successful_retcode(sol)
    @test all(isfinite, sol)
end

@testset "Granjon/CaPO4Sim" begin
    p = GranjonParams()
    u₀ = granjon_initial_state(:hypopara)
    @test length(u₀) == 22
    @test GRANJON_STATE_NAMES[5] == :Ca_p

    du = zeros(22)
    granjon_rhs!(du, u₀, p, 0.0)
    @test all(isfinite, du)
    # These values are source-artifact regression derivatives at the public
    # hypoparathyroidism initial state.
    @test du[1] ≈ 117.47668093096183 rtol=1e-12
    @test du[5] ≈ 0.091290379623257145 rtol=1e-12
    @test du[8] ≈ 0.67640804157304046 rtol=1e-12

    sol = simulate_granjon(; case=:hypopara, tspan=(0.0, 30.0), saveat=1.0)
    @test SciMLBase.successful_retcode(sol)
    @test all(isfinite, sol)
    @test all(sol[5, :] .> 0)

    # The C artifact exposes disease/intervention parameters; ensure the
    # Julia API preserves an intervention override without mutating defaults.
    p_no_pth = GranjonParams(k_prod_PTHg=0.0)
    @test p.k_prod_PTHg == 4.192
    @test p_no_pth.k_prod_PTHg == 0.0
end

@testset "Peterson/Khurana BMD extension" begin
    p = PetersonKhuranaParams()
    u₀ = peterson_khurana_bmd_initial_state(p)
    prob = peterson_khurana_bmd_problem(p; u₀, tspan=(0.0, 1.0))
    du = zeros(33)
    peterson_khurana_bmd_rhs!(du, u₀, prob.p, 0.0)
    @test length(u₀) == 33
    @test du[33] ≈ 0.0 atol=1e-14

    # The added equation is anabolic/catabolic through the source proxies.
    u_catabolic = copy(u₀)
    u_catabolic[16] = 2.0 * prob.p.CTx0
    peterson_khurana_bmd_rhs!(du, u_catabolic, prob.p, 0.0)
    @test du[33] < 0

    sol = simulate_peterson_khurana_bmd(p; tspan=(0.0, 48.0), with_events=false)
    @test SciMLBase.successful_retcode(sol)
    @test all(isfinite, sol)
    @test sol[33, 1] == 100.0
end
