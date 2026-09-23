using Test
using CalciumBoneModels
using SciMLBase

@testset "Komarova" begin
    p = KomarovaParams()
    eq = komarova_equilibrium(p)
    du = zeros(3)
    komarova_rhs!(du, [eq.x̄₁, eq.x̄₂, p.z₀], p, 0.0)
    @test du[1] ≈ 0 atol=1e-12
    @test du[2] ≈ 0 atol=1e-12
    @test du[3] ≈ 0 atol=1e-12

    sol = simulate_komarova(p; tspan=(0.0, 20.0))
    @test SciMLBase.successful_retcode(sol)
    @test all(sol[1, :] .> 0)
    @test all(sol[2, :] .> 0)
    @test sol[3, 1] ≈ p.z₀
end

@testset "Lemaire" begin
    p = LemaireParams()
    sol = simulate_lemaire(p; tspan=(0.0, 10.0))
    @test SciMLBase.successful_retcode(sol)
    @test all(sol[1, :] .> 0)
    @test all(sol[2, :] .> 0)
    @test all(sol[3, :] .> 0)
end

@testset "Pivonka" begin
    p = PivonkaParams()
    u = [6.194e-4, 5.584e-4, 8.070e-4, p.BV₀]
    a_early = pivonka_aux(u, p, 50.0)
    a_late = pivonka_aux(u, p, 150.0)
    @test a_early.D_OC_p == p.D_OC_p_early
    @test a_late.D_OC_p == p.D_OC_p_late
    @test isfinite(a_early.OPG)
    @test isfinite(a_early.RANKL)

    sol = simulate_pivonka(p; tspan=(0.0, 20.0))
    @test SciMLBase.successful_retcode(sol)
    @test all(sol[1, :] .> 0)
    @test all(sol[2, :] .> 0)
    @test all(sol[3, :] .> 0)
    @test all(isfinite, sol[4, :])
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

@testset "Peterson/Khurana local R port" begin
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
    # These values are the direct R-core derivative at the public
    # hypoparathyroidism initial state, to guard the cross-language port.
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

@testset "Nelson surgical-menopause model" begin
    p = NelsonParams()
    @test NELSON_STATE_NAMES[7] == :Bd
    @test nelson_estrogen(0.0, p; surgical=false) == 1.0
    @test nelson_estrogen(p.t_m, p; surgical=false) == 1.0
    @test nelson_estrogen(p.t_m + 365.0, p; surgical=true) < 1.0
    @test nelson_estrogen(p.t_m + 1.0e6, p; surgical=true) ≈ p.Eovx atol=1e-8

    u₀ = nelson_initial_state(p)
    @test length(u₀) == 7
    @test u₀[7] == 1.0
    du = zeros(7)
    nelson_rhs!(du, u₀, p, 0.0)
    @test maximum(abs, du[1:6]) < 1e-9
    @test isfinite(du[7])

    sol = simulate_nelson(p; tspan=(0.0, 365.0), saveat=30.0)
    @test SciMLBase.successful_retcode(sol)
    @test all(isfinite, sol)
    @test all(sol[7, :] .> 0)

    # The optional effects are dormant before the menopause time in the
    # public default parameterization, but are active after a t_m=0 override.
    p_effects = NelsonParams(t_m=0.0, eta_ovx=1.0, tau=0.01, omega_ovx=1.0)
    u₁ = nelson_initial_state(p_effects; new_effects=true)
    du₁ = zeros(7)
    nelson_rhs!(du₁, u₁, p_effects, 365.0; surgical=true, new_effects=true)
    @test isfinite(du₁[3])
    @test isfinite(du₁[5])
end
