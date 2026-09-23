module CalciumBoneModels

using DifferentialEquations

include("lemaire.jl")
include("kroll.jl")
include("raposo.jl")
include("peterson_khurana.jl")
include("granjon.jl")

export LemaireParams,
       lemaire_rhs!,
       lemaire_problem,
       simulate_lemaire,
       KrollParams,
       kroll_aux,
       kroll_rhs!,
       kroll_problem,
       simulate_kroll,
       RaposoParams,
       raposo_initial_state,
       raposo_aux,
       raposo_rhs!,
       raposo_problem,
       simulate_raposo,
       PetersonKhuranaParams,
       peterson_khurana_initial_state,
       peterson_khurana_rhs!,
       peterson_khurana_problem,
       simulate_peterson_khurana,
       PetersonKhuranaBMDParams,
       peterson_khurana_bmd_initial_state,
       peterson_khurana_bmd_rhs!,
       peterson_khurana_bmd_problem,
       simulate_peterson_khurana_bmd,
       GRANJON_STATE_NAMES,
       GranjonParams,
       granjon_initial_state,
       granjon_aux,
       granjon_rhs!,
       granjon_problem,
       simulate_granjon

end
