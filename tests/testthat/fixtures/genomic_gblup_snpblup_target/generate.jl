# Reproducible generator for the Julia-native genomic GBLUP / SNP-BLUP target
# fixture (#49). Run from the repo root:
#
#     julia --project=. test/fixtures/genomic_gblup_snpblup_target/generate.jl
#
# No RNG is used. The supplied allele frequencies are deliberately not re-estimated
# from this tiny sample; they make the VanRaden method-1 G positive definite so the
# precision-route GBLUP can use inv(G) without ridge regularization.
using HSquared
using LinearAlgebra
using Printf

const DIR = @__DIR__

fmt(x) = @sprintf("%.17g", x)

ids = ["g1", "g2", "g3", "g4"]
marker_ids = ["m1", "m2", "m3", "m4", "m5", "m6"]
M = [
    0.0 1 2 1 0 2
    2.0 1 0 1 2 0
    1.0 0 1 2 1 1
    0.0 2 1 0 2 1
]
p = [0.30, 0.45, 0.60, 0.40, 0.55, 0.35]
y = [10.0, 12.0, 11.0, 9.0]
X = ones(length(y), 1)
Z = Matrix{Float64}(I, length(y), length(y))
sigma_g2 = 2.0
sigma_e2 = 1.0

G = genomic_relationship_matrix(M; allele_frequencies = p)
isposdef(Symmetric(G)) || error("fixture G must be positive definite")
Ginv = inv(Symmetric(G))
gblup = fit_gblup(y, X, Z, Ginv, sigma_g2, sigma_e2; ids = ids)
snp = fit_snp_blup(y, X, M, sigma_g2, sigma_e2; allele_frequencies = p, ids = marker_ids)
diff = maximum(abs.(breeding_values(gblup).values .- snp.gebv))
diff < 5e-12 || error("GBLUP/SNP-BLUP target mismatch: $diff")

function write_csv(path, rows)
    open(path, "w") do io
        for row in rows
            println(io, join(row, ","))
        end
    end
end

write_csv(joinpath(DIR, "phenotypes.csv"),
          vcat([["id", "y"]], [[ids[i], fmt(y[i])] for i in eachindex(ids)]))
write_csv(joinpath(DIR, "markers.csv"),
          vcat([[["id"]; marker_ids]], [[ids[i]; [fmt(M[i, j]) for j in axes(M, 2)]] for i in axes(M, 1)]))
write_csv(joinpath(DIR, "allele_frequencies.csv"),
          vcat([["marker", "frequency"]], [[marker_ids[j], fmt(p[j])] for j in eachindex(marker_ids)]))
write_csv(joinpath(DIR, "expected_genomic_relationship.csv"),
          vcat([[["id"]; ids]], [[ids[i]; [fmt(G[i, j]) for j in axes(G, 2)]] for i in axes(G, 1)]))
write_csv(joinpath(DIR, "expected_genomic_precision.csv"),
          vcat([[["id"]; ids]], [[ids[i]; [fmt(Ginv[i, j]) for j in axes(Ginv, 2)]] for i in axes(Ginv, 1)]))
write_csv(joinpath(DIR, "expected_beta.csv"),
          [["effect", "value"], ["Intercept", fmt(only(fixed_effects(gblup)))]])
write_csv(joinpath(DIR, "expected_gebv.csv"),
          vcat([["id", "gblup", "snp_blup"]],
               [[ids[i], fmt(breeding_values(gblup).values[i]), fmt(snp.gebv[i])] for i in eachindex(ids)]))
write_csv(joinpath(DIR, "expected_marker_effects.csv"),
          vcat([["marker", "effect"]], [[marker_ids[j], fmt(snp.marker_effects[j])] for j in eachindex(marker_ids)]))
write_csv(joinpath(DIR, "expected_metadata.csv"),
          [["key", "value"],
           ["sigma_g2", fmt(sigma_g2)],
           ["sigma_e2", fmt(sigma_e2)],
           ["k", fmt(snp.k)],
           ["n_records", string(length(ids))],
           ["n_markers", string(length(marker_ids))],
           ["method", "vanraden1_supplied_frequencies"],
           ["g_positive_definite", "true"],
           ["gblup_snp_blup_max_abs_gebv_diff", fmt(diff)]])

println("wrote fixture CSVs to ", DIR)
