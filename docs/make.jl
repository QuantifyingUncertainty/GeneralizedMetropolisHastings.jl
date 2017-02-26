using Documenter, GeneralizedMetropolisHastings

makedocs(
    modules = [GeneralizedMetropolisHastings],
    format = :html,
    clean = false,
    sitename = "GeneralizedMetropolisHastings.jl",
    authors = "Kris De Meyer",
    pages = Any[
        "Home" => "index.md",
        "Manual" => Any[
            "Guide" => "man/guide.md",
            "man/folderstructure.md",
            "man/data.md",
            "man/noise.md",
            "man/models.md",
            ],
        "Library" => Any[
            "Public" => "lib/public.md",
            "Internal" => "lib/internal.md",
            ],
        ]
)

deploydocs(
    repo = "github.com/QuantifyingUncertainty/GeneralizedMetropolisHastings.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
