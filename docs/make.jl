using Documenter, EMST

makedocs(modules=[EMST],
        doctest=true)

deploydocs(deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/lithom/EMST.git",
    julia  = "0.6.2",
    osname = "Windows")
