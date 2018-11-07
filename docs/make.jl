using Documenter
using ExactWrightFisher

makedocs(
    sitename = "ExactWrightFisher",
    format = :html,
    modules = [ExactWrightFisher]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
