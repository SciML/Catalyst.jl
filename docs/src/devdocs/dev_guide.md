# Catalyst Developer Documentation

## Release Process
Beginning with v15, Catalyst is using a new release process to try to ensure
continuing stability of releases. Before making a release one should

1. Create a new release branch, i.e. "release-15.0.0"
2. On this branch, cap major dependencies to their latest version that works and
   for which tests pass.
   - Caps need to be included in both Project.toml and docs/Project.toml.
   - Do not cap the master branch as this can prevent upstream libraries from
     properly testing against Catalyst, and hide breaking changes that impact
     Catalyst.
3. Check docs build with the capped dependencies. Visually verify via checking
   the artifact in the doc build that the docs actually look ok (since sometimes
   issues can arise that do not lead to actual errors in the doc CI).
4. Release via the [registration
   issue](https://github.com/SciML/JumpProcesses.jl/issues/73) with the
   command:
   
    ```
    @JuliaRegistrator register branch=release-15.0.0
    ```
    
    modifying as appropriate for the version you are releasing.

If there is subsequently a need to increment the version of a dependency, this
should be done via a new release that follows the above process, and modifies
the [patch, minor, or major Catalyst version (as appropriate for the potential
impact of the dependency change on Catalyst users)](https://semver.org/). If the
dependency being updated is a non-breaking release, and would have automatically
been installed by the package resolver had it not been capped, a patch release
should be preferred. If the new release branch is branched from master, *it
needs to ensure Project.toml caps are all ≥ to those listed in the previous
Catalyst release branch*.
