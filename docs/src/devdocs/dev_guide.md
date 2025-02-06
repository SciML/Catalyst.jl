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
   issue](https://github.com/SciML/Catalyst.jl/issues/127) with the
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

## Development advice

### Inspecting documentation of a PR
When updating documentation it is typically useful to view how the updated documentation will look like before completing a PR. This can be done both locally and directly through GitHub. Here, some issues are very easy to see in the actually built html files, but hard to see in the text files from which these are generated.

Whenever a PR to Catalyst is created, CI will create a corresponding documenter build job. If this one passes, you can access the built documentation (which will be the new Catalyst documentation) from it. Follow these steps to view the built docs:
1. Click on the build job (at the bottom of the PR conversation). 
2. Expand the "Upload site as artifact" section (in the large black are).
3. Click on the links following the "Artifact download URL: " text.
4. This will download a zip file containing the documentations. Extract it to a location on your computer and then open the "index.html" file.

To build the Catalyst documentation locally:
1. Run the ".julia/dev/Catalyst/docs/make.jl" file. ALternatively, open a Julia session, activate the "docs" environment, and run the file using `include("make.jl").
2. Open the ".julia/dev/Catalyst/docs/build/index.html" file.


### Spellchecking in your code
Especially when writing documentation, but also in other situation, it can be useful to have a spellchecker running through your code. While code can be copied into a spellchecker and checked there, it can also be very useful to (for users of VSCode) run the [Code Spell Checker](https://marketplace.visualstudio.com/items?itemName=streetsidesoftware.code-spell-checker) extension, which will automatically provide simple spell checks for code and documentation as you write it.


