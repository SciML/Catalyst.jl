# Catalyst Developer Documentation

## Release Process
Beginning with v15, Catalyst is using a new release process to try to ensure
continuing stability of releases. Before making a release one should

1. Create a new release branch, i.e. "release-15.0.0".
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

### Checking doc builds for errors
When updating documentation, Catalyst will run any Julia code provided within example blocks to dynamically create figures and outputs. In addition to automatically creating these for us, it also provides an automatic check that all code in documentation is correct. Here, if any of the documentation code throws an error, the build job will fail. The documentation build job can be found at the bottom of a PRs conversation, here is an example of a failed one:

![Failed builddocs link](../assets/devdocs/failed_builddocs_link.png)

To check what errors were produced, click on the "Details" link of the job. Next, any errors can be found at the bottom of the "Build and deploy" section (which should be opened automatically).

### Inspecting documentation of a PR or branch
When updating documentation it is typically useful to view the updated documentation in HTML format (which is the format users will see). Here, some errors are much easier to spot in html format as compared with the raw text files from which these are generated. There are two primary ways to view updated documentation, either by downloading them from the PR or by building the docs locally.

Whenever a PR to Catalyst is created, CI will create a corresponding documenter build job. If the build job passes, you can access the built documentation (which will be the new Catalyst documentation if the PR is merged) through the following steps:
1. Click on "Details" in the documentation build job (at the bottom of the PR conversation tab). 
2. Expand the "Upload site as artifact" section.
3. Click on the link at the end (which follows the "Artifact download URL: " text).
4. This will download a zip folder containing the documentation. Extract it to a location on your computer and then open the "index.html" file.

To build the Catalyst documentation locally:
1. Navigate to the ".julia/dev/Catalyst/docs/" folder and run the "make.jl" file using ">julia make.jl". Alternatively, open a Julia session, activate the "docs" environment, and run the file using `include("make.jl").
2. Open the ".julia/dev/Catalyst/docs/build/index.html" file.

### Spellchecking in your code
Especially when writing documentation, but also when writing normal code, it can be useful to have a spellchecker running through your texts. While code can be copied into a spellchecker and checked there (which is still useful to check grammar), it can also be very useful to (for users of VSCode) run the [Code Spell Checker](https://marketplace.visualstudio.com/items?itemName=streetsidesoftware.code-spell-checker) extension. This will automatically provide simple spell-checking for code and documentation as you write it.
