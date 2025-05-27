# Contributing

Contributions are welcome here are some hints to get you started

## Devloping on branches

It is expected that people branch from develop before making changes

After using get the source code from the README:

> `git checkout -b <my-branch-name>`

or even better create an issue describing what you are looking to do create an issue branch and check this out!

If you are modifying the codes `src/cpp-<name>` then also create a branch for them.

Pushes to branches wont automatically trigger workflows check your code locally using the nox testing framework (in development) and the local build detailed in the README. Remember this code is portable across Linux and Mac so you should test on both those systems.
A PR/push to develop will trigger the Workflow to run CIbuildwheel and then push to testpypi if you have privilages (You will need to bump your version number). You should download the build artifacts (the wheels) and use the test-pypi to check your wheel. Test runners will (eventually) be made to handle this automatically.

PRs to the main branch trigger the ciworkflow again and push to PyPi which releases the new version of the package to the world make sure you are versioning correctly.

## Understanding submodules

This package uses submodules to creates separation between the various sub-packages.
See the 'Get the source code' section of README to download these packages.

if you want to edit a subpackage:

> `cd src/subpackage`
>
> `checkout (-b) <branchname>`
>
> make changes and create commits do pushes.
>
> >
> > While developing you can use the local build method from the top directory, but not cibuildwheel.
> >
>
> make a PR to main on the submodule (or merge if you have privilege)
>
> `checkout main`
>
> make sure all is up to date and as expected
>
> `cd ../..` (back to the main repo trm-py)
>
> Check you are on a branch that reflects the changes you have made
>
> `git status` should show new commits in the submodule
>
> `git add src/subpackage`
>
> `git commit -m "useful commit message"`
>
> `git push origin`
>
> The runners should do building and testing for you, however local building and testing via the local build method is encouraged.