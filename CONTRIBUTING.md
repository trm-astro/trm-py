# Contributing

Contributions are welcome here are some hints to get you started

## Devloping on branches



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