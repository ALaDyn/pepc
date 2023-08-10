Contributing to PEPC
====================

All contributions are welcome. PEPC was designed to be flexible and modular as
far as possible. Adding your own interactions should be straightforward.
Bugfixes or extensions to the existing codebase are also more than welcome
though we encourage you to first file an issue for those.

In general, please feel free to submit issues and enhancement requests. Since
development is usually done within projects or by students, please do not expect
an immediate response.

Any **issue**, **commit**, branch, or merge request should indicate the
frontend if it applies only to a single frontend, i.e. please add `[XXX]` to the
respective title of an issue, commit, or MR and `_XXX` to the branch name (where
XXX is the name of the frontend).

Contributing
------------

Please refer to PEPC's coding style when submitting patches and additions. In
general, we follow the "fork-and-pull" Git workflow.

  1. **Fork** the repository
  2. **Clone** the project to your own machine
  3. **Commit** changes to your _own branch_
  4. **Push** your work back up to your fork
  5. Submit a **Pull request** so that we can review your changes
     (ideally referencing issues)

NOTE: Be sure to merge the latest from "upstream" before making a pull request!

Run existing tests and add new tests for new features or encountered errors.

Contributor License Agreement
-----------------------------

Before we can fully accept and incorporate your changes, please fill in and
return our 'Contributor License Agreement' (CLA.pdf). This can be done fully
electronically (filling in the .pdf form with an approriate pdf viewer) adding
your electronic signature with a x.509 certificate. If signing electronically,
please make sure a chain of trust can be verified (i.e. use a verifiable,
official certificate).

In case of any questions, please do not hesitate to get in touch with us.

Coding Style
------------

The codebase conforms to a single, consistent coding style. New submissions will
be rejected if they do not conform to the current style.
   * Apart from MPI specifiers and pre-processor macros, PEPC generally
     makes use of lower case syntax.
   * Instead of `enddo` and `endif` including the additional space is
     preferred, i.e. `end do` and `end if`
   * Subroutines and functions should contain comments documenting their
     purpose and interface.
   * Expect the use of `doxygen` or `FORD` in the future, so write comments
     accordingly and 'guard' comments if necessary.
   * `git hooks` are in place to check the style of source files before any
     commit. Those checks are two-fold:

       1. a whitespace check via `git diff-index --check`
       2. running `fprettify` on changed files

     If either check fails, the commit will be aborted.

Comment your code if not self-describing.

Please also follow the standard commit-message format. The first line is the
subject, and should generally be less than 50 characters. The second line must
be blank. The subsequent lines are the message body, and should generally be
concise but describe the changes that have been made. Commit often and do not
mix different changes (e.g. fixes for more than one issue or changes that are of
a different nature) on one commit.

fprettify
---------

PEPC's source files are best formatted with the help of
[`fprettify`](https://github.com/pseewald/fprettify). Please make sure to use a
recent, bug fixed version to not introduce errors. You can check the
`Dockerfile` for the version used for CI. All settings for `fprettify`
are provided in a configuration file found under `src/.fprettify.rc`. To run a
check without changing anything, perform a
```sh
fprettify -c src/.fprettify.rc -d src/<path_to_file>/<file>
```
from PEPC's root directory. Dropping the `-d` will apply changes.

<!-- vim: set ts=4 sw=4 tw=80 et :-->
