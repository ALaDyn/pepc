Contributing to PEPC
====================

All contributions are welcome. PEPC was designed to be flexible and modular as
far as possible. Adding your own interactions should be straightforward.
Bugfixes or extensions to the existing codebase are also more than welcome
though we encourage you to first file an issue for those.

In general, please feel free to submit issues and enhancement requests. Since
development is usually done within projects or by students, please do not expect
an immediate response.

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

The codebase conforms to a single, consistent coding style (though not formally
defined yet). New submissions will be rejected if they do not conform to the
current style.
   * Apart from MPI specifiers and pre-processor macros, PEPC generally
     makes use of lower case syntax.
   * Instead of `enddo` and `endif` including the additional space is
     preferred, i.e. `end do` and `end if`
   * Subroutines and functions should contain comments documenting their
     purpose and interface.
   * Expect the use of `fprettyfy` and `doxygen` or `FORD` in the future, so
     write comments accordingly and 'guard' comments if necessary.

Comment your code if not self-describing.

Run existing tests and add new tests for new features or encountered errors.

Please also follow the standard commit-message format. The first line is the
subject, and should generally be less than 50 characters. The second line must
be blank. The subsequent lines are the message body, and should generally be
concise but describe the changes that have been made. Commit often and do not
mix different changes (e.g. fixes for more than one issue or changes that are of
a different nature) on one commit.
