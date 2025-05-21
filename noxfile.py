import nox

@nox.session
def lint(session):
    """Run linters."""
    session.install("flake8")
    session.run("flake8", "trm_py/subs/util.py")
    # Add coverage to linting #8

@nox.session
def tests(session):
    """Run tests."""
    pass
