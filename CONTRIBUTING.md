# Contributing to AnnoQC

Thank you for your interest in contributing to AnnoQC! This document provides guidelines for contributing to the project.

## Code of Conduct

This project and everyone participating in it is governed by our Code of Conduct. By participating, you are expected to uphold this code. Please report unacceptable behavior to [joseph.guhlin@gmail.com](mailto:joseph.guhlin@gmail.com).

## How to Contribute

### Reporting Bugs

Before creating bug reports, please check the existing issues to avoid duplicates. When creating a bug report, include:

- **Title**: Clear and descriptive
- **Description**: What happened and what you expected to happen
- **Steps to reproduce**: Minimal reproduction case
- **Environment**:
  - OS and version
  - AnnoQC version (`annoqc --version`)
  - Rust version (`rustc --version`)
  - External tool versions (DIAMOND, HMMER, MAFFT)
- **Logs**: Relevant error messages or output
- **Files**: Minimal input files to reproduce (if applicable)

### Suggesting Enhancements

Enhancement suggestions are welcome! Please provide:

- **Use case**: What problem would this solve?
- **Proposed solution**: How should it work?
- **Alternatives**: What other approaches did you consider?
- **Impact**: Who would benefit and how?

### Pull Requests

#### Development Setup

1. **Fork and clone**:
   ```bash
   git clone https://github.com/YOUR_USERNAME/AnnoQC.git
   cd AnnoQC
   ```

2. **Install pre-commit hooks** (recommended):
   ```bash
   pip install pre-commit
   ./scripts/install_hooks.sh
   ```

3. **Install development dependencies**:
   ```bash
   # Via pixi (recommended)
   pixi install

   # Or via cargo
   cargo build
   ```

4. **Run tests**:
   ```bash
   cargo test --all-features
   ```

#### Making Changes

1. **Create a branch**:
   ```bash
   git checkout -b feature/your-feature-name
   # or
   git checkout -b fix/your-bug-fix
   ```

2. **Make your changes** following the coding standards below

3. **Test your changes**:
   ```bash
   # Run all tests
   cargo test --all-features

   # Run clippy
   cargo clippy --all-targets -- -D warnings

   # Check formatting
   cargo fmt --all -- --check
   ```

4. **Update documentation** if needed:
   - Add/update doc comments
   - Update README.md if user-facing
   - Update relevant mdBook chapters

5. **Commit your changes**:
   ```bash
   git add .
   git commit -m "Brief description of changes"
   ```

   Commit message format:
   - Start with a verb (Add, Fix, Remove, Update, Refactor)
   - Be specific and concise
   - Reference issues: `Fixes #123` or `Refs #456`

6. **Push and create PR**:
   ```bash
   git push origin feature/your-feature-name
   ```

   Then create a pull request on GitHub.

#### PR Guidelines

- **Keep it focused**: One PR per feature/fix
- **Small PRs**: Easier to review and merge
- **Include tests**: All new code should have tests
- **Update docs**: Update documentation for user-facing changes
- **Pass CI**: All checks must pass (tests, clippy, formatting)

### Coding Standards

#### Rust Style

- Follow standard Rust style guidelines
- Use `cargo fmt` for formatting (enforced by pre-commit hooks)
- Avoid `unwrap()` in production code; use proper error handling
- Prefer `?` operator for error propagation
- Add doc comments for public API

#### Documentation

```rust
/// Brief description of what this does.
///
/// More detailed explanation if needed.
///
/// # Arguments
///
/// * `arg1` - Description of argument
/// * `arg2` - Description of argument
///
/// # Returns
///
/// Description of return value
///
/// # Examples
///
/// ```
/// use annoqc::function_name;
///
/// let result = function_name(arg1, arg2);
/// assert_eq!(result, expected);
/// ```
pub fn function_name(arg1: Type1, arg2: Type2) -> ReturnType {
    // implementation
}
```

#### Testing

- Write tests for new functionality
- Use descriptive test names
- Test edge cases and error conditions
- Use property-based tests (proptest) for invariants

```rust
#[test]
fn test_edge_case_description() {
    // Arrange
    let input = /* test input */;

    // Act
    let result = function_under_test(input);

    // Assert
    assert_eq!(result, expected);
}
```

## Development Workflow

### Running Tests

```bash
# Run all tests
cargo test

# Run tests with output
cargo test -- --nocapture

# Run specific test
cargo test test_name

# Run tests in parallel
cargo test -- --test-threads=4
```

### Checking Code Quality

```bash
# Format code
cargo fmt

# Check formatting without making changes
cargo fmt --all -- --check

# Run linter
cargo clippy --all-targets -- -D warnings

# Build documentation
cargo doc --no-deps --open
```

### Benchmarking

```bash
# Run benchmarks
cargo bench

# Compare with baseline
cargo install critcmp
critcmp baseline new
```

## Project Structure

```
AnnoQC/
├── src/                # Source code
│   ├── main.rs        # CLI entry point and orchestration
│   ├── ecs.rs         # Bevy ECS systems
│   ├── diamond.rs     # DIAMOND integration
│   ├── scoring.rs     # Scoring algorithms
│   └── ...
├── tests/             # Integration tests
├── book/              # Documentation (mdBook)
├── examples/          # Example configurations
├── scripts/           # Utility scripts
└── Cargo.toml         # Project metadata
```

### Adding Features

1. **New module**: Add to `src/` with `mod.rs` if needed
2. **CLI options**: Add to `main.rs` using `clap`
3. **Config options**: Add to config structs and update `config.example.toml`
4. **Output fields**: Update JSONL/CSV output schemas
5. **Tests**: Add unit tests in module, integration tests in `tests/`
6. **Documentation**: Update relevant mdBook chapters

### External Tool Integration

When integrating new external tools:

1. Add version check to `preflight.rs`
2. Add to `pixi.toml` dependencies
3. Update INSTALLATION.md
4. Update run manifest (tool version tracking)
5. Add tests with mocked tool output

## Release Process

Maintainers follow this process for releases:

1. Update version in `Cargo.toml`
2. Update `CHANGELOG.md`
3. Create git tag: `git tag -a v0.2.0 -m "Release v0.2.0"`
4. Push tag: `git push origin v0.2.0`
5. CI builds multi-platform binaries
6. Create GitHub Release with artifacts
7. Deploy documentation

## Community Guidelines

### Communication

- **GitHub Issues**: For bugs, features, questions
- **Discussions**: For design discussions, RFCs
- **Email**: [joseph.guhlin@gmail.com](mailto:joseph.guhlin@gmail.com) for security issues

### Getting Help

1. Check existing documentation and issues
2. Search codebase for similar implementations
3. Ask questions in GitHub Discussions
4. Join community channels (if available)

### Review Process

- Maintainers review PRs as time allows
- Be patient; review may take time
- Address review feedback promptly
- Ask for clarification if needed

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

## Recognition

Contributors are recognized in:
- `CONTRIBUTORS.md` file
- Release notes for significant contributions
- GitHub contributor statistics

Thank you for contributing to AnnoQC!
