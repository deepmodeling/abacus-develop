# ABACUS Copilot Instructions

Before generating code or reviewing pull requests in this repository, follow
the ABACUS development baseline in:

- `AGENTS.md`
- `docs/developers_guide/agent_governance.md`

Treat `AGENTS.md` as the short entry point and
`docs/developers_guide/agent_governance.md` as the complete rule source. Review
the pull request diff before general style feedback, and separate blocking
diff-scoped issues from historical or advisory observations.

Use the ABACUS review finding format for actionable findings:

```markdown
Rule: <rule name>
Severity: error | warning | info
Location: <file:line>
Reason: <specific trigger>
Suggested action: <concrete next step>
Exception: allowed | not allowed | human approval required
```
