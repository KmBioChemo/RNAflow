# Security Policy

## Reporting a vulnerability

If you discover a security issue in RNAflow, please report it
**privately** rather than opening a public issue. Email the maintainer
at **<karim.matmat@unibas.ch>** with a description and, where possible,
steps to reproduce. You can expect an acknowledgement within a few
working days.

## Things to be aware of

RNAflow is an interactive analysis application meant to run locally, or
on a server you control. In particular:

- **AI-assisted interpretation is opt-in and sends data to a third
  party.** When you enable it and provide an API key, RNAflow sends gene
  names, fold-changes, FDRs and enrichment terms — **but not the raw
  count matrix** — to the Anthropic Claude API. Do **not** enable it for
  confidential or unpublished data that you cannot share with an
  external API, and be aware that on a shared/remote deployment this
  data leaves the host.
- **Project files execute
  [`readRDS()`](https://rdrr.io/r/base/readRDS.html).** Opening a
  `.rnaflow.rds` project deserializes an R object; only open project
  files from sources you trust.
- **API keys** are read from the `ANTHROPIC_API_KEY` environment
  variable or entered in the session and are never written to saved
  projects or reports.

## Supported versions

Security fixes are applied to the latest release on the `main` branch.
Please update to the latest version before reporting an issue.
